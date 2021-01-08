import pandas
import itertools
import time
import sys
import json
from collections import defaultdict
import numpy as np
import os
import argparse
import logging
import allel
from multiprocessing import Process, Queue, cpu_count


logger = logging.getLogger()
defaultLogLevel = "DEBUG"

# you must define the PYENSEMBL_CACHE_DIR before importing the pyensembl module
logger.info('setting pyensembl dir to /var/tmp/pyensembl-cache')
os.environ['PYENSEMBL_CACHE_DIR'] = '/var/tmp/pyensembl-cache'
import pyensembl

# p2 = P(VUS is pathogenic and patient carries a pathogenic variant in trans) (arbitrarily set by goldgar et al)
# Integrated Evaluation of DNA Sequence Variants of Unknown Clinical Significance: Application to BRCA1 and BRCA2
brca1_p2 = 0.0001
brca2_p2 = 0.001

classStrings = { 'Pathogenic':[ 'Pathogenic' ], 'Benign':[ 'Benign', 'Likely benign' ],
                 'Unknown': [ 'Uncertain significance', '-']}
sigColName = 'Clinical_significance_ENIGMA'
coordinateColumnBase = 'Genomic_Coordinate_hg'
alleleFrequencyName = 'Allele_frequency_ExAC'


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, set):
            return list(obj)
        elif isinstance(obj, tuple):
            return list(obj)
        else:
            return super(NpEncoder, self).default(obj)

class NpDecoder(json.JSONDecoder):
    def default(self, obj):
        if isinstance(obj, list):
            return set(obj)
        else:
            return super(NpDecoder, self).default(obj)

def parseArgs():
    parser = argparse.ArgumentParser(usage="cooccurrenceFinder args [options]")
    parser.add_argument("--vcf", dest="vcf", help="name of file containing VCF data, default=None", default=None)
    parser.add_argument("--ipv", dest="ipv", help="ipv file name, default=ipv.json", default='ipv.json')
    parser.add_argument("--vpi", dest="vpi", help="vpi file name, default=vpi.json", default='vpi.json')
    parser.add_argument("--out", dest="out", help="output file name, default=out.json", default='out.json')
    parser.add_argument("--all", dest="all", help="all vars file name, default=all.json", default='all.json')
    parser.add_argument("--anno", dest="anno", help="annotation file name, default=None", default=None)
    parser.add_argument("--h", dest="h", help="Human genome version (37 or 38). Default=None", default=None)
    parser.add_argument("--e", dest="e", help="Ensembl version - 75 (for 37) or 99 (for 38). Default=None", default=None)
    parser.add_argument("--c", dest="c", help="Chromosome of interest. Default=None", default=None)
    parser.add_argument("--g", dest="g", help="Gene of interest. Default=None", default=None)
    parser.add_argument("--p", dest="p", help="Phased (boolean). Default=False", default='True')
    parser.add_argument("--s", dest="s", help="Save variants per individual. Default=True", default='True')
    parser.add_argument("--i", dest="i", help="Include pathog vars per VUS in report. Default=True", default='True')
    parser.add_argument("--n", dest="n", help="Number of processes. Default=1", default=cpu_count())
    parser.add_argument("--b", dest="b", help="BRCA variants file. Default=brca-variants", default=None)
    parser.add_argument("--d", dest="d", help="directory containing pyensembl-cache. Default=/var/tmp/pyensembl-cache", default='/var/tmp/pyensembl-cache')
    parser.add_argument("--r", dest="r", help="Rare frequency cutoff. Default=0.01", default=0.01)
    parser.add_argument("--f", dest="f", help="Topmed freeze. Default=0", default=0)
    parser.add_argument("--log", dest="logLevel", help="Logging level. Default=%s" % defaultLogLevel, default=defaultLogLevel)
    return parser.parse_args()

def configureLogging(options):
    # Parse the log level
    numeric_level = getattr(logging, options.logLevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.logLevel)
    # Setup a logger
    logger.setLevel(numeric_level)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(numeric_level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.debug("Established logger")

def main():
    # main just parses the CLI and then calls the run() method with appropriate args
    options = parseArgs()

    configureLogging(options)

    print(options)


    run(int(options.h), int(options.e), options.c, options.g, bool(eval(options.p)), options.vcf, bool(eval(options.s)),
        int(options.n), options.b, options.d, float(options.r), options.ipv, options.vpi, options.all, options.anno,
        options.out, int(options.f))

def run(hgVersion, ensemblRelease, chromosome, gene, phased, vcfFileName, saveVarsPerIndivid, numProcs,
        brcaFileName, pyensemblDir, rareCutoff, ipvFileName, vpiFileName, allVariantsFileName, annoFileName,
        outputFileName, freeze):


    logger.info('setting pyensembl dir to ' + pyensemblDir)
    os.environ['PYENSEMBL_CACHE_DIR'] = '/var/tmp/pyensembl-cache'

    logger.info('reading annotation data from ' + annoFileName)
    with open(annoFileName, 'r') as f:
        annoDF = pandas.read_csv(annoFileName, header=0, sep='\t')
    f.close()

    logger.info('reading BRCA data from ' + brcaFileName)
    t = time.time()
    brcaDF, pathogenicVariants, benignVariants, unknownVariants = findVariantsInBRCA(brcaFileName, classStrings, hgVersion)
    logger.info('elapsed time in findVariantsInBRCA() ' + str(time.time() -t))

    logger.info('number of pathogenic variants is ' + str(len(pathogenicVariants)))
    logger.info('number of benign variants is ' + str(len(benignVariants)))
    logger.info('number of vus variants is ' + str(len(unknownVariants)))

    logger.info('reading VCF file ' + vcfFileName)
    t = time.time()
    vcf = readVCFFile(vcfFileName)
    logger.info('elapsed time in readVCFFile() ' + str(time.time() -t))

    logger.info('finding variants per individual in ' + vcfFileName)
    t = time.time()
    q = Queue()
    w = Queue()
    processList = list()
    for i in range(numProcs):
        p = Process(target=findVarsPerIndividual, args=(q, w, vcf, benignVariants, pathogenicVariants, chromosome, gene,
                                                        ensemblRelease, annoDF, i, numProcs, freeze, ))
        p.start()
        processList.append(p)
    logger.info('joining results from forked threads')
    variantsPerIndividual = dict()
    for i in range(numProcs):
        variantsPerIndividual.update(q.get())
    for i in range(numProcs):
        processList[i].join()
    logger.info('elapsed time in findVariantsPerIndividual() ' + str(time.time() -t))

    t = time.time()
    cohortSize = len(variantsPerIndividual)
    logger.info('number of samples is ' + str(cohortSize))
    individualsPerVariant = findIndividualsPerVariant(variantsPerIndividual, vcf, chromosome,brcaDF, hgVersion,
                                                        ensemblRelease, cohortSize)
    logger.info('number of records is ' + str(len(individualsPerVariant)))
    logger.info('elapsed time in updateIndividualsPerVariant() ' + str(time.time() -t))

    logger.info('saving vpi to ' + vpiFileName)
    with open(vpiFileName, 'w') as f:
        json.dump(variantsPerIndividual, f, cls=NpEncoder)
    f.close()

    logger.info('saving ipv to ' + ipvFileName)
    with open(ipvFileName, 'w') as f:
        json.dump(individualsPerVariant, f, cls=NpEncoder)
    f.close()

    logger.info('finding homozygous individuals per vus')
    t = time.time()
    homozygousPerVus = countHomozygousPerVus(variantsPerIndividual, brcaDF, hgVersion, ensemblRelease, gene, rareCutoff)
    logger.info('elapsed time in countHomozygousPerVus() ' + str(time.time() -t))

    #logger.info('finding homozygous individuals per benign')
    #t = time.time()
    #homozygousPerBenign = countHomozygousPerBenign(variantsPerIndividual, brcaDF, hgVersion, ensemblRelease, gene, rareCutoff)
    #logger.info('elapsed time in countHomozygousPerBenign() ' + str(time.time() - t))

    logger.info('finding individuals per cooc')
    t = time.time()
    individualsPerPathogenicCooccurrence, n, k = findIndividualsPerCooccurrence(variantsPerIndividual, ensemblRelease,
                                                                                phased, gene)
    logger.info('elapsed time in findIndividualsPerCooccurrence() ' + str(time.time() -t))


    # TODO check this math!
    logger.info('calculating p1')
    # p1 = P(VUS is benign and patient carries a path variant in trans) = 0.5 * overall freq of path muts in cohort
    # calculate total number of benign, pathogenic, and vus variants in cohort
    logger.info('getting all variants for cohort')
    allVariants = getAllVariantsPerClass(variantsPerIndividual)
    numPathogenic = len(set(allVariants['pathogenic']))
    p1 =  0.5 * numPathogenic / cohortSize

    logger.info('saving all variants to ' + allVariantsFileName)
    json_dump = json.dumps(allVariants, cls=NpEncoder)
    with open(allVariantsFileName, 'w') as f:
        f.write(json_dump)
    f.close()

    logger.info('putting all the data together per vus')
    dataPerVus = calculateLikelihood(individualsPerPathogenicCooccurrence, p1, n, k, brcaDF, hgVersion, cohortSize)

    #data_set = {"cooccurring vus": dataPerVus, "homozygous vus": homozygousPerVus, "homozygous benign": homozygousPerBenign}
    data_set = {"cooccurring vus": dataPerVus, "homozygous vus": homozygousPerVus}
    json_dump = json.dumps(data_set, cls=NpEncoder)

    logger.info('saving final VUS data  to ' + outputFileName)
    with open(outputFileName, 'w') as f:
        f.write(json_dump)
    f.close()

def findIndividualsPerVariant(variantsPerIndividual, vcf, chromosome, brcaDF, hgVersion, ensemblRelease, cohortSize):
    individualsPerVariant = dict()
    for individual in variantsPerIndividual:
        for b in variantsPerIndividual[individual]['benign']:
            v = str((b[0][0], b[0][1], b[0][2], b[0][3]))
            if not v in individualsPerVariant:
                individualsPerVariant[v] = {'heterozygous individuals': set(),
                                            'homozygous individuals': set()}
            if b[1] == '1' or b[1] == '2':
                individualsPerVariant[v]['heterozygous individuals'].add(individual)
            elif b[1] == '3':
                individualsPerVariant[v]['homozygous individuals'].add(individual)
            else:
                logger.warning('hmm - didnt add this ben ' + v)
        for p in variantsPerIndividual[individual]['pathogenic']:
            v = str((p[0][0], p[0][1], p[0][2], p[0][3]))
            if not v in individualsPerVariant:
                individualsPerVariant[v] = {'heterozygous individuals': set(),
                                            'homozygous individuals': set()}
            if p[1] == '1' or p[1] == '2':
                individualsPerVariant[v]['heterozygous individuals'].add(individual)
            elif p[1] == '3':
                individualsPerVariant[v]['homozygous individuals'].add(individual)
            else:
                logger.warning('hmm - didnt add this path ' + v)
        for vus in variantsPerIndividual[individual]['vus']:
            v = str((vus[0][0], vus[0][1], vus[0][2], vus[0][3]))
            if not v in individualsPerVariant:
                individualsPerVariant[v] = {'heterozygous individuals': set(),
                                            'homozygous individuals': set()}
            if vus[1] == '1' or vus[1] == '2':
                individualsPerVariant[v]['heterozygous individuals'].add(individual)
            elif vus[1] == '3':
                individualsPerVariant[v]['homozygous individuals'].add(individual)
            else:
                logger.warning('hmm - didnt add this vus ' + v)
    individualsPerVariant = addVariantInfo(individualsPerVariant, vcf, chromosome, ['FIBC_I', 'FIBC_P'], brcaDF,
                                           hgVersion, cohortSize, ensemblRelease)

    return individualsPerVariant


def isExonic(ensemblRelease, chrom, pos):
    ensembl = pyensembl.EnsemblRelease(release=ensemblRelease)
    try:
        exons = ensembl.exons_at_locus(contig=int(chrom), position=int(pos))
    except Exception as e:
        logger.error('exception: ' + str(e))
        return None
    return len(exons) > 0

def addVariantInfo(individualsPerVariant, vcf, chromosome, infoList, brcaDF, hgVersion, cohortSize, ensemblRelease):
    # add infoList stuff from INFO field
    for variant in range(len(vcf['calldata/GT'])):
        if int(vcf['variants/CHROM'][variant].replace('chr', '')) != int(chromosome):
            continue
        c = int(vcf['variants/CHROM'][variant].replace('chr', ''))
        p = int(vcf['variants/POS'][variant])
        r = str(vcf['variants/REF'][variant])
        a = str(vcf['variants/ALT'][variant][0])
        v = str((c,p,r,a))
        if v in individualsPerVariant:
            '''for info in infoList:
                individualsPerVariant[v][info] = vcf['variants/' + info][variant]'''
            maxPop, maxFreq, minPop, minFreq = getGnomadData(brcaDF, eval(v), hgVersion)
            individualsPerVariant[v]['maxPop'] = maxPop
            individualsPerVariant[v]['maxFreq'] = maxFreq
            individualsPerVariant[v]['minPop'] = minPop
            individualsPerVariant[v]['minFreq'] = minFreq
            individualsPerVariant[v]['cohortFreq'] = float(len(individualsPerVariant[v]['homozygous individuals']) + \
                len(individualsPerVariant[v]['heterozygous individuals']) ) / float(cohortSize)
            individualsPerVariant[v]['exonic'] = isExonic(ensemblRelease, c, p)
        else:
            logger.warning('variant ' + str(v) + ' not in ipv dict?')
    return individualsPerVariant

def getAllVariantsPerClass(vpi):
    allVariants = dict()
    allVariants['benign'] = list()
    allVariants['pathogenic'] = list()
    allVariants['vus'] = list()

    for i in vpi:
        for b in vpi[i]['benign']:
            allVariants['benign'].append(b)
        for p in vpi[i]['pathogenic']:
            allVariants['pathogenic'].append(p)
        for v in vpi[i]['vus']:
            allVariants['vus'].append(v)

    return allVariants

def findVariantsInBRCA(fileName, classStrings, hgVersion):
    brcaDF = pandas.read_csv(fileName, sep='\t', header=0, dtype=str)
    # Genomic_Coordinate_hg37
    # chr13:g.32972575:G>T
    pathVars = set()
    benignVars = set()
    vusVars = set()
    for i in range(len(brcaDF)):
        # TODO use HGVS? problem is indel representation
        # TODO VCF has standard of using left-most pos where HGVS has standard of using right-most pos for indel
        # if cDNA (3' side or 5'?) => standard is using 5' strand
        coord = brcaDF.loc[i, coordinateColumnBase + str(hgVersion)].split(':')
        chrom = int(coord[0].split('chr')[1])
        pos = int(coord[1].split('g.')[1])
        ref, alt = coord[2].split('>')
        var = (chrom, pos, ref, alt)
        if str(brcaDF.loc[i, sigColName]) in classStrings['Pathogenic']:
            pathVars.add(var)
        elif str(brcaDF.loc[i, sigColName]) in classStrings['Benign']:
            benignVars.add(var)
        elif str(brcaDF.loc[i, sigColName]) in classStrings['Unknown']:
            vusVars.add(var)
        else:
            continue

    return brcaDF, pathVars, benignVars, vusVars

def getGnomadData(brcaDF, vus, hgVersion):
    # TODO write a unit test
    # 13:g.32393468:C>CT
    #hgString = 'chr' + str(vus[0][0]) + ':g.' + str(vus[0][1]) + ':' + str(vus[0][2]) + '>' + str(vus[0][3])
    hgString = 'chr' + str(vus[0]) + ':g.' + str(vus[1]) + ':' + str(vus[2]) + '>' + str(vus[3])

    # first, get list of columns for GnomAD allleles
    gnomad = [v for v in list(brcaDF.columns) if 'GnomAD' in v]
    alleleFrequencies = [v for v in gnomad if 'Allele_frequency' in v]

    # second, get frequencies across exomes and genomes to determine max
    # return population, frequency, count, and number
    # replace "frequency" with "count" and "number" in Allele_frequency_genome_AFR_GnomAD

    maxFrequency = 0.0
    maxPopulation = None
    minFrequency = 1.0
    minPopulation = None
    for af in alleleFrequencies:
        freq=0.0
        alleleFreqList = brcaDF[brcaDF[coordinateColumnBase + str(hgVersion)] == hgString][af].tolist()
        if alleleFreqList:
            try:
                freq = float(alleleFreqList[0])
            except ValueError:
                continue
            if freq > maxFrequency:
                maxFrequency = freq
                maxPopulation = af
            if freq < minFrequency:
                minFrequency = freq
                minPopulation = af

    return (maxPopulation, maxFrequency, minPopulation, minFrequency)

def countHomozygousPerBenign(variantsPerIndividual, brcaDF, hgVersion, ensemblRelease, geneOfInterest, rareCutoff):
    homozygousPerBenign = dict()

    for individual in variantsPerIndividual:
        for ben in variantsPerIndividual[individual]['benign']:
            if (ben[1] == '3') and (getGenesForVariant(ben[0], ensemblRelease, geneOfInterest)):
                if str(ben[0]) not in homozygousPerBenign:
                    homozygousPerBenign[str(ben[0])] = dict()
                    homozygousPerBenign[str(ben[0])]['count'] = 0
                    maxPop, maxPopFreq, minPop, minPopFreq = getGnomadData(brcaDF, ben[0], hgVersion)
                    homozygousPerBenign[str(ben[0])]['maxPop'] = maxPop
                    homozygousPerBenign[str(ben[0])]['maxPopFreq'] = maxPopFreq
                    homozygousPerBenign[str(ben[0])]['minPop'] = minPop
                    homozygousPerBenign[str(ben[0])]['minPopFreq'] = minPopFreq
                homozygousPerBenign[str(ben[0])]['count'] += 1

    cohortSize = len(variantsPerIndividual)
    for ben in homozygousPerBenign:
        maxPopFreq = homozygousPerBenign[ben]['maxPopFreq']
        cohortFreq = float(homozygousPerBenign[ben]['count'])/ float(cohortSize)
        homozygousPerBenign[ben]['cohortFreq'] = float(cohortFreq)
        if cohortFreq < rareCutoff or (maxPopFreq != 0 and maxPopFreq < rareCutoff):
            homozygousPerBenign[ben]['RARE']= True
        else:
            homozygousPerBenign[ben]['RARE'] = False


    return homozygousPerBenign


def countHomozygousPerVus(variantsPerIndividual, brcaDF, hgVersion, ensemblRelease, geneOfInterest, rareCutoff):
    homozygousPerVus = dict()

    for individual in variantsPerIndividual:
        for vus in variantsPerIndividual[individual]['vus']:
            if (vus[1] == '3') and (getGenesForVariant(vus[0], ensemblRelease, geneOfInterest)):
                if str(vus[0]) not in homozygousPerVus:
                    homozygousPerVus[str(vus[0])] = dict()
                    homozygousPerVus[str(vus[0])]['count'] = 0
                    maxPop, maxPopFreq, minPop, minPopFreq = getGnomadData(brcaDF, vus[0], hgVersion)
                    homozygousPerVus[str(vus[0])]['maxPop'] = maxPop
                    homozygousPerVus[str(vus[0])]['maxPopFreq'] = maxPopFreq
                    homozygousPerVus[str(vus[0])]['minPop'] = minPop
                    homozygousPerVus[str(vus[0])]['minPopFreq'] = minPopFreq
                homozygousPerVus[str(vus[0])]['count'] += 1

    cohortSize = len(variantsPerIndividual)
    for vus in homozygousPerVus:
        #maxPopFreq = homoZygousPerVus[vus][1][1]
        maxPopFreq = homozygousPerVus[vus]['maxPopFreq']
        cohortFreq = float(homozygousPerVus[vus]['count'])/ float(cohortSize)
        homozygousPerVus[vus]['cohortFreq'] = float(cohortFreq)
        if cohortFreq < rareCutoff or (maxPopFreq != 0 and maxPopFreq < rareCutoff):
            homozygousPerVus[vus]['RARE']= True
        else:
            homozygousPerVus[vus]['RARE'] = False


    return homozygousPerVus

def calculateLikelihood(pathCoocs, p1, n, k, brcaDF, hgVersion, cohortSize):

    # vus coocs data: {(vus1, vus2):[individuals]}
    # "([10, 89624243, 'A', 'G'], [10, 89624304, 'C', 'T')]": ["0000057940", "0000057950"],
    # "([10, 89624304, 'C', 'T'], [10, 89624243, 'A', 'G')]": ["0000057940", "0000057950"],
    # "([10, 89624243, 'A', 'G'], [10, 89624298, 'C', 'T')]": ["0000057950"],
    # "([10, 89624304, 'C', 'T'], [10, 89624298, 'C', 'T')]": ["0000057950"],
    # "([10, 89624298, 'C', 'T'], [10, 89624243, 'A', 'G')]": ["0000057950"],
    # "([10, 89624298, 'C', 'T'], [10, 89624304, 'C', 'T')]": ["0000057950"]}

    # new version: {(vus, path):[individuals]}
    # {
    # ((13, 32906579, 'A', 'C'), (13, 32913597, 'CAGAA', 'C')): ['0000061858', '000061999'],
    # ((13, 32906579, 'A', 'C'), (13, 32907420, 'GA', 'G')): ['0999937461']
    # }

    # now calculate log likelihood ratios!
    likelihoodRatios = dict()
    for vus in n:
        # decide whether to use brca1_p2 or brca2_p2
        if vus[0] == 13:
            p2 = brca2_p2
        elif vus[0] == 17:
            p2 = brca1_p2
        else:
            logger.error("unknown chromosome: " + str(vus[0]))
            continue

        n_ = n[vus]
        k_ = k[vus]
        denom = ((p1 ** k_) * (1 - p1) ** (n_ - k_))
        if k_ == 0:
            continue
        elif denom == 0:
            likelihoodRatios[vus] = sys.float_info.min
        else:
            likelihoodRatios[vus] = ((p2 ** k_) * (1 - p2) ** (n_ - k_)) / ((p1 ** k_) * (1 - p1) ** (n_ - k_))


    # find all the pathogenic variants this vus co-occurred with
    pathVarsPerVus = defaultdict(list)
    for cooc in pathCoocs:
        vus = cooc[0]
        pathVarsPerVus[vus].append(cooc[1])

    # put it all together in a single dict
    dataPerVus = dict()
    for vus in likelihoodRatios:
        if vus[0] == 13:
            p2 = brca2_p2
        elif vus[0] == 17:
            p2 = brca1_p2
        else:
            logger.error("unknown chromosome: " + str(vus[0]))
            continue
        maxPop, maxPopFreq, minPop, minPopFreq = getGnomadData(brcaDF, vus, hgVersion)
        cohortFreq = float(n[vus]) / float(cohortSize)
        data = {'likelihood data': {'p1':p1, 'p2':p2, 'n':n[vus], 'k':k[vus], 'likelihood':likelihoodRatios[vus]},
                    'allele frequencies':{'maxPop':maxPop, 'maxPopFreq':maxPopFreq, 'minPop': minPop, 'minPopFreq': minPopFreq,
                                          'cohortFreq':cohortFreq}, 'pathogenic variants': pathVarsPerVus[vus]}


        dataPerVus[str(vus)] = data

    return dataPerVus

def readVCFFile(vcfFileName):
    return allel.read_vcf(vcfFileName, fields=['samples', 'calldata/GT', 'variants/ALT', 'variants/CHROM',
                                               'variants/FILTER_PASS', 'variants/ID', 'variants/POS', 'variants/QUAL',
                                               'variants/REF', 'variants/INFO'])

def divide(n, d):
   res = list()
   qu = int(n/d)
   rm = n%d
   for i in range(d):
       if i < rm:
           res.append(qu + 1)
       else:
           res.append(qu)
   return res

def getStartAndEnd(partitionSizes, threadID):
    start = 0
    for i in range(threadID):
        start += partitionSizes[i]
    end = start + partitionSizes[threadID]

    return start, end

def findVarsPerIndividual(q, w, vcf, benignVariants, pathogenicVariants, chromosome, gene, ensemblRelease, annoDF,
                          threadID, numProcesses, freeze):
    '''infoFields = ['variants/ABE', 'variants/ABZ', 'variants/AC', 'variants/AF',
     'variants/AN', 'variants/ANN', 'variants/AVGDP', 'variants/BETA_IF', 'variants/BQZ',
     'variants/CYZ', 'variants/FIBC_I', 'variants/FIBC_P', 'variants/FILTER_PASS', 'variants/FLT20', 'variants/GC',
     'variants/GN', 'variants/HWEAF_P', 'variants/HWE_SLP_I', 'variants/HWE_SLP_P', 'variants/IBC_I', 'variants/IBC_P',
     'variants/ID', 'variants/IOR', 'variants/MAX_IF', 'variants/MIN_IF', 'variants/NM0', 'variants/NM1',
     'variants/NMZ', 'variants/NS_NREF', 'variants/QUAL', 'variants/STZ', 'variants/SVM']'''

    variantsPerIndividual = dict()
    individuals = list(vcf['samples'])
    n = len(individuals)
    partitionSizes = divide(n, numProcesses)
    start,end = getStartAndEnd(partitionSizes, threadID)
    logger.info('threadID = ' + str(threadID) + ' processing from ' + str(start) + ' to ' + str(end))
    logger.debug('looping through ' + str(len(individuals)) + ' in samples in VCF')
    for i in range(start, end):
        variantsPerIndividual[individuals[i]] = dict()
        variantsPerIndividual[individuals[i]]['benign'] = list()
        variantsPerIndividual[individuals[i]]['pathogenic'] = list()
        variantsPerIndividual[individuals[i]]['vus'] = list()

        for variant in range(len(vcf['calldata/GT'])):
            if int(vcf['variants/CHROM'][variant].replace('chr', '')) != int(chromosome):
                logger.warning('wrong chromosome?')
                continue
            if 1 in vcf['calldata/GT'][variant][i]:
                c = int(vcf['variants/CHROM'][variant].replace('chr', ''))
                p = int(vcf['variants/POS'][variant])
                r = str(vcf['variants/REF'][variant])
                a = str(vcf['variants/ALT'][variant][0])
                if getGenesForVariant([c,p,r,a], ensemblRelease, gene) is None:
                    logger.warning('no gene for variant? ' + str([c,p,r,a]))
                    continue

                genotype = str(int(str(vcf['calldata/GT'][variant][i][0]) + str(vcf['calldata/GT'][variant][i][1]), 2))
                try:
                    if freeze == 5:
                        seqCenter = annoDF[annoDF['sample.id'] == individuals[i]]['CENTER'].iloc[0]
                    elif freeze == 8:
                        seqCenter = annoDF[annoDF['sample.id'] == individuals[i]]['seq_center'].iloc[0]
                    elif freeze == 9:
                        seqCenter = annoDF[annoDF['sample.id'] == individuals[i]]['seq_center'].iloc[0]
                except Exception as e:
                    seqCenter = "NA"
                try:
                    study = annoDF[annoDF['sample.id'] == individuals[i]]['study'].iloc[0]
                except Exception as e:
                    study = "NA"
                if (c, p, r, a) in benignVariants:
                    variantsPerIndividual[individuals[i]]['benign'].append(((c, p, r, a), genotype, seqCenter, study))
                elif (c, p, r, a) in pathogenicVariants:
                    variantsPerIndividual[individuals[i]]['pathogenic'].append(((c, p, r, a), genotype, seqCenter, study))
                # if not a known VUS, it is a VUS now
                else:
                    variantsPerIndividual[individuals[i]]['vus'].append(((c, p, r, a),  genotype, seqCenter, study))

    q.put(variantsPerIndividual)

def getGenesForVariant(variant, ensemblRelease, geneOfInterest):
    ensembl = pyensembl.EnsemblRelease(release=ensemblRelease)
    chrom = variant[0]
    if type(chrom) is str:
        chrom = chrom.split('chr')[1]
    pos = variant[1]
    try:
        #exons = ensembl.exons_at_locus(contig=int(chrom), position=int(pos))
        genes = ensembl.gene_names_at_locus(contig=int(chrom), position=int(pos))
        # TODO could get BRCA and other gene like ZAR1L?
        g_of_i = set()
        g_of_i.add(geneOfInterest)
        g = set(genes)
        intersectingGenes = g_of_i.intersection(g)
        if len(intersectingGenes) != 1:
            return None
        else:
            # this is pythonic way of returning member of singleton set
            return (intersectingGenes,)
    except Exception as e:
        logger.error('exception: ' + str(e))
        return None

def findIndividualsPerCooccurrence(variantsPerIndividual, ensemblRelease, phased, gene):

    individualsPerPathogenicCooccurrence = defaultdict(list)
    n = defaultdict(int)
    k = defaultdict(int)

    for individual in variantsPerIndividual:
        vusVarList = list(variantsPerIndividual[individual]['vus'])
        pathVarList = list(variantsPerIndividual[individual]['pathogenic'])

        for v in vusVarList:
            n[tuple(v[0])] += 1

        vusCrossPath = list(itertools.product(vusVarList, pathVarList))
        for cross in vusCrossPath:
            if sameGeneSameParent(cross[0], cross[1], phased, ensemblRelease, gene):
                k[tuple(cross[0][0])] += 1
                individualsPerPathogenicCooccurrence[(tuple(cross[0][0]), tuple(cross[1][0]))].append(individual)

    return individualsPerPathogenicCooccurrence, n, k

def sameGeneSameParent(vus, path, phased, ensemblRelease, gene):

    if not phased:
        return getGenesForVariant(vus[0], ensemblRelease, gene) == \
               getGenesForVariant(path[0], ensemblRelease, gene)
    else:
        # looking for vus in cis with path
        # if vus is 1|1, then it's either in cis or both in cis and in trans with path
        # else if vus and path are on opposite chromosomes
        return (getGenesForVariant(vus[0], ensemblRelease, gene) == getGenesForVariant(path[0], ensemblRelease, gene)) \
               and \
                ((vus[1] == '3') or ((vus[1] == '2' and path[1] == '1') or (vus[1] == '1' and path[1] == '2')))

                #((vus[1] == '1|1') or ((vus[1] == '1|0' and path[1] == '0|1') or (vus[1] == '0|1' and path[1] == '1|0')))

if __name__ == "__main__":
    main()


