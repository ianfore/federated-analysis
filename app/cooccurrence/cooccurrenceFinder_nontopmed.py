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

def main():
    # main just parses the CLI and then calls the run() method with appropriate args

    parser = argparse.ArgumentParser(usage="cooccurrenceFinder args [options]")
    parser.add_argument("--vcf", dest="vcf", help="name of file containing VCF data, default=None", default=None)
    parser.add_argument("--ipv", dest="ipv", help="ipv file name, default=ipv.json", default='ipv.json')
    parser.add_argument("--vpi", dest="vpi", help="vpi file name, default=vpi.json", default='vpi.json')
    parser.add_argument("--out", dest="out", help="output file name, default=out.json", default='out.json')
    parser.add_argument("--all", dest="all", help="all vars file name, default=all.json", default='all.json')
    parser.add_argument("--h", dest="h", help="Human genome version (37 or 38). Default=None", default=None)
    parser.add_argument("--e", dest="e", help="Ensembl version - 75 (for 37) or 99 (for 38). Default=None", default=None)
    parser.add_argument("--c", dest="c", help="Chromosome of interest. Default=None", default=None)
    parser.add_argument("--g", dest="g", help="Gene of interest. Default=None", default=None)
    parser.add_argument("--p", dest="p", help="Phased (boolean). Default=False", default='True')
    parser.add_argument("--n", dest="n", help="Number of processes. Default=cpu_count", default=cpu_count())
    parser.add_argument("--b", dest="b", help="BRCA variants file. Default=brca-variants", default=None)
    parser.add_argument("--d", dest="d", help="directory containing pyensembl-cache. Default=/var/tmp/pyensembl-cache", default='/var/tmp/pyensembl-cache')
    parser.add_argument("--r", dest="r", help="Rare frequency cutoff. Default=0.01", default=0.01)
    parser.add_argument("--pf", dest="pf", help="Pathology input file. Default=None", default=None)
    parser.add_argument("--log", dest="logLevel", help="Logging level. Default=%s" % defaultLogLevel, default=defaultLogLevel)
    options = parser.parse_args()

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

    g_options = options.g
    b_options = options.b
    v_options = options.vcf
    ipv_options = options.ipv
    vpi_options = options.vpi
    all_options = options.all
    out_options = options.out
    c_options = options.c
    d_options = options.d
    p_options = bool(eval(options.p))
    h_options = int(options.h)
    e_options = int(options.e)
    n_options = int(options.n)
    r_options = float(options.r)
    pf_options = options.pf

    print(options)


    run(h_options, e_options, c_options, g_options, p_options, v_options, n_options, b_options, d_options, r_options,
        ipv_options, vpi_options, all_options, out_options, pf_options)

def run(hgVersion, ensemblRelease, chromosome, gene, phased, vcfFileName, numProcs,
        brcaFileName, pyensemblDir, rareCutoff, ipvFileName, vpiFileName, allVariantsFileName, outputFileName, pathologyFile):

    logger.info('setting pyensembl dir to ' + pyensemblDir)
    os.environ['PYENSEMBL_CACHE_DIR'] = '/var/tmp/pyensembl-cache'

    logger.info('reading BRCA data from ' + brcaFileName)
    t = time.time()
    brcaDF, pathogenicVariants, benignVariants, unknownVariants = findVariantsInBRCA(brcaFileName, classStrings, hgVersion)
    logger.info('elapsed time in findVariantsInBRCA() ' + str(time.time() -t))

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
                                                        ensemblRelease, i, numProcs, ))
        p.start()
        processList.append(p)
    logger.info('joining results from forked threads')
    variantsPerIndividual = dict()
    individualsPerVariant = dict()
    for i in range(numProcs):
        variantsPerIndividual.update(q.get())
        individualsPerVariant.update(w.get())
    for i in range(numProcs):
        processList[i].join()
    logger.info('elapsed time in findVariantsPerIndividual() ' + str(time.time() -t))

    # find individuals per variant
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
        for p in variantsPerIndividual[individual]['pathogenic']:
            v = str((p[0][0], p[0][1], p[0][2], p[0][3]))
            if not v in individualsPerVariant:
                individualsPerVariant[v] = {'heterozygous individuals': set(),
                                            'homozygous individuals': set()}
            if p[1] == '1' or p[1] == '2':
                individualsPerVariant[v]['heterozygous individuals'].add(individual)
            elif p[1] == '3':
                individualsPerVariant[v]['homozygous individuals'].add(individual)
        for vus in variantsPerIndividual[individual]['vus']:
            v = str((vus[0][0], vus[0][1], vus[0][2], vus[0][3]))
            if not v in individualsPerVariant:
                individualsPerVariant[v] = {'heterozygous individuals': set(),
                                            'homozygous individuals': set()}
            if vus[1] == '1' or vus[1] == '2':
                individualsPerVariant[v]['heterozygous individuals'].add(individual)
            elif vus[1] == '3':
                individualsPerVariant[v]['homozygous individuals'].add(individual)
    cohortSize = len(variantsPerIndividual)
    individualsPerVariant = addVariantInfo(individualsPerVariant, vcf, chromosome, brcaDF, hgVersion, cohortSize,
                                           ensemblRelease)

    # commenting out for biobank japan
    #logger.info('saving vpi to ' + vpiFileName)
    #with open(vpiFileName, 'w') as f:
    #    json.dump(variantsPerIndividual, f, cls=NpEncoder)
    #f.close()
    #logger.info('saving ipv to ' + ipvFileName)
    #with open(ipvFileName, 'w') as f:
    #    json.dump(individualsPerVariant, f, cls=NpEncoder)
    #f.close()


    logger.info('counting zygous individuals per vus')
    t = time.time()
    homozygousPerVus = countZygousPerVus(variantsPerIndividual, brcaDF, hgVersion, ensemblRelease, gene)
    logger.info('elapsed time in countZygousPerVus() ' + str(time.time() -t))


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
    numBenign = len(set(allVariants['benign']))
    p1 =  0.5 * numBenign / cohortSize
    logger.info('saving all variants to ' + allVariantsFileName)
    json_dump = json.dumps(allVariants, cls=NpEncoder)
    with open(allVariantsFileName, 'w') as f:
        f.write(json_dump)
    f.close()

    logger.info('putting all the data together per vus')
    dataPerVus = calculateLikelihood(individualsPerPathogenicCooccurrence, p1, n, k, brcaDF, hgVersion, cohortSize, rareCutoff)

    data_set = {"cooccurring vus": dataPerVus, "homozygous vus": homozygousPerVus}
    json_dump = json.dumps(data_set, cls=NpEncoder)

    logger.info('saving final VUS data  to ' + outputFileName)
    with open(outputFileName, 'w') as f:
        f.write(json_dump)
    f.close()

    logger.info('intersecting variants with pathology data')
    intersectionFile = '/data/' + str(chromosome) + '-intersection.json'
    intersectPathology(pathologyFile, data_set, individualsPerVariant, intersectionFile )

def intersectPathology(pathologyFile, data_set, ipv, intersectFile ):
    logger.info('reading data from ' + pathologyFile)
    pathologyDF = pandas.read_csv(pathologyFile, sep='\t', header=0)

    #variantsDF = pandas.DataFrame.from_dict(data_set)
    variantsDF = data_set

    #ipvDF = pandas.DataFrame.from_dict(ipv)
    ipvDF = ipv

    # determine total number of cases and controls for cohort frequency calculations
    numCases = 0
    numControls = 0
    for i in range(len(pathologyDF)):
        if pathologyDF.iloc[i]['Age at onset'] != 0:
            numCases += 1
        else:
            numControls += 1

    pathologyPerCoocIndividual = dict()
    for variant in variantsDF['cooccurring vus']:
        for pathogenicVariant in variantsDF['cooccurring vus'][variant]['pathogenic variants']:
            pv = str(tuple(pathogenicVariant))
            heterozygousIndividuals = ipvDF[pv]['heterozygous individuals']
            pathologyPerCoocIndividual[pv] = dict()
            pathologyPerCoocIndividual[pv]['pathologies'] = list()
            pathologyPerCoocIndividual[pv]['numCases'] = 0
            pathologyPerCoocIndividual[pv]['numControls'] = 0
            pathologies = dict()
            for hi in heterozygousIndividuals:
                hiInt = int(hi)
                row = pathologyDF.loc[pathologyDF['ID'] == hiInt]
                aao = row['Age at onset'].tolist()
                if aao:
                    pathologies['Age at onset'] = aao[0]
                    pathologyPerCoocIndividual[pv]['numCases']  += 1
                else:
                    pathologies['Age at onset'] = 0.0
                    pathologyPerCoocIndividual[pv]['numControls'] += 1
                pathologies['Ovarian cancer history'] = row['Ovarian cancer history'].tolist()
                pathologies['Bilateral breast cancer'] = row['Bilateral breast cancer'].tolist()
                pathologies['Tissue type (3 groups)'] = row['Tissue type (3 groups)'].tolist()
                pathologies['TMN classification / T'] = row['TMN classification / T'].tolist()
                pathologies['TNM classification / N'] = row['TNM classification / N'].tolist()
                pathologies['TNM classification / M'] = row['TNM classification / M'].tolist()
                pathologies['ER'] = row['ER'].tolist()
                pathologies['PgR'] = row['PgR'].tolist()
                pathologies['HER2'] = row['HER2'].tolist()

                pathologyPerCoocIndividual[pv]['pathologies'].append(pathologies)

            if numCases == 0:
                pathologyPerCoocIndividual[pv]['caseFreq'] = 0
            else:
                pathologyPerCoocIndividual[pv]['caseFreq'] = float(pathologyPerCoocIndividual[pv]['numCases'] )/float(numCases)

            if numControls == 0:
                pathologyPerCoocIndividual[pv]['controlFreq'] = 0
            else:
                pathologyPerCoocIndividual[pv]['controlFreq'] = float(pathologyPerCoocIndividual[pv]['numControls'] )/float(numControls)


    pathologyPerHomoIndividual = dict()
    for variant in variantsDF['homozygous vus']:
        homozygousIndividuals = ipvDF[variant]['homozygous individuals']
        pathologyPerHomoIndividual[variant] = dict()
        pathologyPerHomoIndividual[variant]['pathologies'] = list()
        pathologyPerHomoIndividual[variant]['numCases'] = 0
        pathologyPerHomoIndividual[variant]['numControls'] = 0
        pathologies = dict()
        for hi in homozygousIndividuals:
            hiInt = int(hi)
            row = pathologyDF.loc[pathologyDF['ID'] == hiInt]
            aao = row['Age at onset'].tolist()
            if aao:
                pathologies['Age at onset'] = aao[0]
                pathologyPerHomoIndividual[variant]['numCases'] += 1
            else:
                pathologies['Age at onset'] = 0.0
                pathologyPerHomoIndividual[variant]['numControls'] += 1
            pathologies['Ovarian cancer history'] = row['Ovarian cancer history'].tolist()
            pathologies['Bilateral breast cancer'] = row['Bilateral breast cancer'].tolist()
            pathologies['Tissue type (3 groups)'] = row['Tissue type (3 groups)'].tolist()
            pathologies['TMN classification / T'] = row['TMN classification / T'].tolist()
            pathologies['TNM classification / N'] = row['TNM classification / N'].tolist()
            pathologies['TNM classification / M'] = row['TNM classification / M'].tolist()
            pathologies['ER'] = row['ER'].tolist()
            pathologies['PgR'] = row['PgR'].tolist()
            pathologies['HER2'] = row['HER2'].tolist()

            pathologyPerHomoIndividual[variant]['pathologies'].append(pathologies)

        if numCases == 0:
            pathologyPerHomoIndividual[variant]['caseFreq'] = 0
        else:
            pathologyPerHomoIndividual[variant]['caseFreq'] = float(pathologyPerHomoIndividual[variant]['numCases']) / float(numCases)
        if numControls == 0:
            pathologyPerHomoIndividual[variant]['controlFreq'] = 0
        else:
            pathologyPerHomoIndividual[variant]['controlFreq'] = float(pathologyPerHomoIndividual[variant]['numControls']) / float(
        numControls)

    pathologyPerAllIndividuals = dict()
    #pathologyPerAllIndividuals.update(pathologyPerHomoIndividual)
    #pathologyPerAllIndividuals.update(pathologyPerCoocIndividual)
    pathologyPerAllIndividuals['homozygous'] = pathologyPerHomoIndividual
    pathologyPerAllIndividuals['cooccurring'] = pathologyPerCoocIndividual

    json_dump = json.dumps(pathologyPerAllIndividuals, indent=4, sort_keys=True)
    with open(intersectFile, 'w') as f:
        f.write(json_dump)
    f.close()


def isExonic(ensemblRelease, chrom, pos):
    ensembl = pyensembl.EnsemblRelease(release=ensemblRelease)
    try:
        exons = ensembl.exons_at_locus(contig=int(chrom), position=int(pos))
    except Exception as e:
        logger.error('exception: ' + str(e))
        return None
    return len(exons) > 0

def addVariantInfo(individualsPerVariant, vcf, chromosome, brcaDF, hgVersion, cohortSize, ensemblRelease):
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
            maxPop, maxFreq, minPop, minFreq = getGnomadData(brcaDF, eval(v), hgVersion)
            individualsPerVariant[v]['maxPop'] = maxPop
            individualsPerVariant[v]['maxFreq'] = maxFreq
            individualsPerVariant[v]['minPop'] = minPop
            individualsPerVariant[v]['minFreq'] = minFreq
            individualsPerVariant[v]['cohortFreq'] = float(len(individualsPerVariant[v]['homozygous individuals']) + \
                len(individualsPerVariant[v]['heterozygous individuals']) ) / float(cohortSize)
            individualsPerVariant[v]['exonic'] = isExonic(ensemblRelease, c, p)

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


def countZygousPerVus(variantsPerIndividual, brcaDF, hgVersion, ensemblRelease, geneOfInterest):
    homozygousPerVus = dict()

    for individual in variantsPerIndividual:
        for vus in variantsPerIndividual[individual]['vus']:
            if (vus[1] == '3') and (getGenesForVariant(vus[0], ensemblRelease, geneOfInterest)):
                if str(vus[0]) not in homozygousPerVus:
                    homozygousPerVus[str(vus[0])] = dict()
                    homozygousPerVus[str(vus[0])]['nHom'] = 0
                    homozygousPerVus[str(vus[0])]['nHet'] = 0
                    maxPop, maxPopFreq, minPop, minPopFreq = getGnomadData(brcaDF, vus[0], hgVersion)
                    homozygousPerVus[str(vus[0])]['maxPop'] = maxPop
                    homozygousPerVus[str(vus[0])]['maxPopFreq'] = maxPopFreq
                homozygousPerVus[str(vus[0])]['nHom'] += 1
                for vus2 in variantsPerIndividual[individual]['vus']:
                    if vus2[1] == '1' or vus2[1] == '2':
                        homozygousPerVus[str(vus[0])]['nHet'] += 1

    cohortSize = len(variantsPerIndividual)
    for vus in homozygousPerVus:
        cohortFreq = float(homozygousPerVus[vus]['nHom'])/ float(cohortSize)
        homozygousPerVus[vus]['cohortFreq'] = float(cohortFreq)

    return homozygousPerVus

def calculateLikelihood(pathCoocs, p1, n, k, brcaDF, hgVersion, cohortSize, rareCutoff):

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
                'allele frequencies':{'maxPop':maxPop, 'maxPopFreq':maxPopFreq, 'cohortFreq':cohortFreq},
                'pathogenic variants': pathVarsPerVus[vus]}
        dataPerVus[str(vus)] = data

    return dataPerVus

def readVCFFile(vcfFileName):
    return allel.read_vcf(vcfFileName, fields=['samples', 'calldata/GT', 'variants/ALT', 'variants/CHROM',
                                               'variants/ID', 'variants/POS', 'variants/QUAL', 'variants/REF',
                                               'variants/INFO'])

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

def findVarsPerIndividual(q, w, vcf, benignVariants, pathogenicVariants, chromosome, gene, ensemblRelease,
                          threadID, numProcesses):
    '''infoFields = ['variants/ABE', 'variants/ABZ', 'variants/AC', 'variants/AF',
     'variants/AN', 'variants/ANN', 'variants/AVGDP', 'variants/BETA_IF', 'variants/BQZ',
     'variants/CYZ', 'variants/FIBC_I', 'variants/FIBC_P', 'variants/FILTER_PASS', 'variants/FLT20', 'variants/GC',
     'variants/GN', 'variants/HWEAF_P', 'variants/HWE_SLP_I', 'variants/HWE_SLP_P', 'variants/IBC_I', 'variants/IBC_P',
     'variants/ID', 'variants/IOR', 'variants/MAX_IF', 'variants/MIN_IF', 'variants/NM0', 'variants/NM1',
     'variants/NMZ', 'variants/NS_NREF', 'variants/QUAL', 'variants/STZ', 'variants/SVM']'''

    variantsPerIndividual = dict()
    individualsPerVariant = dict()
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
                continue
            if 1 in vcf['calldata/GT'][variant][i]:
                c = int(vcf['variants/CHROM'][variant].replace('chr', ''))
                p = int(vcf['variants/POS'][variant])
                r = str(vcf['variants/REF'][variant])
                a = str(vcf['variants/ALT'][variant][0])
                if getGenesForVariant([c,p,r,a], ensemblRelease, gene) is None:
                    continue

                genotype = str(int(str(vcf['calldata/GT'][variant][i][0]) + str(vcf['calldata/GT'][variant][i][1]), 2))
                if (c, p, r, a) in benignVariants:
                    variantsPerIndividual[individuals[i]]['benign'].append(((c, p, r, a), genotype))
                elif (c, p, r, a) in pathogenicVariants:
                    variantsPerIndividual[individuals[i]]['pathogenic'].append(((c, p, r, a), genotype))
                # if not a known VUS, it is a VUS now
                else:
                    variantsPerIndividual[individuals[i]]['vus'].append(((c, p, r, a),  genotype))

    q.put(variantsPerIndividual)
    w.put(individualsPerVariant)

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


