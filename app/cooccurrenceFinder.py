import pandas
import itertools
import time
import sys
import json
from collections import defaultdict
import pyensembl
import numpy as np
import os
import argparse
import logging
import allel
from multiprocessing import Process, Queue, cpu_count


logger = logging.getLogger()
defaultLogLevel = "DEBUG"

# p2 = P(VUS is pathogenic and patient carries a pathogenic variant in trans) (arbitrarily set by goldgar et al)
# Integrated Evaluation of DNA Sequence Variants of Unknown Clinical Significance: Application to BRCA1 and BRCA2
brca1_p2 = 0.0001
brca2_p2 = 0.001

classStrings = { 'Pathogenic':[ 'Pathogenic' ], 'Benign':[ 'Benign', 'Likely benign' ],
                 'Unknown': [ 'Uncertain significance', '-']}
sigColName = 'Clinical_significance_ENIGMA'
coordinateColumnBase = 'Genomic_Coordinate_hg'
alleleFrequencyName = 'Allele_frequency_ExAC'
#DATA_DIR='/data/'
#brcaFileName = DATA_DIR + 'brca-variants.tsv'
#brcaFileName = 'brca-variants.tsv'
#variantsPerIndividualFileName = DATA_DIR + 'variantsPerIndividual.json'
variantsPerIndividualFileName = 'variantsPerIndividual.json'

#os.environ['PYENSEMBL_CACHE_DIR'] = DATA_DIR + 'pyensembl-cache'
#os.environ['PYENSEMBL_CACHE_DIR'] = 'pyensembl-cache'


COMMON_VARIANT_CUTOFF_FREQUENCY=0.01

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

    parser.add_argument("--v", dest="v", help="name of file containing VCF data, default=None", default=None)

    parser.add_argument("--o", dest="o", help="name of JSON-formatted output file, default=None", default=None)

    parser.add_argument("--h", dest="h", help="Human genome version (37 or 38). Default=None", default=None)

    parser.add_argument("--e", dest="e", help="Ensembl version - 75 (for 37) or 99 (for 38). Default=None", default=None)

    parser.add_argument("--c", dest="c", help="Chromosome of interest. Default=None", default=None)

    parser.add_argument("--g", dest="g", help="Gene of interest. Default=None", default=None)

    parser.add_argument("--p", dest="p", help="Phased (boolean). Default=False", default='True')

    parser.add_argument("--s", dest="s", help="Save variants per individual to file. Default=False", default='False')

    parser.add_argument("--i", dest="i", help="Include pathog vars per VUS in report. Default=False", default='True')

    parser.add_argument("--a", dest="a", help="calculate allele freqs for homozygous. Default=False", default='True')

    parser.add_argument("--n", dest="n", help="Number of processes. Default=1", default=cpu_count())

    parser.add_argument("--b", dest="b", help="BRCA variants file. Default=brca-variants", default=None)

    parser.add_argument("--d", dest="d", help="directory containing pyensembl-cache. Default=None", default=None)

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

    '''if isinstance(options.g, str):
        g_options = list(eval(options.g))
    else:
        g_options = list(options.g)'''
    g_options = options.g
    b_options = options.b
    v_options = options.v
    o_options = options.o
    c_options = options.c
    d_options = options.d
    p_options = bool(eval(options.p))
    s_options = bool(eval(options.s))
    i_options = bool(eval(options.i))
    a_options = bool(eval(options.a))
    h_options = int(options.h)
    e_options = int(options.e)
    n_options = int(options.n)

    print(options)

    '''run(h_options, e_options, c_options, g_options, p_options, DATA_DIR + options.vcf_filename,
        DATA_DIR + options.output_filename, s_options, i_options, a_options, n_options)'''
    run(h_options, e_options, c_options, g_options, p_options, v_options, o_options, s_options, i_options, a_options,
        n_options, b_options, d_options)

def run(hgVersion, ensemblRelease, chromosome, gene, phased, vcfFileName, outputFileName, saveVarsPerIndivid,
        includePaths, calculateAlleleFreqs, numProcs, brcaFileName, pyensemblDir):

    logger.info('setting pyensembl dir to ' + pyensemblDir)
    os.environ['PYENSEMBL_CACHE_DIR'] = pyensemblDir

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
    processList = list()
    for i in range(numProcs):
        p = Process(target=findVarsPerIndividual, args=(q, vcf, benignVariants, pathogenicVariants, chromosome, i, numProcs,))
        p.start()
        processList.append(p)
    logger.info('joining results from forked threads')
    variantsPerIndividual = dict()
    for i in range(numProcs):
        variantsPerIndividual.update(q.get())
    for i in range(numProcs):
        processList[i].join()
    logger.info('elapsed time in findVariantsPerIndividual() ' + str(time.time() -t))


    if calculateAlleleFreqs:
        logger.info('finding homozygous  individuals per vus')
        t = time.time()
        homozygousPerVus = countHomozygousPerVus(variantsPerIndividual, brcaDF, hgVersion, ensemblRelease, gene)
        logger.info('elapsed time in countHomozygousPerVus() ' + str(time.time() -t))

    if saveVarsPerIndivid:
        logger.info('saving variantsPerIndividual to ' + variantsPerIndividualFileName)
        with open(variantsPerIndividualFileName, 'w') as f:
            json.dump(variantsPerIndividual, f, cls=NpEncoder)
        f.close()

    logger.info('finding individuals per cooc')
    t = time.time()
    individualsPerPathogenicCooccurrence, n, k = findIndividualsPerCooccurrence(variantsPerIndividual, ensemblRelease,
                                                                                phased, gene)
    logger.info('elapsed time in findIndividualsPerCooccurrence() ' + str(time.time() -t))


    # TODO check this math!
    logger.info('calculating p1')
    # p1 = P(VUS is benign and patient carries a path variant in trans) = 0.5 * overall freq of path muts in cohort

    # calculate total number of benign, pathogenic, and vus variants in cohort

    '''numBenignWithPath = 0
    for cooc in individualsPerPathogenicCooccurrence:
        numBenignWithPath += len(individualsPerPathogenicCooccurrence[cooc])'''

    totalBenign = set()
    totalPathogenic = set()
    totalVUS = set()

    for i in variantsPerIndividual:
        for b in variantsPerIndividual[i]['benign']:
            totalBenign.add(b)
        for p in variantsPerIndividual[i]['pathogenic']:
            totalPathogenic.add(p)
        for v in variantsPerIndividual[i]['vus']:
            totalVUS.add(v)

    numBenign = len(totalBenign)
    numPathogenic = len(totalPathogenic)
    numVUS = len(totalVUS)

    print('total benign: ' + str(numBenign))
    print('total pathogenic: ' + str(numPathogenic))
    print('total vus: ' + str(numVUS))

    cohortSize = len(variantsPerIndividual)
    p1 =  0.5 * numBenign / cohortSize

    logger.info('putting all the data together per vus')
    dataPerVus = calculateLikelihood(individualsPerPathogenicCooccurrence, p1, n, k, includePaths, brcaDF, hgVersion, cohortSize)

    if calculateAlleleFreqs:
        jsonOutputList = [dataPerVus, homozygousPerVus]
    else:
        jsonOutputList = dataPerVus
    logger.info('saving final VUS data  to ' + outputFileName)
    with open(outputFileName, 'w') as f:
        json.dump(jsonOutputList, f, cls=NpEncoder)
    f.close()

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

    return (maxPopulation, maxFrequency)

def countHomozygousPerVus(variantsPerIndividual, brcaDF, hgVersion, ensemblRelease, geneOfInterest):
    homoZygousPerVus = defaultdict(list)

    for individual in variantsPerIndividual:
        for vus in variantsPerIndividual[individual]['vus']:
            if (vus[1] == '3') and (getGenesForVariant(vus[0], ensemblRelease, geneOfInterest)):
                if str(vus) not in homoZygousPerVus:
                    homoZygousPerVus[str(vus)].append(0)
                    maxPop, maxFreq = getGnomadData(brcaDF, vus[0], hgVersion)
                    homoZygousPerVus[str(   vus)].append([maxPop, maxFreq])
                homoZygousPerVus[str(vus)][0] += 1

    cohortSize = len(variantsPerIndividual)
    for vus in homoZygousPerVus:
        maxPopFreq = homoZygousPerVus[vus][1][1]
        cohortFreq = float(homoZygousPerVus[vus][0])/ float(cohortSize)
        homoZygousPerVus[vus].append(float(cohortFreq))
        if cohortFreq < COMMON_VARIANT_CUTOFF_FREQUENCY or (maxPopFreq != 0 and maxPopFreq < COMMON_VARIANT_CUTOFF_FREQUENCY):
            homoZygousPerVus[vus].append('RARE')

    return homoZygousPerVus

def calculateLikelihood(pathCoocs, p1, n, k, includePathogenicVariants, brcaDF, hgVersion, cohortSize):

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
        maxPop, maxPopFreq = getGnomadData(brcaDF, vus, hgVersion)
        cohortFreq = float(n[vus]) / float(cohortSize)
        if includePathogenicVariants:
            data = [[p1, p2, n[vus], k[vus], likelihoodRatios[vus]], [maxPop, maxPopFreq, cohortFreq], pathVarsPerVus[vus]]
        else:
            data = [[p1, p2, n[vus], k[vus], likelihoodRatios[vus]], [maxPop, maxPopFreq, cohortFreq]]
        if maxPopFreq < COMMON_VARIANT_CUTOFF_FREQUENCY or cohortFreq < COMMON_VARIANT_CUTOFF_FREQUENCY:
            data.append('RARE')
        dataPerVus[str(vus)] = data

    return dataPerVus

def readVCFFile(vcfFileName):
    return allel.read_vcf(vcfFileName)

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

def findVarsPerIndividual(q, vcf, benignVariants, pathogenicVariants, chromosome, threadID, numProcesses):

    variantsPerIndividual = dict()

    individuals = list(vcf['samples'])

    # calculate start and stop samples per threadID
    n = len(individuals)
    partitionSizes = divide(n, numProcesses)
    start,end = getStartAndEnd(partitionSizes, threadID)

    logger.info('threadID = ' + str(threadID) + ' processing from ' + str(start) + ' to ' + str(end))

    logger.debug('looping through ' + str(len(individuals)) + ' in samples in VCF')
    #for i in range(len(individuals)):
    for i in range(start, end):
        #logger.debug('looking at individual ' + str(individuals[i]))
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
                genotype = str(int(str(vcf['calldata/GT'][variant][i][0]) + str(vcf['calldata/GT'][variant][i][1]), 2))
                if (c, p, r, a) in benignVariants:
                    variantsPerIndividual[individuals[i]]['benign'].append(((c, p, r, a), genotype))
                elif (c, p, r, a) in pathogenicVariants:
                    variantsPerIndividual[individuals[i]]['pathogenic'].append(((c, p, r, a), genotype))
                # if not a known VUS, it is a VUS now
                else:
                    variantsPerIndividual[individuals[i]]['vus'].append(((c, p, r, a), genotype))
    #return variantsPerIndividual
    q.put(variantsPerIndividual)

def getGenesForVariant(variant, ensemblRelease, geneOfInterest):
    ensembl = pyensembl.EnsemblRelease(release=ensemblRelease)
    chrom = variant[0]
    if type(chrom) is str:
        chrom = chrom.split('chr')[1]
    pos = variant[1]
    try:
        genes = ensembl.gene_names_at_locus(contig=int(chrom), position=int(pos))
        # TODO could get BRCA and other gene like ZAR1L?
        g_of_i = set(geneOfInterest)
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


