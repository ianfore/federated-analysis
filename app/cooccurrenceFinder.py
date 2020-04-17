import pandas
import itertools
from multiprocessing.pool import ThreadPool
from multiprocessing import cpu_count
import time
import sys
import json
from collections import defaultdict
import pyensembl
import numpy as np
import os
from pandas.api.types import is_numeric_dtype
import argparse
import logging
import ast

logger = logging.getLogger()

defaultLogLevel = "WARNING"

def countColumnsAndMetaRows(fileName):
    '''The header line names the 8 fixed, mandatory columns.These columns are as follows:
    1. # CHROM
    2. POS
    3. ID
    4. REF
    5. ALT
    6. QUAL
    7. FILTER
    8. INFO
    If genotype data is present in the file, these are followed by a FORMAT column header, then an arbitrary number
    of sample IDs.Duplicate sample IDs are not allowed.The header line is tab - delimited.'''
    metaRowCount = 0
    with open(fileName, 'r') as f:
        for line in f:
            if line.startswith('##'):
                metaRowCount += 1
            elif line.startswith('#CHROM'):
                if 'FORMAT' in line:
                    preIDcolumnCount = 9
                else:
                    preIDcolumnCount = 8
                break
    f.close()
    return metaRowCount, preIDcolumnCount


classStrings = { 'Pathogenic':[ 'Pathogenic' ], 'Benign':[ 'Benign', 'Likely benign' ],
                 'Unknown': [ 'Uncertain significance', '-']}
sigColName = 'Clinical_significance_ENIGMA'
coordinateColumnBase = 'Genomic_Coordinate_hg'
alleleFrequencyName = 'Allele_frequency_ExAC'
DATA_DIR='/data/'
brcaFileName = DATA_DIR + 'brca-variants.tsv'
variantsPerIndividualFileName = DATA_DIR + 'variantsPerIndividual.json'

os.environ['PYENSEMBL_CACHE_DIR'] = DATA_DIR + 'pyensembl-cache'



# p2 = P(VUS is pathogenic and patient carries a pathogenic variant in trans) (arbitrarily set by goldgar et al)
# Integrated Evaluation of DNA Sequence Variants of Unknown Clinical Significance: Application to BRCA1 and BRCA2
brca1_p2 = 0.0001
brca2_p2 = 0.001

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

    parser.add_argument("vcf_filename", type=str,
                        help="name of file containing VCF data")

    parser.add_argument("output_filename", type=str,
                        help="name of JSON-formatted output file")

    parser.add_argument("--h", dest="h", help="Human genome version (37 or 38). Default=37", default=37)

    parser.add_argument("--e", dest="e", help="Ensembl version - 75 (for 37) or 99 (for 38). Default=75", default=75)

    parser.add_argument("--c", dest="c", help="List of chromosomes of interest. Default=[13,17]", default=[13,17])

    parser.add_argument("--g", dest="g", help="List of genes of interest. Default=['BRCA1', 'BRCA2']", default=['BRCA1',\
                                                                                                                'BRCA2'])

    parser.add_argument("--p", dest="p", help="Phased (boolean). Default=False", default=False)

    parser.add_argument("--s", dest="s", help="Save variants per individual to file. Default=False", default=False)


    parser.add_argument("--t", dest="t", help="Thread count. Default 1", default=1)


    parser.add_argument("--log", dest="logLevel", help="Logging level. Default=%s" %
                                                       defaultLogLevel, default=defaultLogLevel)

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

    # TODO for some reason the parsing of ['BRCA1','BRCA2'] works on brcaexchange-dev and my mac but not on crims
    # crimson complains it's a malformed string in eval
    # so this is just extra logic for crimson (using eval vs ast.literal_eval to convert from string)
    if isinstance(options.g, str):
        g_options = list(eval(options.g))
    else:
        g_options = list(options.g)
    c_options = list(eval(options.c))
    p_options = bool(eval(options.p))
    s_options = bool(eval(options.s))
    t_options = int(options.t)
    h_options = int(options.h)
    e_options = int(options.e)

    run(h_options, e_options, c_options, g_options, p_options, DATA_DIR + options.vcf_filename,
        DATA_DIR + options.output_filename, s_options, t_options)


def printUsage(args):
    sys.stderr.write('cooccurrenceFinder.py <genome-version> <ensembl-release> "[chr list]" "[gene list] True|False')

def run(hgVersion, ensemblRelease, chromosomes, genes, phased, vcfFileName, outputFileName, saveVarsPerIndivid, threadCount):

    print('hgversion = ' + str(hgVersion))
    print('ensembl = ' + str(ensemblRelease))
    print('chroms = ' + str(chromosomes))
    print('genes = ' + str(genes))
    print('phased = ' + str(phased))

    myVCFMetadataLines, myVCFskipCols = countColumnsAndMetaRows(vcfFileName)

    print('reading BRCA data from ' + brcaFileName)
    t = time.time()
    brcaDF, pathogenicVariants, benignVariants, unknownVariants = \
        findVariantsInBRCA(brcaFileName, classStrings, hgVersion)
    elapsed_time = time.time() - t
    print('elapsed time in findVariantsInBRCA() ' + str(elapsed_time))

    print('reading VCF data from ' + vcfFileName)
    t = time.time()
    vcfData = readVCFFile(vcfFileName, myVCFMetadataLines, chromosomes)
    elapsed_time = time.time() - t
    print('elapsed time in readVCFFile() ' + str(elapsed_time))


    print('finding variants per individual')
    t = time.time()
    results = list()
    variantsPerIndividual = dict()
    pool = ThreadPool(processes=threadCount)
    for i in range(threadCount):
        results.append(pool.apply_async(findVariantsPerIndividual, args=(vcfData, benignVariants, pathogenicVariants,
                                                                 myVCFskipCols, threadCount, i, phased)))
    for result in results:
        result.wait()
        variantsPerIndividual.update(result.get())
    elapsed_time = time.time() - t
    print('elapsed time in findVariantsPerIndividual() ' + str(elapsed_time))

    # find vus with genotype 1|1
    # TODO de-dup these!!
    print('finding homozygous  individuals per vus')
    homozygousPerVus = countHomozygousPerVus(variantsPerIndividual, brcaDF, hgVersion, ensemblRelease, genes)


    if saveVarsPerIndivid:
        print('saving variantsPerIndividual to ' + variantsPerIndividualFileName)
        with open(variantsPerIndividualFileName, 'w') as f:
            json.dump(variantsPerIndividual, f, cls=NpEncoder)
        f.close()

    '''print('reading in variantsPerIndividual from ' + variantsPerIndividualFileName)
    with open(variantsPerIndividualFileName) as f:
        variantsPerIndividual = json.load(f, cls=NpDecoder)
    f.close()'''


    # TODO make this multi-threaded (cpu is pegged at 99% for this method)
    # TODO or figure out vectorization!
    print('finding individuals per cooc')
    t = time.time()
    individualsPerPathogenicCooccurrence, n, k = findIndividualsPerCooccurrence(variantsPerIndividual, ensemblRelease,
                                                                                phased, genes)
    elapsed_time = time.time() - t
    print('elapsed time in findIndividualsPerCooccurrence() ' + str(elapsed_time))

    # TODO make the program idempotent (save/read data frames and take cli arg to figure out)
    # p1 = P(VUS is benign and patient carries a pathogenic variant in trans)
    numBenignWithPath = 0
    for cooc in individualsPerPathogenicCooccurrence:
        numBenignWithPath += len(individualsPerPathogenicCooccurrence[cooc])
    total = len(variantsPerIndividual)
    p1 =  0.5 * numBenignWithPath / total

    print('putting all the data together per vus')
    dataPerVus = calculateLikelihood(individualsPerPathogenicCooccurrence, p1, n, k)

    jsonOutputList = [dataPerVus, homozygousPerVus]
    print('saving final VUS data  to ' + outputFileName)
    with open(outputFileName, 'w') as f:
        json.dump(jsonOutputList, f, cls=NpEncoder)
    f.close()

def getGnomadData(brcaDF, vus, hgVersion):
    # TODO write a unit test
    # 13:g.32393468:C>CT
    hgString = 'chr' + str(vus[0][0]) + ':g.' + str(vus[0][1]) + ':' + str(vus[0][2]) + '>' + str(vus[0][3])
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

# TODO merge this method with findVariantsPerIndividual()
def countHomozygousPerVus(variantsPerIndividual, brcaDF, hgVersion, ensemblRelease, genesOfInterest):
    homoZygousPerVus = defaultdict(list)
    for individual in variantsPerIndividual:
        for vus in variantsPerIndividual[individual]['vus']:
            if (vus[1] == '1|1' or vus[1] == '1/1') and (getGenesForVariant(vus[0], ensemblRelease, genesOfInterest)):
                if str(vus) not in homoZygousPerVus:
                    homoZygousPerVus[str(vus)].append(0)
                    maxFreq = getGnomadData(brcaDF, vus, hgVersion)
                    homoZygousPerVus[str(vus)].append(maxFreq)
                homoZygousPerVus[str(vus)][0] += 1
    return homoZygousPerVus


def calculateLikelihood(pathCoocs, p1, n, k):

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
            print("unknown chromosome: " + str(vus[0]))
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
        #vus = repr(eval(cooc)[0])
        vus = cooc[0]
        #pathVarsPerVus[vus].append(eval(cooc)[1])
        pathVarsPerVus[vus].append(cooc[1])

    # put it all together in a single dict
    dataPerVus = dict()
    for vus in likelihoodRatios:
        if vus[0] == 13:
            p2 = brca2_p2
        elif vus[0] == 17:
            p2 = brca1_p2
        else:
            print("unknown chromosome: " + str(vus[0]))
            continue
        data = [p1, p2, n[vus], k[vus], likelihoodRatios[vus], pathVarsPerVus[vus]]
        dataPerVus[str(vus)] = data

    return dataPerVus


def readVCFFile(fileName, numMetaDataLines, chromosomes):
    # #CHROM  POS             ID      REF     ALT     QUAL    FILTER  INFO            FORMAT  0000057940      0000057950
    # 10      89624243        .       A       G       .       .       AF=1.622e-05    GT      0/0             0/0
    # 0/0 => does not have variant on either chromosome (homozygous negative)
    # 0/1  => has variant on 1 chromosome (heterozygous positive)
    # 1/1 =>  has variant on both chromosomes (homozygous positive)
    df = pandas.read_csv(fileName, sep='\t', skiprows=numMetaDataLines, dtype={'POS':int}, header=0)
    # this creates a bug: df = df[df.apply(lambda r: r.str.contains('1/1').any() or r.str.contains('0/1').any(), axis=1)]
    # filter chromosomes in CHROMOSOMES here
    df.columns = df.columns.str.replace('#', '')

    #if df.CHROM.dtype is not int:
    if not is_numeric_dtype(df['CHROM']):
        df['CHROM'] = df['CHROM'].str.replace('chr', '')
        df['CHROM'] = pandas.to_numeric(df['CHROM'])
    chromsDF = df[df.CHROM.isin(chromosomes)]
    # the index in the above operation is no longer contiguous from 0, so we need to reset for future operations on df
    chromsDF.index = np.arange(0, len(chromsDF))

    return chromsDF

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


def findVariantsPerIndividual(vcfDF, benignVariants, pathogenicVariants, skipCols, nThreads, threadID, phased):

    # find mutations
    # CHROM  POS             ID      REF     ALT     QUAL    FILTER  INFO            FORMAT  0000057940      0000057950
    # 10      89624243        .       A       G       .       .       AF=1.622e-05    GT      0/0             0/0
    # 0/0 => does not have variant on either strand (homozygous negative)
    # 0/1  => has variant on 1 strand (heterozygous positive)
    # 1/1 =>  has variant on both strands (homozygous positive)

    # 1|0 => paternal has variant, maternal does not
    # 0|1 => paternal does not have variant, maternal does
    # 1|1 => both paternal and maternal have variant
    # 0|0 => neither paternal nor maternal have variant
    variantsPerIndividual = dict()
    numColumns = len(vcfDF.columns)
    numIndividuals = len(vcfDF.columns) - skipCols
    partitionSize = int(numIndividuals / nThreads)
    start = threadID * partitionSize + skipCols
    if threadID == nThreads - 1:
        end = numColumns
    else:
        end = skipCols + start + partitionSize

    # get list of individuals
    individuals = vcfDF.columns[start:end]

    # iterate through columns (not rows! iterows() takes 20x longer b/c pandas are stored column-major)
    for individual in individuals:
        variantsPerIndividual[individual] = dict()
        variantsPerIndividual[individual]['benign'] = list()
        variantsPerIndividual[individual]['pathogenic'] = list()
        variantsPerIndividual[individual]['vus'] = list()

        if phased:
            includeList = ['1|0', '0|1', '1|1']
        else:
            includeList = ['1/0', '0/1', '1/1']

        listOfVariantIndices = list()
        for x in includeList:
            listOfVariantIndices += list(np.where(vcfDF[individual] == x)[0])

        for i in listOfVariantIndices:
            try:
                var = (int(vcfDF.loc[i, 'CHROM']), int(vcfDF.loc[i, 'POS']), str(vcfDF.loc[i, 'REF']), str(vcfDF.loc[i, 'ALT']))
                alleles = vcfDF.loc[i, individual]
                if var in benignVariants:
                    variantsPerIndividual[individual]['benign'].append((var, alleles))
                elif var in pathogenicVariants:
                    variantsPerIndividual[individual]['pathogenic'].append((var, alleles))
                # if not a known VUS, it is a VUS now
                else:
                    variantsPerIndividual[individual]['vus'].append((var, alleles))
            except Exception as e:
                print("exception for index " + str(i) + " of individual " + str(individual))
                continue

    return variantsPerIndividual


def getGenesForVariant(variant, ensemblRelease, genesOfInterest):
    ensembl = pyensembl.EnsemblRelease(release=ensemblRelease)
    chrom = variant[0]
    if type(chrom) is str:
        chrom = chrom.split('chr')[1]
    pos = variant[1]
    try:
        genes = ensembl.gene_names_at_locus(contig=int(chrom), position=int(pos))
        # TODO could get BRCA and other gene like ZAR1L?
        g_of_i = set(genesOfInterest)
        g = set(genes)
        intersectingGenes = g_of_i.intersection(g)
        if len(intersectingGenes) != 1:
            return None
        else:
            # this is pythonic way of returning member of singleton set
            return (intersectingGenes,)
    except Exception as e:
        print('exception: ' + str(e))
        return None

def findIndividualsPerCooccurrence(variantsPerIndividual, ensemblRelease, phased, genes):

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
            if sameGeneSameParent(cross[0], cross[1], phased, ensemblRelease, genes):
                k[tuple(cross[0][0])] += 1
                individualsPerPathogenicCooccurrence[(tuple(cross[0][0]), tuple(cross[1][0]))].append(individual)

    return individualsPerPathogenicCooccurrence, n, k

def sameGeneSameParent(vus, path, phased, ensemblRelease, genes):

    if not phased:
        return getGenesForVariant(vus[0], ensemblRelease, genes) == \
               getGenesForVariant(path[0], ensemblRelease, genes)
    else:
        # looking for vus in cis with path
        # if vus is 1|1, then it's either in cis or both in cis and in trans with path
        # else if vus and path are on opposite chromosomes
        return (getGenesForVariant(vus[0], ensemblRelease, genes) == getGenesForVariant(path[0], ensemblRelease, genes)) \
               and \
               ((vus[1] == '1|1') or ((vus[1] == '1|0' and path[1] == '0|1') or (vus[1] == '0|1' and path[1] == '1|0')))

if __name__ == "__main__":
    main()


