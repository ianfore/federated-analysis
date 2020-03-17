import pandas
import itertools
from multiprocessing.pool import ThreadPool
from multiprocessing import cpu_count
import time
import sys
import json
from collections import defaultdict
import ast
import pyensembl
import numpy as np
from functools import partial
import os
import ast
from pandas.api.types import is_numeric_dtype


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

#nThreads = cpu_count()
nThreads=1
classStrings = { 'Pathogenic':[ 'Pathogenic' ], 'Benign':[ 'Benign', 'Likely benign' ],
                 'Unknown': [ 'Uncertain significance', '-']}
sigColName = 'Clinical_significance_ENIGMA'
#DATA_DIR='/Users/jcasaletto/PycharmProjects/BIOBANK/federated-analysis/data'
DATA_DIR='/data'
brcaFileName = DATA_DIR + '/brca-variants.tsv'
vcfFileName = DATA_DIR + '/13-BreastCancer.shuffle.vcf'
#vcfFileName = DATA_DIR + '/topmed-test.vcf'
variantsPerIndividualFileName = DATA_DIR + '/variantsPerIndividual.json'
pathogenicCooccurrencesFileName = DATA_DIR + '/pathogenicCooccurrences.json'
benignCooccurrencesFileName = DATA_DIR + '/benignCooccurrences.json'
vusCooccurrencesFileName = DATA_DIR + '/vusCooccurrences.json'
vusFinalDataFileName = DATA_DIR + '/vusFinalData.json'
os.environ['PYENSEMBL_CACHE_DIR'] = DATA_DIR + '/pyensembl-cache'

myVCFMetadataLines, myVCFskipCols = countColumnsAndMetaRows(vcfFileName)


# p2 = P(VUS is pathogenic and patient carries a pathogenic variant in trans) (arbitrarily set by tavtigian et al)
p2 = 0.0001

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
        else:
            return super(NpEncoder, self).default(obj)

class NpDecoder(json.JSONDecoder):
    def default(self, obj):
        if isinstance(obj, list):
            return set(obj)
        else:
            return super(NpDecoder, self).default(obj)

def main():
    if len(sys.argv) != 5:
        printUsage(sys.argv)
        sys.exit(1)
    GENOME_VERSION=sys.argv[1]
    ENSEMBL_RELEASE=sys.argv[2]
    CHROMOSOMES = ast.literal_eval(sys.argv[3])
    GENES = ast.literal_eval(sys.argv[4])
    for i in  range(len(CHROMOSOMES)):
        try:
            CHROMOSOMES[i] = int(CHROMOSOMES[i])
        except Exception as e:
            print('exception parsing chromosome list: ' + str(e))
            sys.exit(2)

    combo(GENOME_VERSION, ENSEMBL_RELEASE, CHROMOSOMES, GENES)


def printUsage(args):
    sys.stderr.write('myCooccurrenceFinder.py <ensembl-release> "[chr list]" "[gene list]')

def combo(hgVersion, ensemblRelease, chromosomes, genes):

    print('reading BRCA data from ' + brcaFileName)
    t = time.time()
    count, pathogenicVariants, benignVariants, unknownVariants = \
        findVariantsInBRCA(brcaFileName, classStrings, sigColName, hgVersion)
    elapsed_time = time.time() - t
    print('elapsed time in findVariantsInBRCA() ' + str(elapsed_time))

    print('reading VCF data from ' + vcfFileName)
    t = time.time()
    vcfData = readVCFFile(vcfFileName, myVCFMetadataLines, chromosomes)
    elapsed_time = time.time() - t
    print('elapsed time in readVCFFile() ' + str(elapsed_time))

    print('slicing up VCF data by chromosome')
    t = time.time()
    variantsPerChromosome = findVariantsPerChromosome(vcfData, chromosomes)
    elapsed_time = time.time() - t
    print('elapsed time in variantsPerChromosome() ' + str(elapsed_time))

    print('slicing up VCF data by gene')
    t = time.time()
    benignPerGene, pathogenicPerGene, vusPerGene = findVariantsPerGene(variantsPerChromosome,
                                                                       benignVariants, pathogenicVariants,
                                                                       unknownVariants, ensemblRelease, genes)
    elapsed_time = time.time() - t
    print('elapsed time in variantsPerGene() ' + str(elapsed_time))

    print('finding variants per individual')
    t = time.time()
    results = list()
    variantsPerIndividual = dict()
    pool = ThreadPool(processes=nThreads)
    for i in range(nThreads):
        results.append(pool.apply_async(findVariantsPerIndividual, args=(vcfData, benignVariants, pathogenicVariants,
                                                                 unknownVariants, myVCFskipCols, nThreads, i)))
    for result in results:
        result.wait()
        variantsPerIndividual.update(result.get())
    elapsed_time = time.time() - t
    print('elapsed time in findVariantsPerIndividual() ' + str(elapsed_time))

    print('saving variantsPerIndividual to ' + variantsPerIndividualFileName)
    with open(variantsPerIndividualFileName, 'w') as f:
        json.dump(variantsPerIndividual, f, cls=NpEncoder)
    f.close()

    '''print('reading in variantsPerIndividual from ' + variantsPerIndividualFileName)
    with open(variantsPerIndividualFileName) as f:
        variantsPerIndividual = json.load(f, cls=NpDecoder)
    f.close()'''

    # TODO make this multi-threaded (cpu is pegged at 99% for this method)
    print('finding individuals per cooc')
    t = time.time()
    individualsPerBenignCooccurrence, individualsPerPathogenicCooccurrence, individualsPerVUSCooccurrence = \
        findIndividualsPerCooccurrence(variantsPerIndividual, benignPerGene, pathogenicPerGene, vusPerGene, ensemblRelease)
    elapsed_time = time.time() - t
    print('elapsed time in findIndividualsPerCooccurrence() ' + str(elapsed_time))



    # TODO make the program idempotent (save/read data frames and take cli arg to figure out)
    # p1 = P(VUS is benign and patient carries a pathogenic variant in trans)
    p1 = 0.5 * len(benignVariants) / count

    print('putting all the data together per vus')
    dataPerVus = calculateLikelihood(individualsPerPathogenicCooccurrence, individualsPerBenignCooccurrence,
                                     individualsPerVUSCooccurrence, p1)


    print('saving final VUS data  to ' + vusFinalDataFileName)
    with open(vusFinalDataFileName, 'w') as f:
        json.dump(dataPerVus, f, cls=NpEncoder)
    f.close()

    # now do the report: only report the "interesting" ones that pass the tests
    # data = [p1, p2, nk[vus][0], nk[vus][1], likelihoodRatios[vus], pathVarsPerVus[vus]]



def calculateLikelihood(pathCoocs, benCoocs, vusCoocs, p1):

    # vus coocs data: {(vus1, vus2):[individuals]}
    # "([10, 89624243, 'A', 'G'], [10, 89624304, 'C', 'T')]": ["0000057940", "0000057950"],
    # "([10, 89624304, 'C', 'T'], [10, 89624243, 'A', 'G')]": ["0000057940", "0000057950"],
    # "([10, 89624243, 'A', 'G'], [10, 89624298, 'C', 'T')]": ["0000057950"],
    # "([10, 89624304, 'C', 'T'], [10, 89624298, 'C', 'T')]": ["0000057950"],
    # "([10, 89624298, 'C', 'T'], [10, 89624243, 'A', 'G')]": ["0000057950"],
    # "([10, 89624298, 'C', 'T'], [10, 89624304, 'C', 'T')]": ["0000057950"]}

    #vusIndividuals = defaultdict(set)
    vusIndividuals = dict()
    for cooc in vusCoocs.keys():
        vus = repr(eval(cooc)[0])
        if vus not in vusIndividuals:
            vusIndividuals[vus] = set()
        for i in vusCoocs[cooc]:
            vusIndividuals[vus].add(i)

    for cooc in pathCoocs:
        vus = repr(eval(cooc)[0])
        if vus not in vusIndividuals:
            vusIndividuals[vus] = set()
        for i in pathCoocs[cooc]:
            vusIndividuals[vus].add(i)

    for cooc in benCoocs:
        vus = repr(eval(cooc)[0])
        if vus not in vusIndividuals:
            vusIndividuals[vus] = set()
        for i in benCoocs[cooc]:
            vusIndividuals[vus].add(i)

    # len of list associated with each vus is n (total number of times it was observed)
    vusList = list(vusIndividuals.keys())

    # path coocs data: {(vus, path):[individuals]}
    # "((10, 89624243, 'A', 'G'), (10, 89624245, 'GA', 'G'))": ["0000057940", "0000057950"],
    # "((10, 89624243, 'A', 'G'), (10, 89624248, 'A', 'G'))": ["0000057940"],
    # "((10, 89624304, 'C', 'T'), (10, 89624245, 'GA', 'G'))": ["0000057940", "0000057950"],
    # "((10, 89624304, 'C', 'T'), (10, 89624248, 'A', 'G'))": ["0000057940", "0000057960"],
    # "((10, 89624298, 'C', 'T'), (10, 89624245, 'GA', 'G'))": ["0000057950"]
    # }

    # now we calculate k (total number of times it co-occurred with a pathogenic variant)
    vusPathIndividuals = dict()
    for cooc in pathCoocs:
        vus = repr(eval(cooc)[0])
        if vus not in vusPathIndividuals:
            vusPathIndividuals[vus] = set()
        for i in pathCoocs[cooc]:
            vusPathIndividuals[vus].add(i)


    # per vus, calculate total number of times observed (n) and number of times observed with path variant (k)
    nk = dict()
    for vus in vusList:
        if vus not in vusIndividuals:
            n = 0
        else:
            n = len(vusIndividuals[vus])
        if vus not in vusPathIndividuals:
            k = 0
        else:
            k = len(vusPathIndividuals[vus])
        nk[vus] = (n, k)

    # now calculate log likelihood ratios!
    likelihoodRatios = dict()
    for vus in nk:
        n = nk[vus][0]
        k = nk[vus][1]
        denom = ((p1 ** k) * (1 - p1) ** (n - k))
        if k == 0:
            continue
        elif denom == 0:
            likelihoodRatios[vus] = sys.float_info.min
        else:
            likelihoodRatios[vus] = ((p2 ** k) * (1 - p2) ** (n - k)) / ((p1 ** k) * (1 - p1) ** (n - k))


    # find all the pathogenic variants this vus co-occurred with
    pathVarsPerVus = defaultdict(list)
    for cooc in pathCoocs:
        vus = repr(eval(cooc)[0])
        pathVarsPerVus[vus].append(eval(cooc)[1])

    # put it all together in a single dict
    dataPerVus = dict()
    for vus in likelihoodRatios:
        data = [p1, p2, nk[vus][0], nk[vus][1], likelihoodRatios[vus], pathVarsPerVus[vus]]
        dataPerVus[vus] = data

    return dataPerVus

def findVariantsPerGene(variantsPerChromosome, benignVariants, pathogenicVariants, unknownVariants, ensemblRelease, genesOfInterest):
    vusPerGene = defaultdict(list)
    pathogenicPerGene = defaultdict(list)
    benignPerGene = defaultdict(list)

    # TODO download ensembl db into docker container and build image from that (sometimes pyensembl install command fails)

    for chrom in variantsPerChromosome:
        x = variantsPerChromosome[chrom]
        for i in range(len(x)):
            chrom = int(chrom)
            try:
                pos = int(x.iloc[i][1])
                #pos = int(x.loc[i, 1])
                ref = str(x.iloc[i][3])
                #ref = str(x.loc[i, 3])
                alt = str(x.iloc[i][4])
                #alt = str(x.loc[i, 4])
                v = (chrom, pos, ref, alt)
            except Exception as e:
                print('exception at column ' + str(i))
                continue
            #print(v)


            # see if variant is pathogenic
            if v in pathogenicVariants:
                gene = getGeneForVariant(v, ensemblRelease)
                if gene is None or gene not in genesOfInterest:
                    continue
                else:
                    pathogenicPerGene[gene].append(v)

            # see if variant is benign
            elif v in benignVariants:
                gene = getGeneForVariant(v, ensemblRelease)
                if gene is None or gene not in genesOfInterest:
                    continue
                else:
                    benignPerGene[gene].append(v)

            # see if variant is vus? or just call it a vus b/c it's not a known benign or pathogenic var?
            #if v in unknownVariants:
            else:
                gene = getGeneForVariant(v, ensemblRelease)
                if gene is None or gene not in genesOfInterest:
                    continue
                else:
                    vusPerGene[gene].append(v)

    return benignPerGene, pathogenicPerGene, vusPerGene

def findVariantsPerChromosome(variants, chromosomes):
    variantsPerChromosome = defaultdict(dict)
    for chrom in variants['CHROM'].unique():
        if chrom in chromosomes:
            variantsPerChromosome[chrom] = variants.loc[variants['CHROM']==chrom]
    return variantsPerChromosome



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

    return chromsDF

def findVariantsInBRCA(fileName, classStrings, sigColName, hgVersion):
    brcaDF = pandas.read_csv(fileName, sep='\t', header=0, dtype=str)
    # Genomic_Coordinate_hg37
    # chr13:g.32972575:G>T
    pathVars = set()
    benignVars = set()
    vusVars = set()
    for i in range(len(brcaDF)):
        coord = brcaDF.loc[i, 'Genomic_Coordinate_hg' + hgVersion].split(':')
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

    return len(brcaDF), pathVars, benignVars, vusVars


def findVariantsPerIndividual(vcfDF, benignVariants, pathogenicVariants, unknownVariants, skipCols, nThreads, threadID):

    # find mutations
    # CHROM  POS             ID      REF     ALT     QUAL    FILTER  INFO            FORMAT  0000057940      0000057950
    # 10      89624243        .       A       G       .       .       AF=1.622e-05    GT      0/0             0/0
    # 0/0 => does not have variant on either strand (homozygous negative)
    # 0/1  => has variant on 1 strand (heterozygous positive)
    # 1/1 =>  has variant on both strands (homozygous positive)

    #variantsPerIndividual = defaultdict(lambda: defaultdict(list))
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

        includeList = ['1/0', '0/1', '1/1', '1|0', '0|1', '1|1']
        listOfVariantIndices = list()
        for x in includeList:
            listOfVariantIndices += list(np.where(vcfDF[individual] == x)[0])

        for i in listOfVariantIndices:
            try:
                '''record = vcfDF.iloc[i]
                varList = [record['CHROM'], int(record['POS']), record['REF'], record['ALT']]'''
                varList = (int(vcfDF.loc[i, 'CHROM']), int(vcfDF.loc[i, 'POS']), vcfDF.loc[i, 'REF'], vcfDF.loc[i, 'ALT'] )
                '''c = int(vcfDF.loc[i, 'CHROM'])
                p = int(vcfDF.loc[i, 'POS'])
                r = vcfDF.loc[i, 'REF']
                a = vcfDF.loc[i, 'ALT']
                varList = (c, p, r, a)'''
                #varList = [ int(vcfDF.loc[i, 'CHROM']), int(vcfDF.loc[i, 'POS']), vcfDF.loc[i, 'REF'], vcfDF.loc[i, 'ALT'] ]
                if varList in benignVariants:
                    variantsPerIndividual[individual]['benign'].append(varList)
                elif varList in pathogenicVariants:
                    variantsPerIndividual[individual]['pathogenic'].append(varList)
                # if not a known VUS, should we call it unknown here?
                #elif varList in unknownVariants:
                else:
                    variantsPerIndividual[individual]['vus'].append(varList)
            except Exception as e:
                print("exception for index " + str(i) + " of individual " + str(individual))
                print("exception: " + str(e))
                continue

    return variantsPerIndividual


def getGeneForVariant(variant, ensemblRelease):
    ensembl = pyensembl.EnsemblRelease(release=ensemblRelease)
    chrom = variant[0]
    if type(chrom) is str:
        # TODO this isn't graceful at all
        chrom = chrom.split('chr')[1]
    pos = variant[1]
    try:
        genes = ensembl.gene_names_at_locus(contig=int(chrom), position=int(pos))
        if len(genes) == 1:
            return genes[0]
        else:
            return None
    except Exception as e:
        print('exception: ' + str(e))
        return None

def findIndividualsPerCooccurrence(variantsPerIndividual, benPerGene, pathsPerGene, vusPerGene, ensemblRelease):
    individualsPerBenignCooccurrence = defaultdict(list)
    individualsPerPathogenicCooccurrence = defaultdict(list)
    individualsPerVUSCooccurrence = defaultdict(list)


    for individual in variantsPerIndividual:
        vusVarList = list(variantsPerIndividual[individual]['vus'])
        pathVarList = list(variantsPerIndividual[individual]['pathogenic'])
        benignVarList = list(variantsPerIndividual[individual]['benign'])
        for i in range(len(vusVarList)):
            vus_gene = getGeneForVariant(vusVarList[i], ensemblRelease)
            if vus_gene is None:
                continue
            for pathVar in pathVarList:
                path_gene = getGeneForVariant(pathVar, ensemblRelease)
                if vus_gene != path_gene or path_gene is None:
                    continue
                elif pathVar in pathsPerGene[path_gene] and vusVarList[i] in vusPerGene[vus_gene]:
                    individualsPerPathogenicCooccurrence[repr((vusVarList[i], pathVar))].append(individual)
            for benVar in benignVarList:
                ben_gene = getGeneForVariant(benVar, ensemblRelease)
                if vus_gene != ben_gene or ben_gene is None:
                    continue
                elif benVar in benPerGene[ben_gene] and vusVarList[i] in vusPerGene[vus_gene]:
                    individualsPerBenignCooccurrence[repr((vusVarList[i],benVar))].append(individual)
            for j in range(len(vusVarList)):
                if i == j:
                    continue
                vus_gene_2 = getGeneForVariant(vusVarList[j], ensemblRelease)
                if vus_gene != vus_gene_2 or vus_gene_2 is None:
                    continue
                elif vusVarList[j] in vusPerGene[vus_gene_2] and vusVarList[i] in vusPerGene[vus_gene]:
                    individualsPerVUSCooccurrence[repr((vusVarList[i], vusVarList[j]))].append(individual)

    return individualsPerBenignCooccurrence, individualsPerPathogenicCooccurrence, individualsPerVUSCooccurrence



if __name__ == "__main__":
    main()


