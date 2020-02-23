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

clinvarVCFMetadataLines = 27
myVCFMetadataLines = 8
myVCFskipCols = 9
nThreads = cpu_count()
#nThreads=1
classStrings = { 'Pathogenic':[ 'Pathogenic' ], 'Benign':[ 'Benign', 'Likely benign' ],
                 'Unknown': [ 'Uncertain significance', '-']}
sigColName = 'Clinical_significance_ENIGMA'
brcaFileName = '/Users/jcasaletto/PycharmProjects/BIOBANK/federated-analysis/data/brca-variants.tsv'
vcfFileName = '/Users/jcasaletto/PycharmProjects/BIOBANK/federated-analysis/data/BreastCancer.shuffle.vcf'
variantsPerIndividualFileName = '/Users/jcasaletto/PycharmProjects/BIOBANK/federated-analysis/data/variantsPerIndividual.json'
cooccurrencesFileName = '/Users/jcasaletto/PycharmProjects/BIOBANK/federated-analysis/data/cooccurrences.json'
vusFileName = '/Users/jcasaletto/PycharmProjects/BIOBANK/federated-analysis/data/vus.json'
pathVarsFileName = '/Users/jcasaletto/PycharmProjects/BIOBANK/federated-analysis/data/pathogenicVariants.json'
benignVarsFileName = '/Users/jcasaletto/PycharmProjects/BIOBANK/federated-analysis/data/benignVariants.json'
vusToPathogenicVariantsFileName = '/Users/jcasaletto/PycharmProjects/BIOBANK/federated-analysis/data/vusToPathogenicVariants.json'
vusToBenignVariantsFileName = '/Users/jcasaletto/PycharmProjects/BIOBANK/federated-analysis/data/vusToBenignVariants.json'
vusToVusFileName ='/Users/jcasaletto/PycharmProjects/BIOBANK/federated-analysis/data/vusToVus.json'
ensemblRelease=75


def main():
    if len(sys.argv) != 2:
        printUsage(sys.argv)
        sys.exit(1)
    if sys.argv[1] == '-p':
        produceOutputFiles()
    elif sys.argv[1] == '-c':
        consumeOutputFiles()
    else:
        printUsage(sys.argv)
        sys.exit(1)


def printUsage(args):
    sys.stderr.write("use -p to produce output files and -c to consume them")

def produceOutputFiles():
    print("producing output files!")


    print('reading BRCA data from ' + brcaFileName)
    t = time.time()
    count, pathogenicVariants, benignVariants, unknownVariants = \
        findVariantsInBRCA(brcaFileName, classStrings, sigColName)
    elapsed_time = time.time() - t
    print('elapsed time in findVariantsInBRCA() ' + str(elapsed_time))

    print('reading VCF data from ' + vcfFileName)
    t = time.time()
    vcfData = readVCFFile(vcfFileName, myVCFMetadataLines)
    elapsed_time = time.time() - t
    print('elapsed time in readVCFFile() ' + str(elapsed_time))

    print('slicing up VCF data by chromosome')
    t = time.time()
    variantsPerChromosome = findVariantsPerChromosome(vcfData)
    elapsed_time = time.time() - t
    print('elapsed time in variantsPerChromosome() ' + str(elapsed_time))

    print('slicing up VCF data by gene')
    t = time.time()
    vusPerGene, pathogenicPerGene = findVariantsPerGene(variantsPerChromosome, unknownVariants, pathogenicVariants, ensemblRelease)
    elapsed_time = time.time() - t
    print('elapsed time in variantsPerGene() ' + str(elapsed_time))

    print('finding non-benign variants per individual')
    t = time.time()
    results = list()
    nonBenignVariantsPerIndividual = dict()
    pool = ThreadPool(processes=nThreads)
    for i in range(nThreads):
        results.append(pool.apply_async(findNonBenignVariantsPerIndividual, args=(vcfData, benignVariants, myVCFskipCols, nThreads, i)))
    for result in results:
        result.wait()
        nonBenignVariantsPerIndividual.update(result.get())
    elapsed_time = time.time() - t
    print('elapsed time in findVariantsPerIndividual() ' + str(elapsed_time))

    print('finding vus-path co-occurring variants per individual')
    t = time.time()
    # for each person, find list of pairs of vus and path that co-occur on same gene
    coocsPerIndividual = findCooccurrencesPerIndividual(nonBenignVariantsPerIndividual, pathogenicVariants, unknownVariants)
    elapsed_time = time.time() - t
    print('elapsed time in findVariantsPerIndividual() ' + str(elapsed_time))

    for i in coocsPerIndividual:
        print('individual: ' + str(i))
        print('coocs = ' + str(coocsPerIndividual[i]))

    # then for each co-occurrence pair, add individual to list
    # length of list is number of co-coccurrences for vus-path pair
    # now you can you do the math!

def findVariantsPerGene(variantsPerChromosome, unknownVariants, pathogenicVariants, ensemblRelease):
    ensembl_db = pyensembl.database.Database(gtf_path='/Users/jcasaletto/Library/Caches/pyensembl/GRCh37.75.gtf.db')
    vusPerGene = defaultdict(list)
    pathogenicPerGene = defaultdict(list)
    for chrom in variantsPerChromosome:
        x = variantsPerChromosome[chrom]
        for i in range(len(x.iloc[1])):
            chrom = int(chrom)
            try:
                pos = int(x.iloc[i][1])
                ref = str(x.iloc[i][3])
                alt = str(x.iloc[i][4])
                v = (chrom, pos, ref, alt)
            except Exception as e:
                print('exception at column ' + str(i))
                #print('exception getting variant out of x.iloc[' + str(i) + ']: ' + str(x.iloc[i]))
                continue
            #print(v)
            # see if variant is vus
            if v in unknownVariants:
                genes = getGenesForVariant(v)
                for gene in genes:
                    vusPerGene[gene].append((chrom, pos, ref, alt))

            # see if variant is pathogenic
            elif v in pathogenicVariants:
                genes = getGenesForVariant(v)
                for gene in genes:
                    pathogenicPerGene[gene].append((chrom, pos, ref, alt))

    return vusPerGene, pathogenicPerGene

def findVariantsPerChromosome(variants):
    variantsPerChromosome = defaultdict(dict)
    for chrom in variants['#CHROM'].unique():
        variantsPerChromosome[chrom] = variants.loc[variants['#CHROM']==chrom]
    return variantsPerChromosome

def consumeOutputFiles():
    # do the math!
    print('doing the math!')

def readVCFFile(fileName, numMetaDataLines):
    # #CHROM  POS             ID      REF     ALT     QUAL    FILTER  INFO            FORMAT  0000057940      0000057950
    # 10      89624243        .       A       G       .       .       AF=1.622e-05    GT      0/0             0/0
    # 0/0 => does not have variant on either strand (homozygous negative)
    # 0/1  => has variant on 1 strand (heterozygous positive)
    # 1/1 =>  has variant on both strands (homozygous positive)
    df = pandas.read_csv(fileName, sep='\t', skiprows=numMetaDataLines, dtype={'#CHROM':int, 'POS':int}, header=0)
    # this creates a bug: df = df[df.apply(lambda r: r.str.contains('1/1').any() or r.str.contains('0/1').any(), axis=1)]
    return df

def findVariantsInBRCA(fileName, classStrings, sigColName):
    brcaDF = pandas.read_csv(fileName, sep='\t', header=0, dtype=str)
    # Genomic_Coordinate_hg37
    # chr13:g.32972575:G>T
    pathVars = list()
    benignVars = list()
    vusVars = list()

    for index, row in brcaDF.iterrows():
        coord = row['Genomic_Coordinate_hg37']
        coord = coord.split(':')
        chrom = int(coord[0].split('chr')[1])
        pos = int(coord[1].split('g.')[1])
        ref, alt = coord[2].split('>')
        tup = (chrom, pos, ref, alt)

        if str(row[sigColName]) in classStrings['Pathogenic']:
            pathVars.append(tup)
        elif str(row[sigColName]) in classStrings['Benign']:
            benignVars.append(tup)
        elif str(row[sigColName]) in classStrings['Unknown']:
            vusVars.append(tup)

    return len(brcaDF), pathVars, benignVars, vusVars


def findNonBenignVariantsPerIndividual(vcfDF, benignVariants, skipCols, nThreads, threadID):

    # find mutations
    # #CHROM  POS             ID      REF     ALT     QUAL    FILTER  INFO            FORMAT  0000057940      0000057950
    # 10      89624243        .       A       G       .       .       AF=1.622e-05    GT      0/0             0/0
    # 0/0 => does not have variant on either strand (homozygous negative)
    # 0/1  => has variant on 1 strand (heterozygous positive)
    # 1/1 =>  has variant on both strands (homozygous positive)

    variantsPerIndividual = defaultdict(list)
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

        listOfVariantIndices = list(vcfDF[(vcfDF[individual] != '0/0')].index)

        for variantIndex in listOfVariantIndices:
            try:
                record = vcfDF.iloc[variantIndex]
                varTuple = (record['#CHROM'], record['POS'], record['REF'], record['ALT'])
                if varTuple not in benignVariants:
                    variantsPerIndividual[individual].append(varTuple)
            except Exception as e:
                print("exception for index " + str(variantIndex) + " of individual " + str(individual))
                print("exception: " + str(e))
                continue

    return variantsPerIndividual


def getGenesForVariant(variant):
    ensembl = pyensembl.EnsemblRelease(release=75)
    chrom = variant[0]
    pos = variant[1]
    try:
        return ensembl.gene_names_at_locus(contig=chrom, position=pos)
    except Exception as e:
        return None

def findCooccurrencesPerIndividual(nonBenignVariantsPerIndividual, pathsPerGene, vusPerGene):
    coocsPerIndividual = defaultdict(list)

    for i in nonBenignVariantsPerIndividual:
        variantList = nonBenignVariantsPerIndividual[i]
        for j in range(len(variantList)-1):
            for k in range(j+1, len(variantList)):
                gene_j = getGenesForVariant(variantList[j])
                gene_k = getGenesForVariant(variantList[k])
                if gene_j != gene_k or gene_j is None:
                    continue
                elif variantList[j] in pathsPerGene[gene_j] and variantList[k] in vusPerGene[gene_j]:
                    coocsPerIndividual[i].append((variantList[j], variantList[k]))

                elif variantList[k] in pathsPerGene[gene_k] and variantList[j] in vusPerGene[gene_k]:
                    coocsPerIndividual[i].append((variantList[k], variantList[j]))

    return coocsPerIndividual

def findCooccurrences(patientsWithPathogenicVars):
    cooccurrences = dict()
    for a, b in itertools.combinations(patientsWithPathogenicVars.keys(), 2):
        intersection = patientsWithPathogenicVars[a].intersection(patientsWithPathogenicVars[b])
        if len(intersection) > 0:
            #d1 = {eval(a):patientsWithPathogenicVars[a]}
            #d2 = {eval(b):patientsWithPathogenicVars[b]}
            d1 = {a: patientsWithPathogenicVars[a]}
            d2 = {b: patientsWithPathogenicVars[b]}
            key = (d1, d2)
            cooccurrences[str(key)] = str(intersection)
    return cooccurrences


if __name__ == "__main__":
    main()


