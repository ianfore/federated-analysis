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

clinvarVCFMetadataLines = 27
myVCFMetadataLines = 8
myVCFskipCols = 9
nThreads = cpu_count()
#nThreads=1
classStrings = { 'Pathogenic':[ 'Pathogenic' ], 'Benign':[ 'Benign', 'Likely benign' ],
                 'Unknown': [ 'Uncertain significance', '-']}
sigColName = 'Clinical_significance_ENIGMA'
brcaFileName = '/Users/jcasaletto/PycharmProjects/BIOBANK/federated-analysis/data/variants-test.tsv'
vcfFileName = '/Users/jcasaletto/PycharmProjects/BIOBANK/federated-analysis/data/bc100.vcf'
variantsPerIndividualFileName = '/Users/jcasaletto/PycharmProjects/BIOBANK/federated-analysis/data/variantsPerIndividual.json'
pathogenicCooccurrencesFileName = '/Users/jcasaletto/PycharmProjects/BIOBANK/federated-analysis/data/pathogenicCooccurrences.json'
benignCooccurrencesFileName = '/Users/jcasaletto/PycharmProjects/BIOBANK/federated-analysis/data/benignCooccurrences.json'
vusCooccurrencesFileName = '/Users/jcasaletto/PycharmProjects/BIOBANK/federated-analysis/data/vusCooccurrences.json'
vusFinalDataFileName = '/Users/jcasaletto/PycharmProjects/BIOBANK/federated-analysis/data/vusFinalData.json'
ensemblRelease=75

p2 = 0.0001


def main():
    if len(sys.argv) != 2:
        printUsage(sys.argv)
        sys.exit(1)
    if sys.argv[1] == '-p':
        produceOutputFiles()
    elif sys.argv[1] == '-c':
        consumeOutputFiles()
    elif sys.argv[1] == '-b':
        combo()
    else:
        printUsage(sys.argv)
        sys.exit(1)


def printUsage(args):
    sys.stderr.write("use -p to produce output files and -c to consume them")

def combo():
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
    benignPerGene, pathogenicPerGene, vusPerGene = findVariantsPerGene(variantsPerChromosome,
                                                                       benignVariants, pathogenicVariants,
                                                                       unknownVariants, ensemblRelease)
    elapsed_time = time.time() - t
    print('elapsed time in variantsPerGene() ' + str(elapsed_time))

    print('finding variants per individual')
    t = time.time()
    results = list()
    variantsPerIndividual = defaultdict(lambda: defaultdict(list))
    pool = ThreadPool(processes=nThreads)
    for i in range(nThreads):
        results.append(pool.apply_async(findVariantsPerIndividual, args=(vcfData, benignVariants, pathogenicVariants,
                                                                 unknownVariants, myVCFskipCols, nThreads, i)))
    for result in results:
        result.wait()
        variantsPerIndividual.update(result.get())
    elapsed_time = time.time() - t
    print('elapsed time in findVariantsPerIndividual() ' + str(elapsed_time))

    print('finding individuals per cooc')
    t = time.time()
    individualsPerBenignCooccurrence, individualsPerPathogenicCooccurrence, individualsPerVUSCooccurrence = \
        findIndividualsPerCooccurrence(variantsPerIndividual, benignPerGene, pathogenicPerGene, vusPerGene)
    elapsed_time = time.time() - t
    print('elapsed time in findIndividualsPerCooccurrence() ' + str(elapsed_time))

    # now you can you do the math!
    p1 = 0.5 * len(pathogenicVariants) / count

    dataPerVus = calculateLikelihood(individualsPerPathogenicCooccurrence, individualsPerBenignCooccurrence,
                                     individualsPerVUSCooccurrence, p1)
    print('saving final VUS data  to ' + vusFinalDataFileName)
    with open(vusFinalDataFileName, 'w') as f:
        json.dump(dataPerVus, f)
    f.close()
    print(dataPerVus)


def consumeOutputFiles():
    # do the math!
    print('doing the math!')

    print('reading BRCA data from ' + brcaFileName)
    t = time.time()
    count, pathogenicVariants, benignVariants, unknownVariants = \
        findVariantsInBRCA(brcaFileName, classStrings, sigColName)
    elapsed_time = time.time() - t
    print('elapsed time in findVariantsInBRCA() ' + str(elapsed_time))

    p1 = 0.5 * len(pathogenicVariants) / count

    with open(pathogenicCooccurrencesFileName) as f:
        individualsPerPathogenicCooccurrence = json.load(f)
    f.close()

    with open(benignCooccurrencesFileName) as f:
        individualsPerBenignCooccurrence = json.load(f)
    f.close()


    with open(vusCooccurrencesFileName) as f:
        individualsPerVUSCooccurrence = json.load(f)
    f.close()


    dataPerVus = calculateLikelihood(individualsPerPathogenicCooccurrence, individualsPerBenignCooccurrence, individualsPerVUSCooccurrence, p1)
    print('saving final VUS data  to ' + vusFinalDataFileName)
    with open(vusFinalDataFileName, 'w') as f:
        json.dump(dataPerVus, f)
    f.close()



def calculateLikelihood(pathCoocs, benCoocs, vusCoocs, p1):

    # vus coocs data: {(vus1, vus2):[individuals]}
    # "((10, 89624243, 'A', 'G'), (10, 89624304, 'C', 'T'))": ["0000057940", "0000057950"],
    # "((10, 89624304, 'C', 'T'), (10, 89624243, 'A', 'G'))": ["0000057940", "0000057950"],
    # "((10, 89624243, 'A', 'G'), (10, 89624298, 'C', 'T'))": ["0000057950"],
    # "((10, 89624304, 'C', 'T'), (10, 89624298, 'C', 'T'))": ["0000057950"],
    # "((10, 89624298, 'C', 'T'), (10, 89624243, 'A', 'G'))": ["0000057950"],
    # "((10, 89624298, 'C', 'T'), (10, 89624304, 'C', 'T'))": ["0000057950"]}

    vusIndividuals = defaultdict(set)
    for cooc in vusCoocs.keys():
        vus = eval(cooc)[0]
        for i in vusCoocs[cooc]:
            vusIndividuals[vus].add(i)

    for cooc in pathCoocs:
        vus = eval(cooc)[0]
        for i in pathCoocs[cooc]:
            vusIndividuals[vus].add(i)

    for cooc in benCoocs:
        vus = eval(cooc)[0]
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
    vusPathIndividuals = defaultdict(set)
    for cooc in pathCoocs:
        vus = eval(cooc)[0]
        for i in pathCoocs[cooc]:
            vusPathIndividuals[vus].add(i)


    # per vus, calculate total number of times observed (n) and number of times observed with path variant (k)
    nk = defaultdict(set)
    for vus in vusList:
        nk[vus] = (len(vusIndividuals[vus]), len(vusPathIndividuals[vus]))

    # now calculate log likelihood ratios!
    likelihoodRatios = defaultdict(float)
    for vus in nk:
        n = nk[vus][0]
        k = nk[vus][1]
        likelihoodRatios[vus] = ( (p2**k) * (1-p2)**(n-k)) / ((p1**k) * (1-p1)**(n-k))


    # find all the pathogenic variants this vus co-occurred with
    pathVarsPerVus = defaultdict(set)
    for cooc in pathCoocs:
        vus = eval(cooc)[0]
        pathVarsPerVus[vus].add(eval(cooc)[1])

    # put it all together in a single dict
    dataPerVus = dict()
    for vus in likelihoodRatios:
        data = (p1, p2, nk[vus][0], nk[vus][1], likelihoodRatios[vus], pathVarsPerVus[vus])
        dataPerVus[repr(vus)] = repr(data)

    return dataPerVus

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
    benignPerGene, pathogenicPerGene, vusPerGene = findVariantsPerGene(variantsPerChromosome,
                                                    benignVariants, pathogenicVariants, unknownVariants,ensemblRelease)
    elapsed_time = time.time() - t
    print('elapsed time in variantsPerGene() ' + str(elapsed_time))

    print('finding variants per individual')
    t = time.time()
    '''results = list()
    variantsPerIndividual = defaultdict(lambda: defaultdict(list))
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
        json.dump(variantsPerIndividual, f)
    f.close()'''

    with open(variantsPerIndividualFileName) as f:
        variantsPerIndividual = json.load(f)
    f.close()
    elapsed_time = time.time() - t
    print('elapsed time in findVariantsPerIndividual() ' + str(elapsed_time))




    print('finding individuals per cooc')
    t = time.time()
    individualsPerBenignCooccurrence, individualsPerPathogenicCooccurrence, individualsPerVUSCooccurrence = \
        findIndividualsPerCooccurrence(variantsPerIndividual, benignPerGene, pathogenicPerGene, vusPerGene)
    elapsed_time = time.time() - t
    print('elapsed time in findIndividualsPerCooccurrence() ' + str(elapsed_time))

    print('saving individuals per pathogenic cooc  to ' + pathogenicCooccurrencesFileName)
    with open(pathogenicCooccurrencesFileName, 'w') as f:
        json.dump(individualsPerPathogenicCooccurrence, f)
    f.close()
    print('saving individuals per benign cooc  to ' + benignCooccurrencesFileName)
    with open(benignCooccurrencesFileName, 'w') as f:
        json.dump(individualsPerBenignCooccurrence, f)
    f.close()
    print('saving individuals per vus cooc  to ' + vusCooccurrencesFileName)
    with open(vusCooccurrencesFileName, 'w') as f:
        json.dump(individualsPerVUSCooccurrence, f)
    f.close()

    # now you can you do the math!

def findVariantsPerGene(variantsPerChromosome, benignVariants, pathogenicVariants, unknownVariants, ensemblRelease):
    ensembl_db = pyensembl.database.Database(gtf_path='/Users/jcasaletto/Library/Caches/pyensembl/GRCh37.75.gtf.db')
    vusPerGene = defaultdict(list)
    pathogenicPerGene = defaultdict(list)
    benignPerGene = defaultdict(list)

    for chrom in variantsPerChromosome:
        x = variantsPerChromosome[chrom]
        for i in range(len(x)):
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
                if genes is None:
                    continue
                for gene in genes:
                    vusPerGene[gene].append((chrom, pos, ref, alt))

            # see if variant is pathogenic
            elif v in pathogenicVariants:
                genes = getGenesForVariant(v)
                if genes is None:
                    continue
                for gene in genes:
                    pathogenicPerGene[gene].append((chrom, pos, ref, alt))

            # see if variant is benign
            elif v in benignVariants:
                genes = getGenesForVariant(v)
                if genes is None:
                    continue
                for gene in genes:
                    benignPerGene[gene].append((chrom, pos, ref, alt))

    return benignPerGene, pathogenicPerGene, vusPerGene

def findVariantsPerChromosome(variants):
    variantsPerChromosome = defaultdict(dict)
    for chrom in variants['#CHROM'].unique():
        variantsPerChromosome[chrom] = variants.loc[variants['#CHROM']==chrom]
    return variantsPerChromosome



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


def findVariantsPerIndividual(vcfDF, benignVariants, pathogenicVariants, unknownVariants, skipCols, nThreads, threadID):

    # find mutations
    # #CHROM  POS             ID      REF     ALT     QUAL    FILTER  INFO            FORMAT  0000057940      0000057950
    # 10      89624243        .       A       G       .       .       AF=1.622e-05    GT      0/0             0/0
    # 0/0 => does not have variant on either strand (homozygous negative)
    # 0/1  => has variant on 1 strand (heterozygous positive)
    # 1/1 =>  has variant on both strands (homozygous positive)

    variantsPerIndividual = defaultdict(lambda: defaultdict(list))
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

        includeList = ['1/0', '0/1', '1/1']
        listOfVariantIndices = list()
        for x in includeList:
            listOfVariantIndices += list(np.where(vcfDF[individual] == x)[0])

        for variantIndex in listOfVariantIndices:
            try:
                record = vcfDF.iloc[variantIndex]
                varTuple = (record['#CHROM'], record['POS'], record['REF'], record['ALT'])
                if varTuple in benignVariants:
                    variantsPerIndividual[individual]['benign'].append(varTuple)
                elif varTuple in pathogenicVariants:
                    variantsPerIndividual[individual]['pathogenic'].append(varTuple)
                else:
                    variantsPerIndividual[individual]['vus'].append(varTuple)
            except Exception as e:
                print("exception for index " + str(variantIndex) + " of individual " + str(individual))
                print("exception: " + str(e))
                continue

    return variantsPerIndividual


def getGenesForVariant(variant):
    ensembl_db = pyensembl.database.Database(gtf_path='/Users/jcasaletto/Library/Caches/pyensembl/GRCh37.75.gtf.db')
    ensembl = pyensembl.EnsemblRelease(release=75)
    chrom = variant[0]
    pos = variant[1]
    try:
        genes = ensembl.gene_names_at_locus(contig=int(chrom), position=int(pos))
        if genes:
            return genes
        else:
            return None
    except Exception as e:
        return None

def findIndividualsPerCooccurrence(variantsPerIndividual, benPerGene, pathsPerGene, vusPerGene):
    individualsPerBenignCooccurrence = defaultdict(list)
    individualsPerPathogenicCooccurrence = defaultdict(list)
    individualsPerVUSCooccurrence = defaultdict(list)


    for individual in variantsPerIndividual:
        vusVarList = list(variantsPerIndividual[individual]['vus'])
        pathVarList = list(variantsPerIndividual[individual]['pathogenic'])
        benignVarList = list(variantsPerIndividual[individual]['benign'])
        #for i in range(len(vusVarList) - 1):
        for i in range(len(vusVarList)):
            vus_gene = getGenesForVariant(vusVarList[i])
            if vus_gene is None:
                continue
            for pathVar in pathVarList:
                path_gene = getGenesForVariant(pathVar)
                if vus_gene != path_gene or path_gene is None:
                    continue
                elif pathVar in pathsPerGene[path_gene[0]] and vusVarList[i] in vusPerGene[vus_gene[0]]:
                    individualsPerPathogenicCooccurrence[str((vusVarList[i], pathVar))].append(individual)
            for benVar in benignVarList:
                ben_gene = getGenesForVariant(benVar)
                if vus_gene != ben_gene or ben_gene is None:
                    continue
                elif benVar in benPerGene[ben_gene[0]] and vusVarList[i] in vusPerGene[vus_gene[0]]:
                    individualsPerBenignCooccurrence[str((vusVarList[i],benVar))].append(individual)
            #for j in range(i+1, len(vusVarList)):
            for j in range(len(vusVarList)):
                if i == j:
                    continue
                vus_gene_2 = getGenesForVariant(vusVarList[j])
                if vus_gene != vus_gene_2 or vus_gene_2 is None:
                    continue
                elif vusVarList[j] in vusPerGene[vus_gene_2[0]] and vusVarList[i] in vusPerGene[vus_gene[0]]:
                    individualsPerVUSCooccurrence[str((vusVarList[i], vusVarList[j]))].append(individual)

    return individualsPerBenignCooccurrence, individualsPerPathogenicCooccurrence, individualsPerVUSCooccurrence

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


