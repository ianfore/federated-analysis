import pandas
import itertools
from multiprocessing.pool import ThreadPool
from multiprocessing import cpu_count
import time
import sys
import json
from collections import defaultdict
import ast

clinvarVCFMetadataLines = 27
myVCFMetadataLines = 8
myVCFskipCols = 9
#nThreads = cpu_count()
nThreads=1
classStrings = { 'Pathogenic':[ 'Pathogenic' ], 'Benign':[ 'Benign', 'Likely benign' ],
                 'Unknown': [ 'Uncertain significance', '-']}
sigColName = 'Clinical_significance_ENIGMA'
brcaFileName = '/data/variants.tsv'
vcfFileName = '/data/BreastCancer.shuffle.vcf'
#vcfFileName = '/data/BreastCancer.shuffle-test.vcf'
#vcfFileName = '/data/bc-100.vcf'
variantsPerIndividualFileName = '/data/variantsPerIndividual.json'
cooccurrencesFileName = '/data/cooccurrences.json'
vusFileName = '/data/vus.json'
pathVarsFileName = '/data/pathogenicVariants.json'
vusToPathogenicVariantsFileName = '/data/vusToPathogenicVariants.json'
vusToBenignVariantsFileName = '/data/vusToBenignVariants.json'
vusToVusFileName ='/data/vusToVus.json'


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
    count, pathogenicVariants, benignVariants, unknownVariants = \
        findPathogenicVariantsInBRCA(brcaFileName, classStrings, sigColName)

    print('saving VUS to ' + vusFileName)
    with open(vusFileName, 'w') as f:
        json.dump(str(unknownVariants), f)
    f.close()

    print('saving pathogenic vars to ' + pathVarsFileName)
    with open(pathVarsFileName, 'w') as f:
        json.dump(str(pathogenicVariants), f)
    f.close()

    print('reading variant data from ' + vcfFileName)
    variants = readVariants(vcfFileName, myVCFMetadataLines)

    print('finding variants per individual')
    t = time.time()
    results = list()
    variantsPerIndividual = dict()
    pool = ThreadPool(processes=nThreads)
    for i in range(nThreads):
        results.append(pool.apply_async(findVariantsPerIndividual, args=(variants, myVCFskipCols, nThreads, i,)))
    for result in results:
        result.wait()
        variantsPerIndividual.update(result.get())
    elapsed_time = time.time() - t
    print('elapsed time is ' + str(elapsed_time))

    print('saving variants per individual to ' + variantsPerIndividualFileName)
    with open(variantsPerIndividualFileName, 'w') as f:
        json.dump(str(variantsPerIndividual), f)
    f.close()

    print('finding individuals with 1 or more pathogenic variant')
    pathogenicIndividuals = findIndividualsWithPathogenicVariant(variantsPerIndividual, pathogenicVariants)
    print('found ' + str(len(pathogenicIndividuals)) + ' individuals with a pathogenic variant')

    print('finding cooccurrences of pathogenic variant with any other variant')
    cooccurrences = findCooccurrences(pathogenicIndividuals)

    print('saving cooccurrences of pathogenic variant with any other variant')
    with open(cooccurrencesFileName, 'w') as f:
        json.dump(cooccurrences, f)
    f.close()

    print('mapping cooccurrences of pathogenic and benign variants with vus')
    vusToPathogenicVariants, vusToBenignVariants , vusToVus= mapVusToVariants(cooccurrences, unknownVariants, pathogenicVariants)
    print('saving cooccurrences of vus vars with pathogenic variants to ' + vusToPathogenicVariantsFileName)
    with open(vusToPathogenicVariantsFileName, 'w') as f:
        json.dump(vusToPathogenicVariants, f)
    f.close()
    print('saving cooccurrences of vus vars with benign variants to ' + vusToBenignVariantsFileName)
    with open(vusToBenignVariantsFileName, 'w') as f:
        json.dump(vusToBenignVariants, f)
    f.close()
    print('saving cooccurrences of vus vars with vus vars to ' + vusToVusFileName)
    with open(vusToVusFileName, 'w') as f:
        json.dump(vusToVus, f)
    f.close()

def consumeOutputFiles():
    # do the math!
    print('doing the math!')

def mapVusToVariants(cooccurrences, unknownVariants, pathogenicVariants):
    for c in cooccurrences:
        cooccurrences[c] = cooccurrences[c].replace('),', ');').replace('{', '').replace('}', '')

    vusToPath = defaultdict(list)
    vusToBenign = defaultdict(list)
    vusToVus = defaultdict(list)
    for c in cooccurrences.values():
        variants = c.split(';')
        for v in variants:
            # reconstruct variant
            v_split = v.strip().replace('(', '').replace(')', '').split(',')
            chrom = v_split[0].strip()
            pos = v_split[1].strip()
            ref = ast.literal_eval(v_split[2].strip())
            alt = ast.literal_eval(v_split[3].strip())
            tup = (int(chrom), int(pos), ref, alt)

            # look up variant to see if it's vus
            if tup in unknownVariants:
                # now generate list of pathogenic vars that cooccur with vus
                for cv in variants:
                    if cv in pathogenicVariants:
                        vusToPath[tup].append(cv)
                    elif cv not in unknownVariants:
                        vusToBenign[tup].append(cv)
                    else:
                        vusToVus[tup].append(cv)
    return vusToPath, vusToBenign, vusToVus

def readVariants(fileName, numMetaDataLines):
    # #CHROM  POS             ID      REF     ALT     QUAL    FILTER  INFO            FORMAT  0000057940      0000057950
    # 10      89624243        .       A       G       .       .       AF=1.622e-05    GT      0/0             0/0
    # 0/0 => does not have variant on either strand (homozygous negative)
    # 0/1  => has variant on 1 strand (heterozygous positive)
    # 1/1 =>  has variant on both strands (homozygous positive)
    df = pandas.read_csv(fileName, sep='\t', skiprows=numMetaDataLines, dtype={'#CHROM':int, 'POS':int}, header=0)
    # this creates a bug: df = df[df.apply(lambda r: r.str.contains('1/1').any() or r.str.contains('0/1').any(), axis=1)]
    return df

def findPathogenicVariantsInBRCA(fileName, classStrings, sigColName):
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


def findVariantsPerIndividual(df, skipCols, nThreads, threadID):

    # find mutations
    # #CHROM  POS             ID      REF     ALT     QUAL    FILTER  INFO            FORMAT  0000057940      0000057950
    # 10      89624243        .       A       G       .       .       AF=1.622e-05    GT      0/0             0/0
    # 0/0 => does not have variant on either strand (homozygous negative)
    # 0/1  => has variant on 1 strand (heterozygous positive)
    # 1/1 =>  has variant on both strands (homozygous positive)

    variantsPerIndividual = dict()
    numColumns = len(df.columns)
    numIndividuals = len(df.columns) - skipCols
    partitionSize = int(numIndividuals / nThreads)
    start = threadID * partitionSize + skipCols
    if threadID == nThreads - 1:
        end = numColumns
    else:
        end = skipCols + start + partitionSize

    # get list of individuals
    individuals = df.columns[start:end]

    # iterate through columns (not rows! iterows() takes 20x longer b/c pandas are stored column-major)
    for individual in individuals:
        variantsPerIndividual[individual] = set()

        # transform column values from strings to ints
        #from sklearn.preprocessing import LabelEncoder
        #enc = LabelEncoder()
        #enc.fit(df[individual])
        #df[individual] = enc.transform(df[individual])
        #listOfVariantIndices = list(df[(df[individual] != 0)].index)

        listOfVariantIndices = list(df[(df[individual] != '0/0')].index)

        for variantIndex in listOfVariantIndices:
            try:
                record = df.iloc[variantIndex]
                varTuple = (record['#CHROM'], record['POS'], record['REF'], record['ALT'])
                variantsPerIndividual[individual].add(varTuple)
            except Exception as e:
                print("exception for index " + str(variantIndex) + " of individual " + str(individual))
                print("exception: " + str(e))
                continue

    return variantsPerIndividual

def findIndividualsWithPathogenicVariant(variantsPerIndividual, pathogenicVars):
    # variantsPerIndividua[userID] = [(chrom, pos, ref, alt, qual), ... ]

    # 10748 = [(10, 89624243, 'A', 'G', '.'), (13, 32910721, 'T', 'C'), (13, 32910842, 'A', 'G')]
    # 27089 = [(10, 89624245, 'GA', 'G', '.'), (13, 32906729, 'A', 'C'), (13, 32910721, 'T', 'C')]
    # => both have (13, 32910721, 'T', 'C')

    individualsWithPathogenicVariant = dict()

    for individual in variantsPerIndividual.keys():
        for variant in variantsPerIndividual[individual]:
            if variant in pathogenicVars:
                individualsWithPathogenicVariant[repr(individual) + '_' + str(variant)] = variantsPerIndividual[individual]
                break
    return individualsWithPathogenicVariant

def findCooccurrences(patientsWithPathogenicVars):
    cooccurrences = dict()
    for a, b in itertools.combinations(patientsWithPathogenicVars.keys(), 2):
        intersection = patientsWithPathogenicVars[a].intersection(patientsWithPathogenicVars[b])
        if len(intersection) > 0:
            d1 = {eval(a):patientsWithPathogenicVars[a]}
            d2 = {eval(b):patientsWithPathogenicVars[b]}
            key = (d1, d2)
            cooccurrences[str(key)] = str(intersection)
    return cooccurrences


if __name__ == "__main__":
    main()


