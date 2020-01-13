import pandas
import itertools
from multiprocessing.pool import ThreadPool
from multiprocessing import cpu_count
import time

clinvarVCFMetadataLines = 27
myVCFMetadataLines = 8
myVCFskipCols = 9
nThreads = cpu_count()
#nThreads=1
classStrings = { 'Pathogenic':[ 'Pathogenic' ], 'Benign':[ 'Benign', 'Likely benign' ],
                 'Unknown': [ 'Uncertain significance'], 'Unclassified': [ '-']}
sigColName = 'Clinical_significance_ENIGMA'
brcaFileName = '/data/variants.tsv'
vcfFileName = '/data/BreastCancer.shuffle.vcf'
#vcfFileName = '/data/BreastCancer.shuffle-test.vcf'
#vcfFileName = '/data/bc-100.vcf'
variantsPerIndividualFileName = '/data/variantsPerIndividual.txt'

def main():

    print('reading variants from ' + brcaFileName)
    count, pathogenicVariants, benignVariants, unknownVariants, unclassifiedVariants = findPathogenicVariantsInBRCA(brcaFileName, classStrings, sigColName)
    print('of ' + str(count) + ' variants, found ' + str(len(pathogenicVariants)) + ' pathogenic variants in brca')
    print('of ' + str(count) + ' variants, found ' + str(len(benignVariants)) + ' benign variants in brca')
    print('of ' + str(count) + ' variants, found ' + str(len(unknownVariants)) + ' unknown variants in brca')
    print('of ' + str(count) + ' variants, found ' + str(len(unclassifiedVariants)) + ' unclassified variants in brca')


    print('reading variants from ' + vcfFileName)
    variants = readVariants(vcfFileName, myVCFMetadataLines)
    print('found ' + str(len(variants)) + ' variants')
    print('found ' + str(len(variants.columns) - myVCFskipCols) + ' individuals')

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

    print('saving dictionary to ' + variantsPerIndividualFileName)
    with open(variantsPerIndividualFileName, 'w') as file:
        file.write(str(variantsPerIndividual))

    # then to read it back in
    '''import ast
    with open(variantsPerIndividualFileName, 'r') as f:
        my_set = ast.literal_eval(f.read())'''

    elapsed_time = time.time() - t
    print('elapsed time is ' + str(elapsed_time))

    print('finding individuals with 1 or more pathogenic variant')
    pathogenicIndividuals = findIndividualsWithPathogenicVariant(variantsPerIndividual, pathogenicVariants)
    print('found ' + str(len(pathogenicIndividuals)) + ' individuals with a pathogenic variant')

    print('finding cooccurrences of pathogenic variant with any other variant')
    cooccurrences = findCooccurrences(pathogenicIndividuals)
    for k in cooccurrences.keys():
        print('intersection of ' + str(k) + ' is ' + str(cooccurrences[k]))

    # TODO find "interesting" cooccurrences

def readVariants(fileName, numMetaDataLines):
    # #CHROM  POS             ID      REF     ALT     QUAL    FILTER  INFO            FORMAT  0000057940      0000057950
    # 10      89624243        .       A       G       .       .       AF=1.622e-05    GT      0/0             0/0
    # 0/0 => does not have variant on either strand (homozygous negative)
    # 0/1  => has variant on 1 strand (heterozygous positive)
    # 1/1 =>  has variant on both strands (homozygous positive)
    df = pandas.read_csv(fileName, sep='\t', skiprows=numMetaDataLines, dtype={'#CHROM':int, 'POS':int})
    df = df[df.apply(lambda r: r.str.contains('1/1').any() or r.str.contains('0/1').any(), axis=1)]
    return df



def findPathogenicVariantsInBRCA(fileName, classStrings, sigColName):
    brcaDF = pandas.read_csv(fileName, sep='\t', header=0, dtype=str)
    # Genomic_Coordinate_hg37
    # chr13:g.32972575:G>T
    pathVars = list()
    vusVars = list()
    benignVars = list()
    unclassifiedVars = list()

    for index, row in brcaDF.iterrows():
        coord = row['Genomic_Coordinate_hg37']
        coord = coord.split(':')
        chr = int(coord[0].split('chr')[1])
        pos = int(coord[1].split('.')[1])
        ref = coord[2].split('>')[0]
        alt = coord[2].split('>')[1]
        tuple = (chr, pos, ref, alt)

        if str(row[sigColName]) in classStrings['Pathogenic']:
            #pathVars.append((int(row['Chr']), int(row['Pos']), row['Ref'], row['Alt']))
            pathVars.append(tuple)
        elif str(row[sigColName]) in classStrings['Benign']:
            #benignVars.append((int(row['Chr']), int(row['Pos']), row['Ref'], row['Alt']))
            benignVars.append(tuple)
        elif str(row[sigColName]) in classStrings['Unknown']:
            #vusVars.append((int(row['Chr']), int(row['Pos']), row['Ref'], row['Alt']))
            vusVars.append(tuple)
        elif str(row[sigColName]) in classStrings['Unclassified']:
            #unclassifiedVars.append((int(row['Chr']), int(row['Pos']), row['Ref'], row['Alt']))
            unclassifiedVars.append(tuple)

    return len(brcaDF), pathVars, benignVars, vusVars, unclassifiedVars


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

    #
    for individual in individuals:
        #variantsPerIndividual[individual] = set(df[df[individual] != '0/0'].index)
        variantsPerIndividual[individual] = set()
        #isVariant = df[individual] != '0/0'
        #variantsPerIndividualDF = df[isVariant][['#CHROM', 'POS', 'REF', 'ALT']]


        listOfVariantIndices = list(df[(df[individual] != '0/0')].index)
        for variantIndex in listOfVariantIndices:
            try:
                record = df.iloc[variantIndex]
                varTuple = (record['#CHROM'], record['POS'], record['REF'], record['ALT'])
                variantsPerIndividual[individual].add(varTuple)
            except:
                print('exception: ' + str(variantIndex) + ' from ' + str(listOfVariantIndices))
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
                individualsWithPathogenicVariant[repr(individual)] = variantsPerIndividual[individual]
                break
    return individualsWithPathogenicVariant

def findCooccurrences(patientsWithPathogenicVars):
    cooccurrences = dict()
    for a, b in itertools.combinations(patientsWithPathogenicVars.values(), 2):
        #print('comparing ' + repr(a) + ' and ' + repr(b))
        intersection = a.intersection(b)
        if len(intersection) > 0 and (repr(a), repr(b)) not in cooccurrences and (repr(b), repr(a)) not in cooccurrences:
            print('found match b/w ' + repr(a) + ' and ' + repr(b))
            cooccurrences[(repr(a), repr(b))] = intersection
    return cooccurrences


if __name__ == "__main__":
    main()


