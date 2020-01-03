import pandas
import itertools
from multiprocessing.pool import ThreadPool
import time

clinvarVCFMetadataLines = 27
myVCFMetadataLines = 8
myVCFskipCols = 9
nThreads = 4

def main():


    print('finding pathogenic variants from brca')
    count, pathogenicVariants = findPathogenicVariantsInBRCA('/data/variants.tsv', 'Pathogenic', 'Clinical_significance_ENIGMA')
    print('of ' + str(count) + ' variants, found ' + str(len(pathogenicVariants)) + ' pathogenic variants in brca')

    print('pre-processing data: removing variants with all 0/0')
    count, nonNullVariants = readOccurringVariants('/data/BreastCancer.shuffle.vcf', myVCFMetadataLines)
    print('of ' + str(count) + ' variants, found ' + str(len(nonNullVariants)) + ' that are occurring')

    print('pre-processing data: removing individuals with all 0/0')
    count, nonNullIndividuals = removeNonmutatedIndividuals(nonNullVariants, myVCFskipCols)
    print('of ' + str(count) + ' individuals, found ' + str(len(nonNullIndividuals.columns) - myVCFskipCols) + ' that are mutated')

    print('indexing individuals with variants')
    t = time.process_time()
    results = list()
    variantsPerIndividual = dict()
    with ThreadPool(processes=nThreads) as pool:
        for i in range(nThreads):
            results.append(pool.apply_async(findVariantsPerIndividual, args=(nonNullIndividuals, myVCFskipCols, nThreads, i,)))
        for result in results:
            result.wait()
            variantsPerIndividual.update(result.get())
    elapsed_time = time.process_time() - t
    print('elapsed time is ' + str(elapsed_time))

    print('of ' + str(len(nonNullIndividuals.columns) - myVCFskipCols) + ' individuals, found ' + str(len(variantsPerIndividual)) + ' with some type of variant')

    print('finding individuals with pathogenic variant')
    pathogenicIndividuals = findIndividualsWithPathogenicVariant(variantsPerIndividual, pathogenicVariants)
    print('of ' + str(len(variantsPerIndividual)) + ' individuals with variants, found ' + str(len(pathogenicIndividuals)) + ' with a pathogenic variant')

    print('finding cooccurrences of pathogenic variant with any other variant')
    cooccurrences = findCooccurrences(pathogenicIndividuals)
    for k in cooccurrences.keys():
        print('intersection of ' + str(k) + ' is ' + str(cooccurrences[k]))

def readOccurringVariants(fileName, numMetaDataLines):
    # #CHROM  POS             ID      REF     ALT     QUAL    FILTER  INFO            FORMAT  0000057940      0000057950
    # 10      89624243        .       A       G       .       .       AF=1.622e-05    GT      0/0             0/0
    # 0/0 => does not have variant on either strand (homozygous negative)
    # 0/1  => has variant on 1 strand (heterozygous positive)
    # 1/1 =>  has variant on both strands (homozygous positive)
    df = pandas.read_csv(fileName, sep='\t', skiprows=numMetaDataLines, dtype={'#CHROM':int, 'POS':int})
    count = len(df)
    df = df[df.apply(lambda r: r.str.contains('1/1').any() or r.str.contains('0/1').any(), axis=1)]
    return count, df

def removeNonmutatedIndividuals(df, skipcols):
    # #CHROM  POS             ID      REF     ALT     QUAL    FILTER  INFO            FORMAT  0000057940      0000057950
    # 10      89624243        .       A       G       .       .       AF=1.622e-05    GT      0/0             0/0
    # 0/0 => does not have variant on either strand (homozygous negative)
    # 0/1  => has variant on 1 strand (heterozygous positive)
    # 1/1 =>  has variant on both strands (homozygous positive)
    count = len(df.columns) - skipcols
    df = df.loc[:, (df != '0/0').any(axis=0)]
    return count, df

def findPathogenicVariantsInBRCA(fileName, significanceString, sigColName):
    brcaDF = pandas.read_csv(fileName, sep='\t', header=0, dtype=str)
    pathVars = list()

    for index, row in brcaDF.iterrows():
        if significanceString.lower() == str(row[sigColName]).lower():
            pathVars.append((int(row['Chr']), int(row['Pos']), row['Ref'], row['Alt']))
    return len(brcaDF), pathVars


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
    start = threadID * partitionSize
    if threadID == nThreads - 1:
        end = numColumns
    else:
        end = skipCols + start + partitionSize

    for index, record in df.iterrows():
        for individual in range(skipCols + start, end):
            id = list(df.columns)[individual]
            if record[individual] == '0/1' or record[individual] == '1/1':
                if id not in variantsPerIndividual:
                    variantsPerIndividual[id] = set()
                varTuple = (record['#CHROM'], record['POS'], record['REF'], record['ALT'])
                variantsPerIndividual[id].add(varTuple)
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
        if len(intersection) > 0:
            cooccurrences[(repr(a), repr(b))] = intersection
    return cooccurrences


if __name__ == "__main__":
    main()


