import allel
import sys
import numpy as np
from collections import defaultdict
import json
import pandas as pd


def main():
    vcfFileName = sys.argv[1]
    outputDir = sys.argv[2]

    vcf = allel.read_vcf(vcfFileName, fields=['samples', 'calldata/GT', 'variants/ALT', 'variants/CHROM',
            'variants/FILTER_PASS', 'variants/ID', 'variants/POS', 'variants/QUAL','variants/REF', 'variants/INFO'])

    genoArray = allel.GenotypeArray(vcf['calldata/GT'])

    alleleFrequency = genoArray.count_alleles().to_frequencies()

    expectedHeterozygosity = allel.heterozygosity_expected(alleleFrequency, ploidy=2)
    print('expected: ' + str(expectedHeterozygosity))

    observedHeterozygosity = allel.heterozygosity_observed(genoArray)
    print('observed: ' + str(observedHeterozygosity))

    inbreedingCoefficient = allel.inbreeding_coefficient(genoArray)
    print('ibc: ' + str(inbreedingCoefficient))

    diff = list()
    for i in range(len(expectedHeterozygosity)):
        diff.append(expectedHeterozygosity[i] - observedHeterozygosity[i])
    print('diff: ' + str(diff))

    numVariants = len(genoArray)
    numSamples = len(vcf['samples'])
    posArray = np.asarray([i for i in range(numVariants)])
    isAccessible = np.asarray([True for i in range(numVariants)])
    runsOfHomozygosity = dict()
    for i in range(numSamples):
        genoVector = genoArray[:,i]
        roh = allel.roh_mhmm(genoVector, posArray, is_accessible=isAccessible)
        if not roh[0].empty:
            runsOfHomozygosity[i] = {'individual': i, 'roh': allel.roh_mhmm(genoVector, posArray, is_accessible=isAccessible)}



    with open(outputDir + '/roh.json', 'w') as f:
        pd.DataFrame(runsOfHomozygosity).to_csv(f)
    f.close()

if __name__ == "__main__":
    main()