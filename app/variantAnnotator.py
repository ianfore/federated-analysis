import allel
import sys
import numpy as np
from collections import defaultdict
import json
import pandas as pd


def main():
    vcfFileName = sys.argv[1]
    rohOutputFileName = sys.argv[2]

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

    n = len(vcf['samples'])
    posArray = np.asarray([i for i in range(n)])
    isAccessible = np.asarray([True for i in range(n)])
    runsOfHomozygosity = defaultdict()
    for i in range(n):
        genoVector = genoArray[i]
        runsOfHomozygosity[i] = allel.roh_mhmm(genoVector, posArray, is_accessible=isAccessible)


    with open(rohOutputFileName, 'w') as f:
        f.write(str(runsOfHomozygosity))
    f.close()

if __name__ == "__main__":
    main()