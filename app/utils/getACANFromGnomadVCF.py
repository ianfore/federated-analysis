import pandas as pd
import logging
import argparse
import json
import matplotlib.pyplot as plt
import numpy
from scipy.stats import ks_2samp
import math

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.WARNING)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', help='vcf input file')
    parser.add_argument('-j', '--json', help='json output file')
    parser.add_argument('-g', '--graph', help='graph output file')
    options = parser.parse_args()
    return options


def main():
    vcfFileName = parse_args().vcf
    outputFileName = parse_args().json
    graphFileName = parse_args().graph


    logger.info('finding variants from ' + vcfFileName)
    vcfDF = pd.read_csv(vcfFileName, delimiter='\t', header=0, dtype=str)
    nontopmedKeys = ['AF-non_topmed-afr', 'AF-non_topmed-amr', 'AF-non_topmed-nfe',
                     'AF-non_topmed-fin', 'AF-non_topmed-eas', 'AF-non_topmed-sas']

    nonTopmedVariants = getNontopmedAlleleFreqs(vcfDF, nontopmedKeys)

    topmedVariants = getTopmedAlleleFreqs(vcfDF)

    logger.info('saving ntm allele freqs')
    with open(outputFileName, 'w') as f:
        json.dump(nonTopmedVariants, f)
    f.close()

    logger.info('saving tm allele freqs')
    with open(outputFileName + 'justtm', 'w') as f:
        json.dump(topmedVariants, f)
    f.close()

    '''logger.info('plotting allele freqs')
    topmedDict, nontopmedDict = createDicts(nontopmedVariants, topmedKeys, nontopmedKeys)
    plotDists(topmedDict, nontopmedDict, topmedKeys, nontopmedKeys, graphFileName)'''

def getTopmedAlleleFreqs(vcfDF):
    ethnicitiesList = ['afr', 'amr', 'eas', 'fin', 'nfe', 'sas']
    # AC-non_topmed-afr, AN-non_topmed-afr, AC-afr, AN-afr
    keys = list()
    for e in ethnicitiesList:
        ntmACKey = 'AC-non_topmed-' + e
        keys.append(ntmACKey)
        ntmANKey = 'AN-non_topmed-' + e
        keys.append(ntmANKey)
        tmACKey = 'AC-' + e
        keys.append(tmACKey)
        tmANKey = 'AN-' + e
        keys.append(tmANKey)

    print(keys)

    # AF-non_topmed-afr AF-afr
    justTopmedVariantsDict = dict()

    # get AC and AN for each ethnicity, then calculate AF
    for i in range(len(vcfDF)):
        fields = vcfDF.iloc[i]['INFO'].split(';')
        chrom = vcfDF.iloc[i]['#CHROM'].split('chr')[1]
        pos = vcfDF.iloc[i]['POS']
        ref = vcfDF.iloc[i]['REF']
        alt = vcfDF.iloc[i]['ALT']
        mykey = "(" + chrom + ", " + pos + ", '" + ref + "', '" + alt + "')"
        justTopmedVariantsDict[mykey] = dict()
        for field in fields:
            if '=' in field:
                key = field.split('=')[0]
                if key in keys:
                    print(key)
                    value = float(field.split('=')[1])
                    justTopmedVariantsDict[mykey][key] = value

    for v in justTopmedVariantsDict:
        for e in ethnicitiesList:
            # get non-topmed ac and an
            ntmACKey = 'AC-non_topmed-' + e
            ntmANKey = 'AN-non_topmed-' + e
            ntmAC = justTopmedVariantsDict[v][ntmACKey]
            ntmAN = justTopmedVariantsDict[v][ntmANKey]

            # get topmed ac and an
            tmACKey = 'AC-' + e
            tmANKey = 'AN-' + e
            tmAC = justTopmedVariantsDict[v][tmACKey]
            tmAN = justTopmedVariantsDict[v][tmANKey]

            # subtract
            ac = tmAC - ntmAC
            an = tmAN - ntmAN
            af = 0.0
            if an != 0:
                af = float(ac) / float(an)
            afKey = 'AF-just_topmed-' + e
            justTopmedVariantsDict[v][afKey] = af


def getNontopmedAlleleFreqs(vcfDF, keys):
    variantsDict = dict()
    # pull AF, AC, and AN and figure out how to deal with this!

    logger.info('getting allele freqs')
    for i in range(len(vcfDF)):
        fields = vcfDF.iloc[i]['INFO'].split(';')
        chrom = vcfDF.iloc[i]['#CHROM'].split('chr')[1]
        pos = vcfDF.iloc[i]['POS']
        ref = vcfDF.iloc[i]['REF']
        alt = vcfDF.iloc[i]['ALT']
        mykey = "(" + chrom + ", " + pos + ", '" + ref + "', '" + alt + "')"
        variantsDict[mykey] = dict()
        for field in fields:
            if '=' in field:
                key = field.split('=')[0]
                if key in keys:
                    value = float(field.split('=')[1])
                    variantsDict[mykey][key] = value

    return variantsDict

def createDicts(variantsDict, topmedKeys, nontopmedKeys):
    # create dict for topmed and non-topmed
    topmedDict = dict()
    nontopmedDict = dict()

    topmedKeys.sort()
    nontopmedKeys.sort()

    for key in topmedKeys:
        topmedDict[key] = list()
    for key in nontopmedKeys:
        nontopmedDict[key] = list()
    for variant in variantsDict:
        for key in topmedKeys:
            topmedDict[key].append(variantsDict[variant][key])
        for key in nontopmedKeys:
            nontopmedDict[key].append(variantsDict[variant][key])
    return topmedDict, nontopmedDict

def plotDists(topmedDict, nontopmedDict, topmedKeys, nontopmedKeys, graphFileName):

    n=len(topmedDict[topmedKeys[0]])

    for i in range(len(topmedKeys)):
        # plot scatter
        tmkey = topmedKeys[i]
        ntmkey = nontopmedKeys[i]

        nontopmedList = nontopmedDict[ntmkey]
        topmedList = topmedDict[tmkey]
        logntmList = list()
        for i in range(len(nontopmedList)):
            if nontopmedList[i] == 0:
                logntmList.append(0.0)
            else:
                logntmList.append(math.log(nontopmedList[i], 10))

        logtmList = list()
        logjusttmList = list()
        for i in range(len(topmedList)):
            if topmedList[i] == 0:
                logtmList.append(0.0)
            else:
                logtmList.append(math.log(topmedList[i], 10))

        lowerBound = min([min(logntmList), min(logtmList)])
        upperBound = max([max(logntmList), max(logtmList)])
        lineNumbers = numpy.arange(lowerBound, upperBound, 0.1)

        # plot all
        plt.scatter(logntmList, logtmList, marker='.', color='black')
        plt.scatter(lineNumbers, lineNumbers, marker='.', color='red')
        plt.ylabel('log10(topmed AF)', fontsize=18)
        plt.xlabel('log10(nontopmed AF)', fontsize=18)
        plt.title(graphFileName + '_' + tmkey + '_vs_' + ntmkey + '_scatter_' + ' n=' + str(n))
        plt.savefig(graphFileName + '_' + tmkey + '_vs_' + ntmkey + '_scatter_' + '_n=' + str(n) + '.png')
        plt.close()

        # plot PDF
        lowerLimit = min(logntmList + logtmList)
        upperLimit = max(logntmList + logtmList)
        binSize = (upperLimit - lowerLimit) / 10
        plt.xlim(lowerLimit, upperLimit)
        bins = numpy.arange(lowerLimit, upperLimit, binSize)
        plt.hist([logntmList, logtmList], label=['log10(topmed AF)', 'log10(nontopmed AF)'], density=True, bins=bins)
        plt.xlabel('log10(AF)')
        plt.ylabel('count')
        plt.title(graphFileName + '_' + tmkey + '_vs_' + ntmkey + 'PDF' + ' n=' + str(n))
        plt.legend(loc="upper right")
        plt.savefig(graphFileName + '_' + tmkey + '_vs_' + ntmkey + '_PDF_' + '_n=' + str(n) + '.png')
        plt.close()


        # run KS test
        ksTest = ks_2samp(topmedList, nontopmedList)
        print('ksTest for ' + tmkey + ' vs ' + ntmkey + ' : ' + str(ksTest))



if __name__ == "__main__":
    main()