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
    variantsDict = dict()
    topmedKeys = ['AF-afr', 'AF-amr', 'AF-nfe', 'AF-eas', 'AF-sas', 'AF-oth']
    nontopmedKeys = ['AF-non_topmed-afr', 'AF-non_topmed-amr', 'AF-non_topmed-nfe',
                     'AF-non_topmed-eas', 'AF-non_topmed-sas', 'AF-non_topmed-oth']
    keys = topmedKeys + nontopmedKeys

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

    logger.info('saving allele freqs')
    with open(outputFileName, 'w') as f:
        json.dump(variantsDict, f)
    f.close()

    logger.info('plotting allele freqs')
    plotDists(variantsDict, topmedKeys, nontopmedKeys, graphFileName)

def plotDists(variantsDict, topmedKeys, nontopmedKeys, graphFileName):

    # plot 'Q-Q'
    lineNumbers = numpy.arange(0, 1, 0.01)
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

    n=len(topmedDict[topmedKeys[0]])


    for i in range(len(topmedKeys)):
        # plot scatter
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        tmkey = topmedKeys[i]
        ntmkey = nontopmedKeys[i]

        nontopmedList = nontopmedDict[ntmkey]
        topmedList = topmedDict[tmkey]
        plt.scatter(nontopmedList, topmedList, marker='.', color='black')
        plt.scatter(lineNumbers, lineNumbers, marker='.', color='red')
        plt.ylabel('topmed AF', fontsize=18)
        plt.xlabel('nontopmed AF', fontsize=18)
        plt.title(graphFileName + '_' + tmkey + '_vs_' + ntmkey + '_QQ_' + ' n=' + str(n))
        plt.savefig(graphFileName + '_' + tmkey + '_vs_' + ntmkey + '_QQ_' + '_n=' + str(n) + '.png')
        plt.close()

        # plot PDF
        logntmList = list()
        for i in range(len(nontopmedList)):
            if nontopmedList[i] == 0:
                logntmList.append(0.0)
            else:
                logntmList.append(math.log(nontopmedList[i]))

        logtmList = list()
        for i in range(len(topmedList)):
            if topmedList[i] == 0:
                logtmList.append(0.0)
            else:
                logtmList.append(math.log(topmedList[i]))

        lowerLimit = min(logntmList + logtmList)
        upperLimit = max(logntmList + logtmList)
        plt.xlim(lowerLimit, upperLimit)
        bins = numpy.arange(lowerLimit, upperLimit)
        #plt.hist([logntmList, logtmList], label=['log-topmed', 'log-nontopmed'], range=(0,1), density=True, bins=bins)
        plt.hist([logntmList, logtmList], label=['log-topmed', 'log-nontopmed'], density=True, bins=bins)
        plt.xlabel('AF')
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