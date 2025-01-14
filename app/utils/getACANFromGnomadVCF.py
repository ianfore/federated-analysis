import pandas as pd
import logging
import argparse
import json
import matplotlib.pyplot as plt
import numpy
from scipy.stats import ks_2samp
import math
from statsmodels.graphics.gofplots import qqplot_2samples
import statsmodels.api as sm

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.INFO)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-gv', '--gnomadVCF', help='gnomad vcf input file')
    parser.add_argument('-tv', '--topmedVCF', help='topmed vcf input file')
    parser.add_argument('-j', '--json', help='json output file')
    parser.add_argument('-gp', '--graphPrefix', help='graph output file prefix')
    options = parser.parse_args()
    return options


def main():
    gnomadVCFFile = parse_args().gnomadVCF
    topmedVCFFile = parse_args().topmedVCF
    outputFileName = parse_args().json
    graphFileName = parse_args().graphPrefix


    logger.info('finding variants from ' + gnomadVCFFile)
    gnomadVCFDF = pd.read_csv(gnomadVCFFile, delimiter='\t', header=0, dtype=str)
    nontopmedKeys = ['AF_non_topmed_afr', 'AF_non_topmed_amr', 'AF_non_topmed_nfe',
                     'AF_non_topmed_fin', 'AF_non_topmed_eas', 'AF_non_topmed_sas']

    topmedKeys = ['AF_just_topmed_afr', 'AF_just_topmed_amr', 'AF_just_topmed_nfe',
                     'AF_just_topmed_fin', 'AF_just_topmed_eas', 'AF_just_topmed_sas']

    logger.info('getting nontopmed allele freqs')
    allVariants = getNontopmedAlleleFreqsFromGnomad(gnomadVCFDF, nontopmedKeys)

    logger.info('getting topmed allele freqs')
    getTopmedAlleleFreqsFromGnomad(gnomadVCFDF, allVariants)

    logger.info('reading data from topmed vcf')
    topmedVCFDF = pd.read_csv(topmedVCFFile, delimiter='\t', header=0, dtype=str)


    logger.info('saving ntm allele freqs')
    with open(outputFileName, 'w') as f:
        json.dump(allVariants, f)
    f.close()

    logger.info('plotting allele freqs')
    topmedDict, nontopmedDict = createDicts(allVariants, topmedKeys, nontopmedKeys)


    plotDists(topmedDict, nontopmedDict, topmedKeys, nontopmedKeys, graphFileName)

def getTopmedAlleleFreqsFromGnomad(vcfDF, allVariants):
    ethnicitiesList = ['afr', 'amr', 'eas', 'fin', 'nfe', 'sas']
    # AC-non_topmed-afr, AN-non_topmed-afr, AC-afr, AN-afr
    keys = ['AC', 'AN', 'AF']
    for e in ethnicitiesList:
        ntmACKey = 'AC_non_topmed_' + e
        keys.append(ntmACKey)
        ntmANKey = 'AN_non_topmed_' + e
        keys.append(ntmANKey)
        tmACKey = 'AC_' + e
        keys.append(tmACKey)
        tmANKey = 'AN_' + e
        keys.append(tmANKey)


    # get AC and AN for each ethnicity, then calculate AF
    for i in range(len(vcfDF)):
        fields = vcfDF.iloc[i]['INFO'].split(';')
        chrom = vcfDF.iloc[i]['#CHROM'].split('chr')[1]
        pos = vcfDF.iloc[i]['POS']
        ref = vcfDF.iloc[i]['REF']
        alt = vcfDF.iloc[i]['ALT']
        v = "(" + chrom + ", " + pos + ", '" + ref + "', '" + alt + "')"
        for field in fields:
            if '=' in field:
                key = field.split('=')[0]
                if key in keys:
                    value = float(field.split('=')[1])
                    allVariants[v][key] = value

    for v in allVariants:
        for e in ethnicitiesList:
            afKey = 'AF_just_topmed_' + e
            # get non-topmed ac and an
            ntmACKey = 'AC_non_topmed_' + e
            ntmANKey = 'AN_non_topmed_' + e
            # odd hack that sometimes AN and/or AC not in VCF for certain variants
            if not ntmANKey in allVariants[v] or not ntmACKey in allVariants[v]:
                allVariants[v][afKey] = 0.0
            else:
                ntmAN = allVariants[v][ntmANKey]
                ntmAC = allVariants[v][ntmACKey]
                # get topmed ac and an
                tmACKey = 'AC_' + e
                tmANKey = 'AN_' + e
                tmAC = allVariants[v][tmACKey]
                tmAN = allVariants[v][tmANKey]

                # subtract
                ac = tmAC - ntmAC
                an = tmAN - ntmAN
                af = 0.0
                if an != 0:
                    af = float(ac) / float(an)
                allVariants[v][afKey] = af


def getNontopmedAlleleFreqsFromGnomad(vcfDF, keys):

    variantsDict = dict()
    for i in range(len(vcfDF)):
        fields = vcfDF.iloc[i]['INFO'].split(';')
        chrom = vcfDF.iloc[i]['#CHROM'].split('chr')[1]
        pos = vcfDF.iloc[i]['POS']
        ref = vcfDF.iloc[i]['REF']
        alt = vcfDF.iloc[i]['ALT']
        variant = "(" + chrom + ", " + pos + ", '" + ref + "', '" + alt + "')"
        variantsDict[variant] = dict()
        for field in fields:
            if '=' in field:
                key = field.split('=')[0]
                if key in keys:
                    value = float(field.split('=')[1])
                    variantsDict[variant][key] = value

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
            # key = AF-just_topmed-eth
            ntmKey = key.replace('just', 'non')
            # odd hack that sometimes a variant doesn't have AC and AN for all ethnicities
            if key in variantsDict[variant] and ntmKey in variantsDict[variant]:
                topmedDict[key].append(variantsDict[variant][key])
        for key in nontopmedKeys:
            # key = AF-non_topmed-eth
            tmKey = key.replace('non', 'just')
            if key in variantsDict[variant] and tmKey in variantsDict[variant]:
                nontopmedDict[key].append(variantsDict[variant][key])
    return topmedDict, nontopmedDict

def plotScatter(logntmList, logtmList, lineNumbers, graphFileName, tmkey, ntmkey):
    n=len(logntmList)
    plt.scatter(logntmList, logtmList, marker='.', color='black')
    plt.scatter(lineNumbers, lineNumbers, marker='.', color='red')
    plt.ylabel('log10(topmed AF)', fontsize=18)
    plt.xlabel('log10(non-topmed AF)', fontsize=18)
    plt.title(graphFileName + '_' + tmkey + '_vs_' + ntmkey + '_scatter_' + ' n=' + str(n))
    plt.savefig(graphFileName + '_' + tmkey + '_vs_' + ntmkey + '_scatter_' + '_n=' + str(n) + '.png')
    plt.close()

def plotHist(logntmList, logtmList, graphFileName, tmkey, ntmkey):
    n=len(logntmList)
    lowerLimit = min(logntmList + logtmList)
    upperLimit = max(logntmList + logtmList)
    binSize = (upperLimit - lowerLimit) / 20
    plt.xlim(lowerLimit, upperLimit)
    bins = numpy.arange(lowerLimit, upperLimit, binSize)
    plt.hist([logntmList, logtmList], label=['log10(non_topmed AF)', 'log10(just_topmed AF)'], bins=bins)
    plt.xlabel('log10(AF)')
    plt.ylabel('count')
    plt.title(graphFileName + '_' + tmkey + '_vs_' + ntmkey + '_PDF_' + ' n=' + str(n))
    plt.legend(loc="upper right")
    plt.savefig(graphFileName + '_' + tmkey + '_vs_' + ntmkey + '_PDF_' + 'n=' + str(n) + '.png')
    plt.close()

def createListsPerEthnicity(nontopmedDict, topmedDict, ntmkey, tmkey):
    nontopmedList = nontopmedDict[ntmkey]
    topmedList = topmedDict[tmkey]
    logntmList = list()
    logtmList = list()

    for i in range(len(topmedList)):
        if topmedList[i] == 0 or nontopmedList[i] == 0:
            continue
        else:
            logtmList.append(math.log(topmedList[i], 10))
            logntmList.append(math.log(nontopmedList[i], 10))

    return logntmList, logtmList

def createListsFromDict(nontopmedDict, topmedDict):
    logntmList = list()
    logtmList = list()
    '''for key in nontopmedDict:
        for j in range(len(nontopmedDict[key])):
            if nontopmedDict[key][j] == 0:
                logntmList.append(0.0)
            else:
                logntmList.append(math.log(nontopmedDict[key][j], 10))
    for key in topmedDict:
        for k in range(len(topmedDict[key])):
            if topmedDict[key][k] == 0:
                logtmList.append(0.0)
            else:
                logtmList.append(math.log(topmedDict[key][k], 10))'''
    ntmKeys = list(nontopmedDict.keys())
    tmKeys = list(topmedDict.keys())
    for i in range(len(tmKeys)):
        tmKey = tmKeys[i]
        ntmKey = ntmKeys[i]
        for j in range(len(nontopmedDict[ntmKey])):
            if nontopmedDict[ntmKey][j] != 0 and topmedDict[tmKey][j] != 0:
                logtmList.append(math.log(topmedDict[tmKey][j], 10))
                logntmList.append(math.log(nontopmedDict[ntmKey][j], 10))


    return logntmList, logtmList


def plotDists(topmedDict, nontopmedDict, topmedKeys, nontopmedKeys, graphFileName):

    logntmList, logtmList = createListsFromDict(nontopmedDict, topmedDict)
    lowerBound = min([min(logntmList), min(logtmList)])
    upperBound = max([max(logntmList), max(logtmList)])
    lineNumbers = numpy.arange(lowerBound, upperBound, 0.1)

    # plot scatter
    plotScatter(logntmList, logtmList, lineNumbers, graphFileName, 'non_topmed', 'just_topmed')

    # plot PDF
    plotHist(logntmList, logtmList, graphFileName, 'non_topmed', 'just_topmed')

    # plot QQ-plot
    n = len(logntmList)
    plt.title('all non-topmed vs all just-topmed QQ n=' + str(n))
    ax = plt.gca()
    qqplot_2samples(data1=sm.ProbPlot(numpy.array(logtmList)), data2=sm.ProbPlot(numpy.array(logntmList)),
                    xlabel='non-topmed', ylabel='just-topmed',
                    line="45",ax = ax)
    plt.savefig(graphFileName + '_' + 'non_topmed' + '_vs_' + 'just_topmed' + '_QQ_n=' + str(n) + '.png')
    plt.close()

    # create non-zero lists
    nonZeroTM = [x for x in logtmList if x != 0]
    nonZeroNTM = [x for x in logntmList if x != 0]

    # run KS test
    # ksTest = ks_2samp(topmedDict[tmkey], nontopmedDict[ntmkey])
    ksTest = ks_2samp(nonZeroNTM, nonZeroTM)

    print('ksTest for non-zero: ' + 'just_topmed' + ' vs ' + 'non_topmed' + ' : ' + str(ksTest))

    for i in range(len(topmedKeys)):
        tmkey = topmedKeys[i]
        ntmkey = nontopmedKeys[i]

        logntmList, logtmList = createListsPerEthnicity(nontopmedDict, topmedDict, ntmkey, tmkey)

        lowerBound = min([min(logntmList), min(logtmList)])
        upperBound = max([max(logntmList), max(logtmList)])
        lineNumbers = numpy.arange(lowerBound, upperBound, 0.1)

        # plot scatter
        plotScatter(logntmList, logtmList, lineNumbers, graphFileName, tmkey, ntmkey)

        # plot PDF
        plotHist(logntmList, logtmList, graphFileName, tmkey, ntmkey)

        # plot QQ-plot
        n = len(logntmList)
        plt.title(graphFileName + '_' + ntmkey + '_vs_' + tmkey + '_QQ_' + 'n=' + str(n) )
        ax = plt.gca()
        ntmData = sm.ProbPlot(numpy.array(logntmList))
        tmData = sm.ProbPlot(numpy.array(logtmList))
        qqplot_2samples(data1=tmData, data2=ntmData,
                        xlabel='non-topmed', ylabel='just-topmed',
                        line="45", ax=ax)
        plt.savefig(graphFileName + '_' + ntmkey + '_vs_' + tmkey + '_QQ_n=' + str(n) + '.png')

        plt.close()

        # create non-zero lists
        nonZeroTM = [x for x in topmedDict[tmkey] if x!= 0 ]
        nonZeroNTM = [x for x in nontopmedDict[ntmkey] if x !=0 ]

        # run KS test
        #ksTest = ks_2samp(topmedDict[tmkey], nontopmedDict[ntmkey])
        ksTest = ks_2samp(nonZeroNTM, nonZeroTM)

        print('ksTest for non-zero: ' + ntmkey + ' vs ' + tmkey + ' : ' + str(ksTest))

        # TODO stratify by < 0.01
        # TODO combine / unstratify by ethnicities
        # TODO fix labels
        # TODO broad v non-broad
        # TODO compare gnomad VC vs topmed VC (singletons, all, ?)

if __name__ == "__main__":
    main()