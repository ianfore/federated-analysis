import matplotlib.pyplot as plt
import logging
import numpy as np
import pandas as pd
import sys
from itertools import groupby


logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

def main():
    if len(sys.argv) != 4:
        print('brca1-report.tsv brca2-report.tsv output-dir')
        sys.exit(1)

    brca1_report = sys.argv[1]
    brca2_report = sys.argv[2]
    outputDir = sys.argv[3]

    df_1 = pd.read_csv(brca1_report, header=0, sep='\t')
    df_2 = pd.read_csv(brca2_report, header=0, sep='\t')

    pList = list(df_1['p']) + list(df_2['p'])
    thePList = list()
    for i in range(len(pList)):
        thePList.append(pList[i] **2)

    qList = list(df_1['q']) + list(df_2['q'])
    theQList = list()
    for i in range(len(qList)):
        theQList.append((qList[i] - 1) **2)

    twoPQList = list()
    for i in range(len(qList)):
        twoPQList.append(2 * pList[i]  * qList[i])

    '''x, y = np.unique(pList,return_counts=True)
    plt.scatter(x, y)
    x, y = np.unique(qList, return_counts=True)
    plt.scatter(x, y)
    x, y = np.unique(twoPQList, return_counts=True)
    plt.scatter(x, y)'''
    '''pList.sort()
    qList.sort()
    twoPQList.sort()'''

    '''hwDF = pd.DataFrame()
    hwDF['q**2'] = qList
    hwDF['p**2'] = pList
    hwDF['2pq'] = twoPQList
    hwDF['q**2'].value_counts(normalize=True).plot()
    ax = hwDF.plot.kde()
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.show()'''
    #plt.savefig(outputDir + '/' + 'hw-obs-dists.png')

    pDF = countBins(thePList, 100)
    pDF['AF'] = pDF.index
    print(pDF)
    #plt.scatter(x=pDF['AF'], y=pDF[0])

    qDF = countBins(theQList, 100)
    qDF['AF'] = qDF.index
    qDF["AF"] = qDF["AF"].values[::-1]
    print(qDF)
    #plt.scatter(x=qDF['AF'], y=qDF[0])

    twoPQDF = countBins(twoPQList, 100)
    twoPQDF['AF'] = twoPQDF.index
    print(twoPQDF)
    plt.scatter(x=twoPQDF['AF'], y=twoPQDF[0])

    plt.show()

def countBins(theList, numBins):
    n = len(theList)
    start = theList[0]
    binCounts = dict()
    binSize = (max(theList) - min(theList)) / numBins
    for i in range(numBins):
        end = start + binSize
        binCounts['%.1f'%(start)] = (float(len([j for j in theList if (j >= start and j < end)])) / float(n)) * numBins/100
        start = end
    return pd.DataFrame.from_dict(binCounts, orient='index')

if __name__ == "__main__":
    main()