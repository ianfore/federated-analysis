import matplotlib.pyplot as plt
import logging
import numpy as np
import pandas as pd
import sys

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
    for i in range(len(pList)):
        pList[i] = pList[i] **2

    qList = list(df_1['q']) + list(df_2['q'])
    for i in range(len(qList)):
        qList[i] = (qList[i]) **2

    twoPQList = list()
    for i in range(len(qList)):
        twoPQList.append(2 * pList[i] * qList[i])

    hwDF = pd.DataFrame()
    hwDF['q**2'] = qList
    hwDF['p**2'] = pList
    hwDF['2pq'] = twoPQList
    ax = hwDF.plot.kde()

    plt.xlim(0, 1)
    plt.ylim(0, 1)
    #plt.show()
    plt.savefig(outputDir + '/' + 'hw-obs-dists.png')


def binPlot(theList, binSize, xlabel, ylabel, dtype, sigDigs, binList, outputDir, imageName):
    customBinList = False
    if binList is None:
        sizeOfRange = 0.5 * (min(theList) + max(theList))
        try:
            binList = np.arange(min(theList), max(theList), round((1/binSize) * sizeOfRange, sigDigs) , dtype=dtype)
        except Exception as e:
            logger.error('error in binPlot np.arange(): ' + str(e))
            return
    else:
        customBinList = True
    bins = list()

    for element in theList:
        for i in range(len(binList) - 1):
            lhs = binList[i]
            rhs = binList[i + 1]
            if element >= lhs and element <= rhs:
                if customBinList:
                    bins.append(rhs)
                else:
                    bins.append(dtype(round(0.5 * (lhs + rhs), sigDigs)))
                break

    binMax = binList[len(binList)-1]
    listMax = max(theList)
    for element in theList:
        if element > binMax:
            bins.append(dtype(round(listMax, sigDigs)))



    df_bins = pd.DataFrame({xlabel: bins})
    if (len(df_bins) != 0):
        fontsize=12
        labelsize=7
        df_bins.groupby(xlabel, as_index=False).size().plot(kind='bar')
        plt.xlabel(xlabel, fontsize=fontsize)
        plt.ylabel(ylabel, fontsize=fontsize)
        plt.rc('xtick', labelsize=labelsize)
        plt.rc('ytick', labelsize=labelsize)
        plt.tight_layout()
        #plt.xlim(start, end)
        #plt.ylim(ymin, ymax)
        plt.show()
        #plt.savefig(outputDir + '/' + imageName)

if __name__ == "__main__":
    main()