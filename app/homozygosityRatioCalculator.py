import sys
import json
import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd



logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)


def main():
    if len(sys.argv) != 3:
        print('13-vpi.json outputDir')
        sys.exit(1)

    vpiFileName = sys.argv[1]
    outputDir = sys.argv[2]

    with open(vpiFileName, 'r') as f:
        vpiDict = json.load(f)
    f.close()

    ratios = getRatios(vpiDict)
    with open(outputDir + '/hrpi.json', 'w') as f:
        json.dump(ratios, f)
    f.close()

    ratioList = list()
    for i in ratios:
        ratioList.append(ratios[i][2])

    binList = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    binPlot(ratioList, 10, "homozygosity density", "number of individuals", float, 3, binList, outputDir, '13-hrpi.png')

    scatterPlot(ratios, outputDir, '13-homovhet.png')

def getRatios(vpiDict):
    ratioDict = dict()
    for i in vpiDict:
        homoSum = 0
        hetSum = 0
        for b in vpiDict[i]['benign']:
            if b[1] == "3":
                homoSum += 1
            else:
                hetSum += 1
        for p in vpiDict[i]['pathogenic']:
            if p[1] == "3":
                homoSum += 1
            else:
                hetSum += 1
        for v in vpiDict[i]['vus']:
            if v[1] == "3":
                homoSum += 1
            else:
                hetSum += 1
        if homoSum + hetSum == 0:
            ratioDict[i] = [0, 0, 0]
        else:
            ratioDict[i] = [homoSum, hetSum, float(homoSum) / float(homoSum + hetSum)]
    return ratioDict

def scatterPlot(ratios, outputDir, imageName):
    plt.style.use('seaborn-whitegrid')
    homoList = list()
    hetList = list()
    for i in ratios:
        homoList.append(ratios[i][0])
        hetList.append(ratios[i][1])
    fig = plt.figure()  # set up plot
    plt.xlabel('homozygosity frequency')
    plt.ylabel('heterozygosity frequency')
    plt.scatter(homoList, hetList, s=10, color='black')
    plt.savefig(outputDir + '/' + imageName)


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
        df_bins.groupby(xlabel, as_index=False).size().plot(kind='bar')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        #plt.xlim(start, end)
        #plt.ylim(ymin, ymax)
        #plt.show()
        plt.savefig(outputDir + '/' + imageName)


if __name__ == "__main__":
    main()