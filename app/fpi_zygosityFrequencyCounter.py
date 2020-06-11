import sys
import json
import logging
import matplotlib.pyplot as plt
import numpy as np



logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)


def main():
    if len(sys.argv) != 4:
        print('13 13-fpi.json outputDir')
        sys.exit(1)

    chrom = sys.argv[1]
    fpiFileName = sys.argv[2]
    outputDir = sys.argv[3]

    with open(fpiFileName, 'r') as f:
        fpiDict = json.load(f)
    f.close()

    freqDict, ratioArray = getFrequenciesAndRatios(fpiDict)
    with open(outputDir + '/' + chrom + '-rpi.json', 'w') as f:
        json.dump(freqDict, f)
    f.close()

    binList = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    benignList = list()
    pathogenicList = list()
    vusList = list()
    for individual in freqDict:
        homo = freqDict[individual]['homoBSum']
        het = freqDict[individual]['hetBSum']
        if homo + het == 0:
            pass
        else:
            benignList.append(float(homo)/float(homo + het))

        homo = freqDict[individual]['homoPSum']
        het = freqDict[individual]['hetPSum']
        if homo + het == 0:
            pass
        else:
            pathogenicList.append(float(homo) / float(homo + het))

        homo = freqDict[individual]['homoVSum']
        het = freqDict[individual]['hetVSum']
        if homo + het == 0:
            pass
        else:
            vusList.append(float(homo) / float(homo + het))



    plt.hist([benignList, pathogenicList, vusList], binList, stacked = True, alpha=0.5, color=['green', 'red', 'blue'], label=['benign', 'pathogenic', 'vus'])
    plt.xlabel('homozygosity ratio')
    plt.ylabel('number of individuals')
    plt.title('ratios for chr ' + chrom)
    plt.xticks(np.arange(0.0, 1.1, 0.1))
    plt.legend(loc='upper right')
    plt.savefig(outputDir + '/' + chrom + '-rpi.png')
    plt.show()



def getFrequenciesAndRatios(fpiDict):
    frequencyDict = dict()
    ratioArray = np.zeros((len(fpiDict), 3))
    print('num individuals = ' + str(len(fpiDict)))
    counter = 0
    for i in fpiDict:
        homoBSum = fpiDict[i]['benign']['homo']
        homoPSum = fpiDict[i]['pathogenic']['homo']
        homoVSum = fpiDict[i]['vus']['homo']
        hetBSum = fpiDict[i]['benign']['hetero']
        hetPSum = fpiDict[i]['pathogenic']['hetero']
        hetVSum = fpiDict[i]['vus']['hetero']

        frequencyDict[i] = {'homoBSum': homoBSum, 'hetBSum': hetBSum,
                        'homoPSum': homoPSum, 'hetPSum': hetPSum,
                        'homoVSum': homoVSum, 'hetVSum': hetVSum }
        if homoBSum + hetBSum == 0:
            #ratioArray[counter][0] = 0
            pass
        else:
            ratioArray[counter][0] = homoBSum / (homoBSum + hetBSum)
        if homoPSum + hetPSum == 0:
            #ratioArray[counter][1] = 0
            pass
        else:
            ratioArray[counter][1] = homoPSum / (homoPSum + hetPSum)
        if homoVSum + hetVSum == 0:
            #ratioArray[counter][2] = 0
            pass
        else:
            ratioArray[counter][2] = homoVSum / (homoVSum + hetVSum)

        counter += 1

    return frequencyDict, ratioArray

def homoVHet(ratios, outputDir, imageName):
    plt.style.use('seaborn-whitegrid')
    homoList = list()
    hetList = list()
    for i in ratios:
        homoList.append(ratios[i][0])
        hetList.append(ratios[i][1])
    fig = plt.figure()
    plt.xlabel('homozygosity frequency')
    plt.ylabel('heterozygosity frequency')
    plt.scatter(homoList, hetList, s=10, color='black')
    plt.savefig(outputDir + '/' + imageName)



if __name__ == "__main__":
    main()