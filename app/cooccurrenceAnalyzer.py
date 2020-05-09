import sys
import json
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import pandas as pd
import math

coordinateColumnBase = 'Genomic_Coordinate_hg'
brcaFileName = '/Users/jcasaletto/PycharmProjects/BIOBANK/federated-analysis/data/brca-variants.tsv'
hgVersion = 38


def main():
    if len(sys.argv) != 4:
        print('provide path to variant, vpi, and brca files as args')
        print('13-out.json 13-vpi.json brca-variants.tsv')
        sys.exit(1)

    variantsFileName =sys.argv[1]
    vpiFileName = sys.argv[2]
    brcaFileName = sys.argv[3]

    with open(vpiFileName, 'r') as f:
        vpiDict = json.load(f)
    f.close()

    brcaDF = findVariantsInBRCA(brcaFileName)

    with open(variantsFileName, 'r') as f:
        variantsDict = json.load(f)
    f.close()

    genotypeCounts, frequenciesPerIndividual = countTotalGenotypesForVariants(vpiDict, 0.001, brcaDF, hgVersion, rare=True)
    plotGenotypeCounts(genotypeCounts, rare=True)
    #plotFrequenciesPerIndividual(frequenciesPerIndividual)




    '''variantCounts = countTotalVariants(vpiDict)
    print('benign counts: ' + str(len(variantCounts['benign'])))
    print('pathogenic counts: ' + str(len(variantCounts['pathogenic'])))
    print('vus counts: ' + str(len(variantCounts['vus'])))'''


    #plotVUSByPosition(variantsDict)

    #brcaDF = findVariantsInBRCA(brcaFileName)
    #plotVUSByFrequency(variantsDict, 'maxPopFreq', brcaDF, hgVersion)
    #plotVUSByFrequency(variantsDict, 'cohortFreq', brcaDF, hgVersion)

    #homoVhetero = countHomoAndHeteroPerIndividual(vpiDict, variantsDict, brcaDF, hgVersion)
    #print(homoVhetero)

    #printHWReport(vpiDict, variantsDict)


def countTotalVariants(vpiDict):
    variants = {'benign': set(), 'pathogenic': set(), 'vus': set()}
    for i in vpiDict:
        for v in vpiDict[i]['benign']:
            variants['benign'].add(tuple(v[0]))
        for v in vpiDict[i]['pathogenic']:
            variants['pathogenic'].add(tuple(v[0]))
        for v in vpiDict[i]['vus']:
            variants['vus'].add(tuple(v[0]))

    return variants

def countHomoAndHeteroPerIndividual(vpiDict, variantsDict, brcaDF, hgVersion):
    individualsPerVariant = defaultdict(list)
    # look at each homo vus
    for homoVUS in variantsDict['homozygous vus']:
        foundOne = False
        maxPopFreq = variantsDict['homozygous vus'][homoVUS]['maxPopFreq']
        if  maxPopFreq > 0.001:
            continue
        for individual in vpiDict:
            # find the individuals who have expressed this homo vus
            for v in vpiDict[individual]['vus']:
                if tuple(v[0]) == eval(homoVUS):

                    for b in vpiDict[individual]['benign']:
                        freq = getGnomadData(brcaDF, b[0], hgVersion)[1]
                        if b[1] == "1" or b[1] == "2" and  freq > 0.01:
                            individualsPerVariant[homoVUS].append(individual)
                            foundOne = True
                            break
                break
            if foundOne:
                break

    # now see if the individuals who have the homo VUS also have common benign hetero SNPs
    # this is evidence of genotype errors


    return individualsPerVariant

def findVariantsInBRCA(fileName):
    return pd.read_csv(fileName, sep='\t', header=0, dtype=str)




def getGnomadData(brcaDF, vus, hgVersion):
    # TODO write a unit test
    # 13:g.32393468:C>CT
    #hgString = 'chr' + str(vus[0][0]) + ':g.' + str(vus[0][1]) + ':' + str(vus[0][2]) + '>' + str(vus[0][3])
    hgString = 'chr' + str(vus[0]) + ':g.' + str(vus[1]) + ':' + str(vus[2]) + '>' + str(vus[3])

    # first, get list of columns for GnomAD allleles
    gnomad = [v for v in list(brcaDF.columns) if 'GnomAD' in v]
    alleleFrequencies = [v for v in gnomad if 'Allele_frequency' in v]

    # second, get frequencies across exomes and genomes to determine max
    # return population, frequency, count, and number
    # replace "frequency" with "count" and "number" in Allele_frequency_genome_AFR_GnomAD

    maxFrequency = 0.0
    maxPopulation = None
    for af in alleleFrequencies:
        freq=0.0
        alleleFreqList = brcaDF[brcaDF[coordinateColumnBase + str(hgVersion)] == hgString][af].tolist()
        if alleleFreqList:
            try:
                freq = float(alleleFreqList[0])
            except ValueError:
                continue
            if freq > maxFrequency:
                maxFrequency = freq
                maxPopulation = af

    return (maxPopulation, maxFrequency)

def printHWReport(vpiDict, variantsDict):
    bVars, pVars, vVars = calculateZygosityFrequenciesPerVariant(vpiDict)
    bVars, pVars, vVars = hardyWeinbergChiSquareTest(bVars, pVars, vVars, len(vpiDict))
    rejectHW = {'benign': 0, 'pathogenic': 0, 'vus': 0}
    acceptHW = {'benign': 0, 'pathogenic': 0, 'vus': 0}
    for b in bVars:
        if bVars[b]['accept hw'] is False:
            rejectHW['benign'] += 1
        else:
            acceptHW['benign'] += 1
    for p in pVars:
        if pVars[p]['accept hw'] is False:
            rejectHW['pathogenic'] += 1
        else:
            acceptHW['pathogenic'] += 1

    rejectVUS = {'cooccurring vus': 0, 'homozygous vus': 0}
    acceptVUS = {'cooccurring vus': 0, 'homozygous vus': 0}

    for v in vVars:
        if vVars[v]['accept hw'] is False:
            rejectHW['vus'] += 1
            if str(v) in variantsDict['cooccurring vus']:
                rejectVUS['cooccurring vus'] += 1
            if str(v) in variantsDict['homozygous vus']:
                rejectVUS['homozygous vus'] += 1
        else:
            acceptHW['vus'] += 1
            if str(v) in variantsDict['cooccurring vus']:
                acceptVUS['cooccurring vus'] += 1
            if str(v) in variantsDict['homozygous vus']:
                acceptVUS['homozygous vus'] += 1

    print('reject HW: ' + str(rejectHW))
    print('accept HW: ' + str(acceptHW))
    print('num co-occurring vus that reject HW: ' + str(rejectVUS['cooccurring vus']))
    print('num co-occurring vus that accept HW:' + str(acceptVUS['cooccurring vus']))
    print('num homozygous vus that reject HW: ' + str(rejectVUS['homozygous vus']))
    print('num homozygous vus that accept HW: ' + str(acceptVUS['homozygous vus']))


def calculateZygosityFrequenciesPerVariant(vpiDict):
    benignVariants = dict()
    pathogenicVariants = dict()
    vusVariants = dict()
    for individual in vpiDict:
        for b in vpiDict[individual]['benign']:
            if b:
                if tuple(b[0]) not in benignVariants:
                    benignVariants[tuple(b[0])] = dict()
                    benignVariants[tuple(b[0])]['aa'] = 0
                    benignVariants[tuple(b[0])]['Aa'] = 0
                if b[1] == '3':
                    benignVariants[tuple(b[0])]['aa'] += 1
                else:
                    benignVariants[tuple(b[0])]['Aa'] += 1
        for p in vpiDict[individual]['pathogenic']:
            if p:
                if tuple(p[0]) not in pathogenicVariants:
                    pathogenicVariants[tuple(p[0])] = dict()
                    pathogenicVariants[tuple(p[0])]['aa'] = 0
                    pathogenicVariants[tuple(p[0])]['Aa'] = 0
                if p[1] == '3':
                    pathogenicVariants[tuple(p[0])]['aa'] += 1
                else:
                    pathogenicVariants[tuple(p[0])]['Aa'] += 1
        for v in vpiDict[individual]['vus']:
            if v:
                if tuple(v[0]) not in vusVariants:
                    vusVariants[tuple(v[0])] = dict()
                    vusVariants[tuple(v[0])]['aa'] = 0
                    vusVariants[tuple(v[0])]['Aa'] = 0
                if v[1] == '3':
                    vusVariants[tuple(v[0])]['aa'] += 1
                else:
                    vusVariants[tuple(v[0])]['Aa'] += 1

    n = len(vpiDict)
    for b in benignVariants:
        benignVariants[b]['AA'] = n - (benignVariants[b]['Aa'] + benignVariants[b]['aa'])
    for p in pathogenicVariants:
        pathogenicVariants[p]['AA'] = n - (pathogenicVariants[p]['Aa'] + pathogenicVariants[p]['aa'])
    for v in vusVariants:
        vusVariants[v]['AA'] = n - (vusVariants[v]['Aa'] + vusVariants[v]['aa'])

    return benignVariants, pathogenicVariants, vusVariants

def hardyWeinbergChiSquareTest(bVars, pVars, vVars, n):
    # https://en.wikipedia.org/wiki/Hardy-Weinberg_principle
    # degrees of freedom = 1
    criticalValue = 3.84

    for b in bVars:
        # 2. calculate p = (2 x Obs(AA) + Obs(Aa)) / (2 x (Obs(AA) + Obs(Aa) + Obs(aa))
        bVars[b]['p'] = (2 * bVars[b]['AA'] + bVars[b]['Aa']) / ( 2 * (bVars[b]['AA'] + bVars[b]['Aa'] + bVars[b]['aa']))

        # 3. calculate q = 1 - p
        bVars[b]['q'] = 1 - bVars[b]['p']

        # 4. calculate Exp(AA) = p**2 x n
        expAA = n * bVars[b]['p'] **2

        # 5. calculate Exp(Aa) = 2 x p * q * n
        expAa = 2 * bVars[b]['p'] * bVars[b]['q'] * n

        # 6. calculate Exp(aa) = q**2 x n
        expaa = n * bVars[b]['q'] **2

        # 7. calculate chi-square = sum[ (O - E)**2 / E ]
        if expAA == 0 or expAa == 0 or expaa == 0:
            bVars[b]['chisquare'] = 0
        else:
            bVars[b]['chisquare'] = 1.0/expAA * (bVars[b]['AA'] - expAA)**2 + \
                                    1.0/expAa * (bVars[b]['Aa'] - expAa)**2 + \
                                    1.0/expaa * (bVars[b]['aa'] - expaa)**2
        # 8. compare against p-value for 1 degree of freedom at 0.05 significance (3.84)
        if bVars[b]['chisquare'] >= criticalValue:
            bVars[b]['accept hw'] = False
        else:
            bVars[b]['accept hw'] = True

    for p in pVars:
        # 2. calculate p = (2 x Obs(AA) + Obs(Aa)) / (2 x (Obs(AA) + Obs(Aa) + Obs(aa))
        pVars[p]['p'] = (2 * pVars[p]['AA'] + pVars[p]['Aa']) / ( 2 * (pVars[p]['AA'] + pVars[p]['Aa'] + pVars[p]['aa']))

        # 3. calculate q = 1 - p
        pVars[p]['q'] = 1 - pVars[p]['p']

        # 4. calculate Exp(AA) = p**2 x n
        expAA = n * pVars[p]['p'] **2

        # 5. calculate Exp(Aa) = 2 x p * q * n
        expAa = 2 * pVars[p]['p'] * pVars[p]['q'] * n

        # 6. calculate Exp(aa) = q**2 x n
        expaa = n * pVars[p]['q'] **2

        # 7. calculate chi-square = sum[ (O - E)**2 / E ]
        if expAA == 0 or expAa == 0 or expaa == 0:
            pVars[p]['chisquare'] = 0
        else:
            pVars[p]['chisquare'] = 1.0/expAA * (pVars[p]['AA'] - expAA)**2 + \
                                    1.0/expAa * (pVars[p]['Aa'] - expAa)**2 + \
                                    1.0/expaa * (pVars[p]['aa'] - expaa)**2
        # 8. compare against p-value for 1 degree of freedom at 0.05 significance (3.84)
        if pVars[p]['chisquare'] >= criticalValue:
            pVars[p]['accept hw'] = False
        else:
            pVars[p]['accept hw'] = True

    for v in vVars:
        # 2. calculate p = (2 x Obs(AA) + Obs(Aa)) / (2 x (Obs(AA) + Obs(Aa) + Obs(aa))
        vVars[v]['p'] = (2 * vVars[v]['AA'] + vVars[v]['Aa']) / ( 2 * (vVars[v]['AA'] + vVars[v]['Aa'] + vVars[v]['aa']))

        # 3. calculate q = 1 - p
        vVars[v]['q'] = 1 - vVars[v]['p']

        # 4. calculate Exp(AA) = p**2 x n
        expAA = n * vVars[v]['p'] **2

        # 5. calculate Exp(Aa) = 2 x p * q * n
        expAa = 2 * vVars[v]['p'] * vVars[v]['q'] * n

        # 6. calculate Exp(aa) = q**2 x n
        expaa = n * vVars[v]['q'] **2

        # 7. calculate chi-square = sum[ (O - E)**2 / E ]
        if expAA == 0 or expAa == 0 or expaa == 0:
            vVars[v]['chisquare'] = 0
        else:
            vVars[v]['chisquare'] = 1.0/expAA * (vVars[v]['AA'] - expAA)**2 + \
                                    1.0/expAa * (vVars[v]['Aa'] - expAa)**2 + \
                                    1.0/expaa * (vVars[v]['aa'] - expaa)**2

        # 8. compare against p-value for 1 degree of freedom at 0.05 significance (3.84)
        if vVars[v]['chisquare'] >= criticalValue:
            vVars[v]['accept hw'] = False
        else:
            vVars[v]['accept hw'] = True

    return bVars, pVars, vVars


def binPlot(theList, binSize, xlabel, ylabel, dtype, sigDigs, binList):
    sizeOfRange = 0.5 * (min(theList) + max(theList))
    customBinList = False
    if binList is None:
        binList = np.arange(min(theList), max(theList), round((1/binSize) * sizeOfRange, sigDigs) , dtype=dtype)
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
        plt.show()


def plotVUSByFrequency(variantsDict, freq, brcaDF, hgVersion):
    homoVUS = variantsDict['homozygous vus']
    popFreqs = list()
    for vus in homoVUS:
        popFreqs.append(homoVUS[vus][freq])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    bp = ax.boxplot(popFreqs)
    plt.xlabel(freq)
    plt.show()
    binList = [0.00001, 0.0001, 0.001, 0.01, 0.1]
    binPlot(popFreqs, 25, freq, "number of homo VUS", float, 3, binList)

    popFreqs = list()
    heteroVUS = variantsDict['cooccurring vus']
    for vus in heteroVUS:
        if vus in homoVUS:
            continue
        elif freq is 'maxPopFreq':
            #pFreq = getGnomadData(brcaDF, eval(vus), hgVersion)[1]
            pFreq = heteroVUS[vus]['allele frequencies']['maxPopFreq']
            popFreqs.append(pFreq)
        else:
            pFreq = heteroVUS[vus]['allele frequencies']['cohortFreq']
            popFreqs.append(pFreq)
    print(popFreqs)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    bp = ax.boxplot(popFreqs)
    plt.xlabel(freq)
    plt.show()
    binList = [0.00001, 0.0001, 0.001, 0.01, 0.1]
    binPlot(popFreqs, 25, freq, "number of hetero VUS", float, 3, binList)


def plotVUSByPosition(variantsDict):
    locations = list()
    homozygousVUS = variantsDict['homozygous vus']
    for vus in homozygousVUS:
        locations.append(int(eval(vus)[1]))
    binPlot(locations, 10000, "chromosome position bins", "number of homo VUS", int, 0, None)

    locations = list()
    heterozygousVUS = variantsDict['cooccurring vus']
    for vus in heterozygousVUS:
        if vus in homozygousVUS:
            continue
        else:
            locations.append(int(eval(vus)[1]))
    binPlot(locations, 10000, "chromosome position bins", "number of het VUS", int, 0, None)



def plotFrequenciesPerIndividual(frequenciesPerIndividual):
    # count the number of individuals per frequency
    homo_ben_counts = list()
    homo_path_counts = list()
    homo_vus_counts = list()
    hetero_ben_counts = list()
    hetero_path_counts = list()
    hetero_vus_counts = list()


    for individual in frequenciesPerIndividual:
        homo_ben_counts.append(frequenciesPerIndividual[individual]['benign']['homo'])
        hetero_ben_counts.append(frequenciesPerIndividual[individual]['benign']['hetero'])
        homo_path_counts.append(frequenciesPerIndividual[individual]['pathogenic']['homo'])
        hetero_path_counts.append(frequenciesPerIndividual[individual]['pathogenic']['hetero'])
        homo_vus_counts.append(frequenciesPerIndividual[individual]['vus']['homo'])
        hetero_vus_counts.append(frequenciesPerIndividual[individual]['vus']['hetero'])

    binPlot(homo_ben_counts, 10, "homozygous benign variant count bins", "number of individuals", int, 0, None)
    binPlot(hetero_ben_counts, 10, "heterozygous benign variant count bins", "number of individuals", int, 0, None)
    #binPlot(homo_path_counts, 10, "homozygous pathogenic variant count bins", "number of individuals", int, 0, None)
    binPlot(hetero_path_counts, 10, "heterozygous pathogenic variant count bins", "number of individuals", int, 0, [0, 1])
    binPlot(homo_vus_counts, 10, "homozygous vus variant count bins", "number of individuals", int, 0, None)
    binPlot(hetero_vus_counts, 10, "heterozygous vus variant count bins", "number of individuals", int, 0, None)



def plotGenotypeCounts(genotypeCounts, rare):
    homoCounts = [0 if genotypeCounts['benign']['homo'] == 0 else genotypeCounts['benign']['homo'],
                 0 if genotypeCounts['pathogenic']['homo'] == 0 else genotypeCounts['pathogenic']['homo'],
                 0 if genotypeCounts['vus']['homo'] == 0 else genotypeCounts['vus']['homo']]
    heteroCounts = [0 if genotypeCounts['benign']['hetero'] == 0 else genotypeCounts['benign']['hetero'],
                 0 if genotypeCounts['pathogenic']['hetero'] == 0 else genotypeCounts['pathogenic']['hetero'],
                 0 if genotypeCounts['vus']['hetero'] == 0 else genotypeCounts['vus']['hetero']]


    # plot bar graph
    '''labels = ['benign', 'pathogenic', 'vus']
    x = np.arange(len(labels))
    width = 0.35
    fig, ax = plt.subplots()
    ax.bar(x - width/2, homoCounts, width, label = 'homozygous')
    ax.bar(x + width/2, heteroCounts, width, label = 'heterozygous')
    ax.set_ylabel('log10(counts)')
    ax.set_title('count per zygosity per classification')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend(loc='upper center')
    fig.tight_layout()
    #plt.ylim(0, 7)
    plt.show()'''

    # plot pie chart
    # Pie chart, where the slices will be ordered and plotted counter-clockwise:
    plt.rcParams['font.size'] = 18
    fig1, ax1 = plt.subplots()

    '''if rare:
        labels = ['Homo VUS',  'Hetero VUS']
        colors = ['green', 'brown']
        sizes = [genotypeCounts['vus']['homo'], genotypeCounts['vus']['hetero']]
        explode = (0, 0)
        ax1.pie(sizes, explode=explode, colors=colors, labels=labels, shadow=False, startangle=90)'''
    if False:
        print('hi')

    else:
        labels = ['Homo benign', 'Homo path', 'Homo VUS', 'Hetero benign', 'Hetero path', 'Hetero VUS']
        sizes = [genotypeCounts['benign']['homo'], genotypeCounts['pathogenic']['homo'],
                genotypeCounts['vus']['homo'], genotypeCounts['benign']['hetero'],
                genotypeCounts['pathogenic']['hetero'], genotypeCounts['vus']['hetero']]
        colors = ['red', 'yellow', 'green', 'orange', 'blue', 'brown']
        explode = (0.1, 0.1, 0.1, 0, 0, 0)
        #ax1.pie(sizes, explode=explode, labels=labels, colors = colors, autopct='%1.3f%%', shadow=False, startangle=90)
        ax1.pie(sizes, explode=explode, labels=labels, colors = colors, shadow=False, startangle=90)

    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    #plt.show()
    plt.savefig('/tmp/genotype.png')



def countTotalGenotypesForVariants(vpiDict, rareThreshold, brcaDF, hgVersion, rare):

    genotypeCounts = {'benign': {'homo':0, 'hetero': 0},
                     'pathogenic': {'homo': 0, 'hetero': 0},
                     'vus': {'homo': 0, 'hetero': 0}}

    frequenciesPerIndividual = dict()

    for individual in vpiDict:
        print(individual)
        frequenciesPerIndividual[individual] = {'benign': {'homo': 0, 'hetero': 0},
                                              'pathogenic': {'homo': 0, 'hetero': 0},
                                              'vus': {'homo': 0, 'hetero': 0}}
        for b in vpiDict[individual]['benign']:
            if b:
                if rare and getGnomadData(brcaDF, tuple(b[0]), hgVersion)[1] > rareThreshold:
                    continue
                elif b[1] == '3':
                    genotypeCounts['benign']['homo'] += 1
                    frequenciesPerIndividual[individual]['benign']['homo'] += 1
                else:
                    genotypeCounts['benign']['hetero'] += 1
                    frequenciesPerIndividual[individual]['benign']['hetero'] += 1
        for p in vpiDict[individual]['pathogenic']:
            if p:
                if rare and getGnomadData(brcaDF, tuple(p[0]), hgVersion)[1] > rareThreshold:
                    continue
                elif p[1] == '3':
                    genotypeCounts['pathogenic']['homo'] += 1
                    frequenciesPerIndividual[individual]['pathogenic']['homo'] += 1
                else:
                    genotypeCounts['pathogenic']['hetero'] += 1
                    frequenciesPerIndividual[individual]['pathogenic']['hetero'] += 1
        for v in vpiDict[individual]['vus']:
            if v:
                if rare and getGnomadData(brcaDF, tuple(v[0]), hgVersion)[1] > rareThreshold:
                    continue
                elif v[1] == '3':
                    genotypeCounts['vus']['homo'] += 1
                    frequenciesPerIndividual[individual]['vus']['homo'] += 1
                else:
                    genotypeCounts['vus']['hetero'] += 1
                    frequenciesPerIndividual[individual]['vus']['hetero'] += 1

    print(genotypeCounts)
    return genotypeCounts, frequenciesPerIndividual

if __name__ == "__main__":
    main()