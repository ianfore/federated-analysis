import sys
import json
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import pandas as pd
import logging
import time
from multiprocessing import Process, Queue, cpu_count
import math
import scipy.stats as stats
import subprocess


coordinateColumnBase = 'Genomic_Coordinate_hg'
brcaFileName = '/Users/jcasaletto/PycharmProjects/BIOBANK/federated-analysis/data/brca-variants.tsv'
hgVersion = 38

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)



def main():
    if len(sys.argv) != 10:
        print('13-out.json 13-vpi.json 13-ipv.json brca-variants.tsv brca-regions.json ancestries.json /tmp/myout 16 filter-freq')
        sys.exit(1)

    variantsFileName =sys.argv[1]
    vpiFileName = sys.argv[2]
    ipvFileName = sys.argv[3]
    brcaFileName = sys.argv[4]
    regionsFileName = sys.argv[5]
    ancestriesFileName = sys.argv[6]
    outputDir = sys.argv[7]
    numProcesses = int(sys.argv[8])
    frequencyFilter = float(sys.argv[9])

    logger.info('reading data from ' + vpiFileName)
    with open(vpiFileName, 'r') as f:
        vpiDict = json.load(f)
    f.close()
    vpiDF = pd.DataFrame(vpiDict)

    logger.info('reading data from ' + ipvFileName)
    with open(ipvFileName, 'r') as f:
        ipvDict = json.load(f)
    f.close()

    logger.info('finding variants from ' + brcaFileName)
    brcaDF = findVariantsInBRCA(brcaFileName)

    with open(variantsFileName, 'r') as f:
        variantsDict = json.load(f)
    f.close()

    logger.info('reading data from ' + regionsFileName)
    with open(regionsFileName, 'r') as f:
        regionsDict = json.load(f)
    f.close()

    logger.info('reading data from ' + ancestriesFileName)
    with open(ancestriesFileName, 'r') as f:
        ancestriesDict = json.load(f)
    f.close()

    logger.info('counting genotypes for variants on ' + str(numProcesses) + ' processes')
    t = time.time()
    q1 = Queue()
    q2 = Queue()
    processList = list()
    if frequencyFilter == 0.0:
        rare = False
    else:
        rare = True
    for i in range(numProcesses):
        p = Process(target=countTotalGenotypesForVariants,
                    args=(q1, q2, vpiDF, frequencyFilter, brcaDF, ancestriesDict, hgVersion, rare, i, numProcesses,))
        p.start()
        processList.append(p)
    logger.info('joining results from forked threads')
    genotypeCounts = dict()
    frequenciesPerIndividual = dict()
    for i in range(numProcesses):
        genotypeCounts.update(q1.get())
        frequenciesPerIndividual.update(q2.get())
    for i in range(numProcesses):
        processList[i].join()
    logger.debug('elapsed time in countTotalGenotypesForVariants() ' + str(time.time() - t))
    genotypeCountsFileName = 'genotypeCounts.json'
    logger.info('saving to ' + outputDir + '/' + genotypeCountsFileName)
    with open(outputDir + '/' + genotypeCountsFileName, 'w') as f:
        json.dump(genotypeCounts, f)
    f.close()
    plotGenotypeCounts(genotypeCounts, False, outputDir)
    fpiFileName = 'fpi.json'
    logger.info('saving to ' + outputDir + '/' + fpiFileName)
    with open(outputDir + '/' + fpiFileName, 'w') as f:
        json.dump(frequenciesPerIndividual, f)
    f.close()
    plotFrequenciesPerIndividual(frequenciesPerIndividual, outputDir)

    '''logger.info('saving to ' + outputDir + '/' + 'topmedNotGnomad.json')
    with open (outputDir + '/' + 'topmedNotGnomad.json', 'w') as f:
        json.dump(topmedNotGnomad)
    f.close()'''

    variantCounts = countTotalVariants(vpiDict)
    print('benign counts: ' + str(len(variantCounts['benign'])))
    print('pathogenic counts: ' + str(len(variantCounts['pathogenic'])))
    print('vus counts: ' + str(len(variantCounts['vus'])))


    plotVUSByPosition(variantsDict, outputDir)

    plotVUSByFrequency(variantsDict, 'maxPopFreq', outputDir)
    plotVUSByFrequency(variantsDict, 'cohortFreq', outputDir)


    ipvDict = getHardyWeinbergStats(vpiDict, variantsDict, ipvDict)

    if frequencyFilter != 0:
        ipvDict = filterOnFrequency(ipvDict, frequencyFilter)

    '''iphv = findIndividualsPerHomozygousVariant(vpiDict, variantsDict, 1.0)
    iphvFileName = 'iphv.json'
    logger.info('saving to ' + outputDir + '/' + iphvFileName)
    with open(outputDir + '/' + iphvFileName, 'w') as f:
        json.dump(iphv, f)
    f.close()'''

    '''inCIdomain = findRegionPerVariant(variantsDict, regionsDict)
    domainFileName = 'domain.json'
    logger.info('saving to ' + outputDir + '/' + domainFileName)
    with open(outputDir + '/' + domainFileName, 'w') as f:
        json.dump(inCIdomain, f)
    f.close()
    plotRegionsPerVariant(inCIdomain, outputDir)'''

    # write ipv dict back out now that it has f-value
    ipvOut = ipvFileName.replace('.json', '') + '-f.json'
    with open(ipvOut, 'w') as f:
        json.dump(ipvDict, f)
    f.close()


def filterOnFrequency(ipvDict, frequencyFilter):
    filteredIPVDict = dict()
    for v in ipvDict:
        if ipvDict[v]['maxFreq'] <= frequencyFilter or ipvDict[v]['cohortFreq'] <= frequencyFilter:
            filteredIPVDict[v] = ipvDict[v]
    return filteredIPVDict



def plotRegionsPerVariant(inCIdomain, outputDir):
    '''"(17, 43048687, 'GT', 'G')": [
        "huntsman",
        "brct"
    ],
    "(17, 43109258, 'A', 'G')": [
        "huntsman",
        "ring"
    ],'''

    # get labels
    counts = dict()
    for v in inCIdomain:
        if inCIdomain[v] != None:
            org = inCIdomain[v][0]
            code = inCIdomain[v][1]
            label = org + '-' + code
        else:
            label = 'None'
        if label not in counts:
            counts[label] = 0
        counts[label] += 1

    # Pie chart, where the slices will be ordered and plotted counter-clockwise:
    plt.rcParams['font.size'] = 12
    fig1, ax1 = plt.subplots()

    labels = counts.keys()
    sizes = list()
    for l in labels:
        sizes.append(counts[l])
    print(sizes)
    ax1.pie(sizes, labels=labels, shadow=False, startangle=90)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    #plt.show()
    plt.savefig(outputDir + '/regionsPerVariant.png')

def findRegionPerVariant(variantsDict, regionsDict):
    inCIdomain = dict()
    for v in variantsDict['homozygous vus']:
        chr = int(eval(v)[0])
        pos = int(eval(v)[1])
        if chr == 13:
            gene = 'brca2'
        elif chr == 17:
            gene = 'brca1'
        else:
            print('unrecognized chr: ' + str(chr))
            return None
        geneCIDomains = gene + 'CIDomains'
        for org in regionsDict[geneCIDomains]:
            for domain in regionsDict[geneCIDomains][org]['domains']:
                name = domain['name']
                if gene == 'brca2':
                    start = int(domain['start'])
                    end = int(domain['end'])
                else:
                    start = int(domain['end'])
                    end = int(domain['start'])
                if pos >= start and pos <= end:
                    inCIdomain[v] = (regionsDict[geneCIDomains][org]['code'], name)
                else:
                    inCIdomain[v] = None
    return inCIdomain

def findIndividualsPerHomozygousVariant(vpiDict, variantsDict, rareCutoffFrequency):
    individualsPerHomozygousVariant = defaultdict(list)
    # look at each homo vus
    for homoVUS in variantsDict['homozygous vus']:
        logger.info('looking at vus ' + str(homoVUS))
        maxPopFreq = variantsDict['homozygous vus'][homoVUS]['maxPopFreq']
        if  maxPopFreq > rareCutoffFrequency:
            continue
        for individual in vpiDict:
            # find the individuals who have expressed this as a homo vus
            for v in vpiDict[individual]['vus']:
                if tuple(v[0]) == eval(homoVUS) and v[1] == '3':
                    individualsPerHomozygousVariant[str(tuple(v[0]))].append(individual)
                    break

    return individualsPerHomozygousVariant

def generateListsForCorrelation(fpv, list_1_string, list_2_string):
    l1 = list()
    l2 = list()
    for v in fpv:
        if fpv[v]['homozygous vus']:
            if list_1_string == 'max' or list_1_string == 'min':
                l1.append(fpv[v][list_1_string]['frequency'])
            elif list_1_string == 'cohort':
                l1.append(fpv[v]['cohort'])
            else:
                # ethnicity
                maxEthnic = max(fpv[v][list_1_string]['Allele_frequency_exome_'+list_1_string+'_GnomAD'],
                              fpv[v][list_1_string]['Allele_frequency_genome_'+list_1_string+'_GnomAD'])
                l1.append(maxEthnic)
            if list_2_string == 'max' or list_2_string == 'min':
                l2.append((fpv[v][list_2_string]['frequency']))
            elif list_2_string == 'cohort':
                l2.append(fpv[v]['cohort'])
            else:
                # ethnicity
                maxEthnic = max(fpv[v][list_2_string]['Allele_frequency_exome_' + list_2_string + '_GnomAD'],
                              fpv[v][list_2_string]['Allele_frequency_genome_' + list_2_string + '_GnomAD'])
                l2.append(maxEthnic)

        else:
            continue

    # now clean up the lists -- there may be some [] from maxEthnic. Remove the corresponding elts from both lists
    # you can't remove elts from a list while you are iterating on it :) so add indices to list and remove them later
    indicesToRemove = list()
    for i in range(len(l1)):
        if l1[i] == [] or l1[i] == '-':
            indicesToRemove.append(i)
        elif l2[i] == [] or l2[i] == '-':
            indicesToRemove.append(i)

    l1_stripped = list()
    l2_stripped = list()
    for i in range(len(l1)):
        # here's where you 'remove' an elt at an index by not adding it to the final list :)
        if i not in indicesToRemove:
            # by the time you get here, you'll have 1-elt lists, so just extract the [0] value for ethnic freq
            if type(l1[i]) is list:
                if type(l1[i][0]) is str:
                    try:
                        x = eval(l1[i][0])
                    except Exception as e:
                        print('exception : ' + str(Exception) + 'x = ' + str(x))
                    l1_stripped.append(x)
                else:
                    l1_stripped.append(l1[i][0])
            else:
                l1_stripped.append(l1[i])
            if type(l2[i]) is list:
                if type(l2[i][0]) is str:
                    try:
                        x = eval(l2[i][0])
                    except Exception as e:
                        print('exception : ' + str(Exception) + 'x = ' + str(x))
                    l2_stripped.append(x)
                else:
                    l2_stripped.append(l2[i][0])
            else:
                l2_stripped.append(l2[i])

    return l1_stripped, l2_stripped

def mean(someList):
    total = 0
    for a in someList:
        total += float(a)
    mean = total/len(someList)
    return mean

def standDev(someList):
    listMean = mean(someList)
    dev = 0.0
    for i in range(len(someList)):
        dev += (someList[i]-listMean)**2
    dev = dev**(1/2.0)
    return dev

def pearsonCorrelation(someList1, someList2):

    # First establish the means and standard deviations for both lists.
    xMean = mean(someList1)
    yMean = mean(someList2)
    xStandDev = standDev(someList1)
    yStandDev = standDev(someList2)
    # r numerator
    rNum = 0.0
    for i in range(len(someList1)):
        rNum += (someList1[i]-xMean)*(someList2[i]-yMean)

    # r denominator
    rDen = xStandDev * yStandDev

    r =  rNum/rDen
    return r

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


def findVariantsInBRCA(fileName):
    return pd.read_csv(fileName, sep='\t', header=0, dtype=str)

def getFrequenciesPerVariant(variantsDict, brcaDF, hgVersion, ethnicity1, ethnicity2):
    freqPerVariant = defaultdict(dict)
    # get min, max, median, and #hom and #hem

    for variant in variantsDict['cooccurring vus']:
        #logger.debug('variant: ' + str(variant))
        freqPerVariant[variant]['cooccurring vus'] = True
        freqPerVariant[variant]['homozygous vus'] = False
        gnomadData1 = getGnomadData(brcaDF, variant, hgVersion, ethnicity1)
        for g in gnomadData1:
            freqPerVariant[variant][g] = gnomadData1[g]
        if ethnicity2:
            gnomadData2 = getGnomadData(brcaDF, variant, hgVersion, ethnicity2)
            for g in gnomadData2:
                freqPerVariant[variant][g] = gnomadData2[g]
        freqPerVariant[variant]['cohort'] = variantsDict['cooccurring vus'][variant]['allele frequencies']['cohortFreq']
    for variant in variantsDict['homozygous vus']:
        #logger.debug('variant: ' + str(variant))
        freqPerVariant[variant]['homozygous vus'] = True
        if 'cooccurring vus' not in freqPerVariant[variant]:
            freqPerVariant[variant]['cooccurring vus'] = False
            gnomadData1 = getGnomadData(brcaDF, variant, hgVersion, ethnicity1)
            for g in gnomadData1:
                freqPerVariant[variant][g] = gnomadData1[g]
            if ethnicity2:
                gnomadData2 = getGnomadData(brcaDF, variant, hgVersion, ethnicity2)
                for g in gnomadData2:
                    freqPerVariant[variant][g] = gnomadData2[g]
            freqPerVariant[variant]['cohort'] = variantsDict['homozygous vus'][variant]['cohortFreq']

    return freqPerVariant

def getGnomadData(brcaDF, vus, hgVersion, ethnicity):
    # TODO write a unit test
    # 13:g.32393468:C>CT
    #hgString = 'chr' + str(vus[0][0]) + ':g.' + str(vus[0][1]) + ':' + str(vus[0][2]) + '>' + str(vus[0][3])
    if type(vus) is str:
        vus = eval(vus)
    hgString = 'chr' + str(vus[0]) + ':g.' + str(vus[1]) + ':' + str(vus[2]) + '>' + str(vus[3])

    # first, get list of columns for GnomAD allleles
    #alleleFrequencyPrefixes = ['Allele_frequency_genome_', 'Allele_frequency_']
    # 'Allele_count_exome_AMR_GnomAD', 'Allele_count_hemi_exome_AMR_GnomAD', 'Allele_count_hom_exome_AMR_GnomAD',
    # 'Allele_number_exome_AMR_GnomAD', 'Allele_frequency_exome_AMR_GnomAD' ,'Allele_frequency_AMR_GnomAD'

    gnomad = [v for v in list(brcaDF.columns) if 'GnomAD' in v]

    # second, get frequencies across exomes and genomes to determine max
    # return population, frequency, count, and number
    # replace "frequency" with "count" and "number" in Allele_frequency_genome_AFR_GnomAD

    allDict = dict()

    alleleFrequencies = [v for v in gnomad if 'Allele_frequency' in v]
    allDict['max'] = getMaxGnomad(brcaDF, hgString, hgVersion, alleleFrequencies)
    allDict['min'] = getMinGnomad(brcaDF, hgString, hgVersion, alleleFrequencies)

    if ethnicity is not None:
        alleleFrequencies = [v for v in gnomad if ethnicity in v]
        x =  getPopulationGnomadData(brcaDF, hgString, hgVersion, alleleFrequencies)
        allDict[ethnicity] = x
    return allDict

def getPopulationGnomadData(brcaDF, hgString, hgVersion, alleleFrequencies):
    ethData = dict()
    alleleFrequencies2 = [v for v in alleleFrequencies if 'Allele_frequency' in v]
    for af in alleleFrequencies2:
        ethData[af] = brcaDF[brcaDF[coordinateColumnBase + str(hgVersion)] == hgString][af].tolist()

    # alleleFrequencies2 = [Allele_frequency_genome_NFE_GnomAD, Allele_frequency_exome_NFE_GnomAD, Allele_frequency_NFE_GnomAD]
    # we want just the genome which is [0] and that is a list of one element, so we want [0] of that so [0][0]

    if not ethData:
        return 0.0

    genomeFreq = alleleFrequencies2[0]

    if not ethData[genomeFreq]:
        return 0.0

    if ethData[genomeFreq][0] == math.nan or ethData[genomeFreq][0] == '-':
        return 0.0
    else:
        return float(ethData[genomeFreq][0])
    '''maxFreq = 0.0
    for af in ethData:
        if ethData[af][0] == math.nan or ethData[af][0] == '-':
            continue
        elif float(ethData[af][0]) > maxFreq:
            maxFreq = float(ethData[af][0])
    return maxFreq'''

def getMinGnomad(brcaDF, hgString, hgVersion, alleleFrequencies):
    minData = {'frequency': 1.1, 'population': None}
    for af in alleleFrequencies:
        freq = 0.0
        alleleFreqList = brcaDF[brcaDF[coordinateColumnBase + str(hgVersion)] == hgString][af].tolist()
        if alleleFreqList:
            try:
                if len(alleleFreqList) > 1:
                    print('more than one')
                # yes, it always returns a list of length 0
                freq = float(alleleFreqList[0])
            except ValueError:
                continue
            if freq < minData['frequency'] and freq != 0.0:
                # TODO get homo and abs counts as well
                minData['frequency'] = freq
                minData['population'] = af
    if minData['frequency'] == 1.1:
        minData['population'] = 0.0
    return minData

def getMaxGnomad(brcaDF, hgString, hgVersion, alleleFrequencies):
    maxData = {'frequency': 0.0, 'population': None}
    for af in alleleFrequencies:
        freq = 0.0
        alleleFreqList = brcaDF[brcaDF[coordinateColumnBase + str(hgVersion)] == hgString][af].tolist()
        if alleleFreqList:
            try:
                if len(alleleFreqList) > 1:
                    print('more than one')
                # yes, it always returns a list of length 0
                freq = float(alleleFreqList[0])
            except ValueError:
                continue
            if freq > maxData['frequency']:
                # TODO get homo and abs counts as well
                maxData['frequency'] = freq
                maxData['population'] = af
    return maxData

def getFisherExact(results, key, allValues):
    # create 2x2 contingency table
    a = results[key]['with'][allValues[0]]
    b = results[key]['without'][allValues[0]]
    c = results[key]['with'][allValues[1]]
    d = results[key]['without'][allValues[1]]

    # run fisher exact test
    oddsRatio, pValue = stats.fisher_exact([[a, b], [c, d]])

    # get confidence interval for odds ratio: CI = e^(ln(OR) +/- [1.96 * sqrt(1/a + 1/b + 1/c + 1/d)])
    '''logOR = log(oddsRatio)
    SE_logOR = sqrt(1/a + 1/b + 1/c + 1/d)
    lowerBound = logOR - 1.96 * SE_logOR
    upperBound = logOR + 1.96 * SE_logOR

    CI95_lower = exp(lowerBound)
    CI95_upper = exp(upperBound)'''

    # use R fisher.test to get confidence interval for LOR
    row1 = str('c(' + str(a) + ',' + str(b) + ')')
    row2 = str('c(' + str(c) + ',' + str(d) + ')')
    CIcmd_1 = str('/usr/bin/Rscript -e "fisher.test(rbind(' + row1 + ',' +  row2 + '))"')
    CIcmd_2 = str('| grep -A1 confidence | sed -n "2,2 p" | xargs')
    CIcmd = str(CIcmd_1 + CIcmd_2)
    lb, ub = subprocess.check_output(CIcmd, shell=True).split()

    # cast result = b'0.002439905 19.594803004\n' to floats
    lb = float(lb)
    ub = float(ub)

    # insert OR and p-value into results dict
    results[key]['OR'] = round(oddsRatio, 3)
    results[key]['P value'] =  round(pValue, 3)
    results[key]['95% CI'] = (round(lb, 2), round(ub, 2))

def getHardyWeinbergStats(vpiDict, variantsDict, individualsPerVariant):
    bVars, pVars, vVars = calculateZygosityFrequenciesPerVariant(vpiDict)
    #bVars, pVars, vVars = hardyWeinbergChiSquareTest(bVars, pVars, vVars, len(vpiDict))
    bVars, pVars, vVars = hardyWeinbergFisherExactTest(bVars, pVars, vVars, len(vpiDict), 0.05)
    bVars, pVars, vVars = hardyWeinbergStatistics(bVars, pVars, vVars, significance=0.01)
    rejectHW = {'benign': 0, 'pathogenic': 0, 'vus': 0}
    acceptHW = {'benign': 0, 'pathogenic': 0, 'vus': 0}
    acceptF = {'benign': 0, 'pathogenic': 0, 'vus': 0}
    rejectF = {'benign': 0, 'pathogenic': 0, 'vus': 0}

    for b in bVars:
        if bVars[b]['accept hw'] is False:
            rejectHW['benign'] += 1
        else:
            acceptHW['benign'] += 1
        if bVars[b]['accept F'] is False:
            rejectF['benign'] += 1
        else:
            acceptF['benign'] += 1
        individualsPerVariant[str(b)]['accept F'] = bVars[b]['accept F']
        individualsPerVariant[str(b)]['F'] = bVars[b]['F']
        individualsPerVariant[str(b)]['Z'] = bVars[b]['Z']
        individualsPerVariant[str(b)]['class'] = 'benign'
        individualsPerVariant[str(b)]['aa'] = bVars[b]['aa']
        individualsPerVariant[str(b)]['Aa'] = bVars[b]['Aa']
        individualsPerVariant[str(b)]['AA'] = bVars[b]['AA']


    for p in pVars:
        if pVars[p]['accept hw'] is False:
            rejectHW['pathogenic'] += 1
        else:
            acceptHW['pathogenic'] += 1
        if pVars[p]['accept F'] is False:
            rejectF['pathogenic'] += 1
        else:
            rejectF['pathogenic'] += 1
        individualsPerVariant[str(p)]['F'] = pVars[p]['F']
        individualsPerVariant[str(p)]['Z'] = pVars[p]['Z']
        individualsPerVariant[str(p)]['accept F'] = pVars[p]['accept F']
        individualsPerVariant[str(p)]['class'] = 'pathogenic'
        individualsPerVariant[str(p)]['aa'] = pVars[p]['aa']
        individualsPerVariant[str(p)]['Aa'] = pVars[p]['Aa']
        individualsPerVariant[str(p)]['AA'] = pVars[p]['AA']




    rejectVUS = {'cooccurring vus': 0, 'homozygous vus': 0}
    acceptVUS = {'cooccurring vus': 0, 'homozygous vus': 0}
    rejectVUS_F = {'cooccurring vus': 0, 'homozygous vus': 0}
    acceptVUS_F = {'cooccurring vus': 0, 'homozygous vus': 0}

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

        if vVars[v]['accept F'] is False:
            rejectF['vus'] += 1
            if str(v) in variantsDict['cooccurring vus']:
                rejectVUS_F['cooccurring vus'] += 1
            if str(v) in variantsDict['homozygous vus']:
                rejectVUS_F['homozygous vus'] += 1
        else:
            acceptF['vus'] += 1
            if str(v) in variantsDict['cooccurring vus']:
                acceptVUS_F['cooccurring vus'] += 1
            if str(v) in variantsDict['homozygous vus']:
                acceptVUS_F['homozygous vus'] += 1
        individualsPerVariant[str(v)]['F'] = vVars[v]['F']
        individualsPerVariant[str(v)]['Z'] = vVars[v]['Z']
        individualsPerVariant[str(v)]['accept F'] = vVars[v]['accept F']
        individualsPerVariant[str(v)]['class'] = 'vus'
        individualsPerVariant[str(v)]['aa'] = vVars[v]['aa']
        individualsPerVariant[str(v)]['Aa'] = vVars[v]['Aa']
        individualsPerVariant[str(v)]['AA'] = vVars[v]['AA']

    # check to see if 654 vus that reject HW are same vus that reject F
    vusRejectingBothHWandF = list()
    for v in vVars:
        if str(v) not in variantsDict['homozygous vus']:
            continue
        elif vVars[v]['accept hw'] == False and vVars[v]['accept F'] == False:
            vusRejectingBothHWandF.append(v)


    print('reject HW: ' + str(rejectHW))
    print('accept HW: ' + str(acceptHW))
    print('accept F:  ' + str(acceptF))
    print('reject F:  ' + str(rejectF))
    print('num co-occurring vus that reject HW: ' + str(rejectVUS['cooccurring vus']))
    print('num co-occurring vus that accept HW:' + str(acceptVUS['cooccurring vus']))
    print('num homozygous vus that reject HW: ' + str(rejectVUS['homozygous vus']))
    print('num homozygous vus that accept HW: ' + str(acceptVUS['homozygous vus']))
    print('num co-occurring vus that reject F: ' + str(rejectVUS_F['cooccurring vus']))
    print('num co-occurring vus that accept F:' + str(acceptVUS_F['cooccurring vus']))
    print('num homozygous vus that reject F: ' + str(rejectVUS_F['homozygous vus']))
    print('num homozygous vus that accept F: ' + str(acceptVUS_F['homozygous vus']))

    print('list of vus rejecting both HW and F: ' + str(vusRejectingBothHWandF))
    print('length list of vus rejecting both HW and F: ' + str(len(vusRejectingBothHWandF)))

    return individualsPerVariant

def calculateZygosityFrequenciesPerVariant(vpiDict):
    benignVariants = dict()
    pathogenicVariants = dict()
    vusVariants = dict()
    for individual in vpiDict:
        # if it's in list of variants for individual, then it must be one of 1|1 (3), 0|1 (1), or 1|0 (2)
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

def hardyWeinbergStatistics(bVars, pVars, vVars, significance):
    # https://en.wikipedia.org/wiki/Hardy Weinberg_principle
    # degrees of freedom = 1
    # The inbreeding coefficient, F (see also F-statistics), is one minus the observed frequency of
    # heterozygotes over that expected from Hardy Weinberg equilibrium.
    # F = [ E(f(Aa)) - O(f(Aa)) ] / [ E(f(Aa)) ] = 1 - O(f(Aa))/E(f(Aa))
    # where E(f(Aa)) = 2pq
    # For two alleles, the chi - squared goodness of fit test for Hardy Weinberg proportions is equivalent to the test
    # for inbreeding, F = 0.
    # The inbreeding coefficient is unstable as the expected value approaches zero, and thus not useful for rare and
    # very common alleles. For: E = 0, O > 0, F = inf and E = 0, O = 0, F is undefined.

    # from the Assortative Mating study on GALA data:
    # Allele frequency differences between groups were calculated
    # using standard chi-square tests. We tested for Hardy Weinberg
    # equilibrium at marker loci by using the Z-statistic
    # Z = (((n2 x n0) - n1^2 ) / (n2 + n1/2)(n0 + n1/2) ) N^1/2
    # where n2 and n0 are the number of homozygotes and n1 the
    # number of heterozygotes observed; N = n2 + n1 + n0. Under
    # the null hypothesis of no within-locus allelic correlation, Z
    # has a normal distribution with mean 0 and variance 1. We
    # chose to use a one-sided test as opposed to a two-sided chisquare
    # test because we specifically were searching for an
    # excess of homozygotes, as predicted by assortative mating.

    for b in bVars:
        # 1. get E(f(Aa))
        bVars[b]['E(f(Aa))'] = 2 * bVars[b]['p'] * bVars[b]['q']
        # 2. O(f(Aa)) = bVars[b]['Aa']
        # 3. calculate F
        try:
            bVars[b]['F'] = ( bVars[b]['E(f(Aa))'] - bVars[b]['Aa'] ) / ( bVars[b]['E(f(Aa))'] )
            n0 = bVars[b]['aa']
            n2 = bVars[b]['AA']
            n1 = bVars[b]['Aa']
            N = n0 + n1 + n2
            bVars[b]['Z'] = ( (n2 * n0 - n1**2) / ((n2 + 0.5*n1) * (n0 + 0.5*n1)) ) * math.sqrt(N)


        except Exception as e:
            print('exception calculating F: ' + str(e))
            continue
        if abs(bVars[b]['F']) <= significance:
            bVars[b]['accept F'] = True
        else:
            bVars[b]['accept F'] = False

    for p in pVars:
        # 1. get E(f(Aa))
        pVars[p]['E(f(Aa))'] = 2 * pVars[p]['p'] * pVars[p]['q']
        # 2. O(f(Aa)) = pVars[p]['Aa']
        # 3. calculate F
        try:
            pVars[p]['F'] = ( pVars[p]['E(f(Aa))'] - pVars[p]['Aa'] ) / ( pVars[p]['E(f(Aa))'] )
            n0 = pVars[p]['aa']
            n2 = pVars[p]['AA']
            n1 = pVars[p]['Aa']
            N = n0 + n1 + n2
            pVars[p]['Z'] = ((n2 * n0 - n1 ** 2) / ((n2 + 0.5 * n1) * (n0 + 0.5 * n1))) * math.sqrt(N)
        except Exception as e:
            print('exception calculating F: ' + str(e))
            continue
        if abs(pVars[p]['F']) <= significance:
            pVars[p]['accept F'] = True
        else:
            pVars[p]['accept F'] = False


    for v in vVars:
        # 1. get E(f(Aa))
        vVars[v]['E(f(Aa))'] = 2 * vVars[v]['p'] * vVars[v]['q']
        # 2. O(f(Aa)) = vVars[v]['Aa']
        # 3. calculate F
        try:
            vVars[v]['F'] = ( vVars[v]['E(f(Aa))'] - vVars[v]['Aa'] ) / ( vVars[v]['E(f(Aa))'] )
            n0 = vVars[v]['aa']
            n2 = vVars[v]['AA']
            n1 = vVars[v]['Aa']
            N = n0 + n1 + n2
            vVars[v]['Z'] = ((n2 * n0 - n1 ** 2) / ((n2 + 0.5 * n1) * (n0 + 0.5 * n1))) * math.sqrt(N)
        except Exception as e:
            print('exception calculating F: ' + str(e))
            vVars[v]['F'] = 0.0
        if abs(vVars[v]['F']) <= significance:
            vVars[v]['accept F'] = True
        else:
            vVars[v]['accept F'] = False

    return bVars, pVars, vVars


def hardyWeinbergFisherExactTest(bVars, pVars, vVars, N, pValue):
    # http://europepmc.org/backend/ptpmcrender.fcgi?accid=PMC1199378&blobtype=pdf
    # p(NAB = nAB|N,nA) = [ (2**nAB x N!) / (nAA!nAB!nBB!) ] x (nA! x nB!)/(2N)!
    # let W = (2**nAB x N!)
    # let X = (nAA!nAB!nBB!)
    # let Y = nA! x nB!
    # let Z = (2N)!
    # then p(NAB = nAB|N,nA) = (W / X) x (Y / Z)

    for b in bVars:
        nA = abs(2 * bVars[b]['aa'] - bVars[b]['Aa'])
        nB = abs(2 * bVars[b]['AA'] - bVars[b]['Aa'])
        nAA = bVars[b]['aa']
        nBB = bVars[b]['AA']
        nAB = bVars[b]['Aa']
        # calculate W
        W = (2**nAB) * math.factorial(N)

        # calculate X
        X = math.factorial(nAA) * math.factorial(nAB) * math.factorial(nBB)

        # calculate Y
        Y = math.factorial(nA) * math.factorial(nB)

        # calculate Z
        Z = math.factorial(2 * N)

        # calculate p
        bVars[b]['fisher'] = (W / X) * (Y / Z)

        # 8. compare against p-value for 1 degree of freedom at 0.05 significance (3.84)
        if bVars[b]['fisher'] >= pValue:
            bVars[b]['accept hw'] = False
        else:
            bVars[b]['accept hw'] = True

    for p in pVars:
        nA = abs(2 * pVars[p]['aa'] - pVars[p]['Aa'])
        nB = abs(2 * pVars[p]['AA'] - pVars[p]['Aa'])
        nAA = pVars[p]['aa']
        nBB = pVars[p]['AA']
        nAB = pVars[p]['Aa']
        # calculate W
        W = (2 ** nAB) * math.factorial(N)

        # calculate X
        X = math.factorial(nAA) * math.factorial(nAB) * math.factorial(nBB)

        # calculate Y
        Y = math.factorial(nA) * math.factorial(nB)

        # calculate Z
        Z = math.factorial(2 * N)

        # calculate p
        pVars[p]['fisher'] = (W / X) * (Y / Z)

        # 8. compare against p-value for 1 degree of freedom at 0.05 significance (3.84)
        if pVars[p]['fisher'] >= pValue:
            pVars[p]['accept hw'] = False
        else:
            pVars[p]['accept hw'] = True

    for v in vVars:

        nA = abs(2 * vVars[v]['aa'] - vVars[v]['Aa'])
        nB = abs(2 * vVars[v]['AA'] - vVars[v]['Aa'])
        nAA = vVars[v]['aa']
        nBB = vVars[v]['AA']
        nAB = vVars[v]['Aa']
        # calculate W
        W = (2 ** nAB) * math.factorial(N)

        # calculate X
        X = math.factorial(nAA) * math.factorial(nAB) * math.factorial(nBB)

        # calculate Y
        Y = math.factorial(nA) * math.factorial(nB)

        # calculate Z
        Z = math.factorial(2 * N)

        # calculate p
        vVars[v]['fisher'] = (W / X) * (Y / Z)

        # 8. compare against p-value for 1 degree of freedom at 0.05 significance (3.84)
        if vVars[v]['fisher'] >= pValue:
            vVars[v]['accept hw'] = False
        else:
            vVars[v]['accept hw'] = True

    return bVars, pVars, vVars

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
        '''if expAA == 0 or expAa == 0 or expaa == 0:
            pVars[p]['chisquare'] = 0
        else:
            pVars[p]['chisquare'] = 1.0/expAA * (pVars[p]['AA'] - expAA)**2 + \
                                    1.0/expAa * (pVars[p]['Aa'] - expAa)**2 + \
                                    1.0/expaa * (pVars[p]['aa'] - expaa)**2'''
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

def plotVUSByFrequency(variantsDict, freq, outputDir):
    homoVUS = variantsDict['homozygous vus']
    popFreqs = list()
    for vus in homoVUS:
        popFreqs.append(homoVUS[vus][freq])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    bp = ax.boxplot(popFreqs)
    plt.xlabel(freq)
    #plt.show()
    binList = [0.00001, 0.0001, 0.001, 0.01, 0.1]
    binPlot(popFreqs, 25, freq, "number of homo VUS", float, 3, binList, outputDir, 'homoVUSFreqs.png')

    popFreqs = list()
    heteroVUS = variantsDict['cooccurring vus']
    for vus in heteroVUS:
        if vus in homoVUS:
            continue
        elif freq is 'maxPopFreq':
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
    #plt.show()
    binList = [0.00001, 0.0001, 0.001, 0.01, 0.1]
    if len(popFreqs) == 0:
        print('empty list')
    else:
        binPlot(popFreqs, 25, freq, "number of hetero VUS", float, 3, binList, outputDir, 'heteroVUSFreqs.png')

def plotVUSByPosition(variantsDict, outputDir):
    locations = list()
    homozygousVUS = variantsDict['homozygous vus']
    for vus in homozygousVUS:
        locations.append(int(eval(vus)[1]))
    if len(locations) == 0:
        print('empty list')
    else:
        binPlot(locations, 10000, "chromosome position bins", "number of homo VUS", int, 0, None, outputDir, 'homoVUSPositions.png')

    locations = list()
    heterozygousVUS = variantsDict['cooccurring vus']
    for vus in heterozygousVUS:
        if vus in homozygousVUS:
            continue
        else:
            locations.append(int(eval(vus)[1]))
    if len(locations) == 0:
        print('empty list')
    else:
        binPlot(locations, 10000, "chromosome position bins", "number of het VUS", int, 0, None, outputDir, 'heteroVUSPositions.png')

def plotZygosityRatiosPerIndividual(frequenciesPerIndividual, outputDir):
    benignRatios = list()
    pathogenicRatios = list()
    vusRatios = list()
    for individual in frequenciesPerIndividual:
        f1 = frequenciesPerIndividual[individual]['benign']['homo']
        f2 = frequenciesPerIndividual[individual]['benign']['hetero']
        if f2 != 0:
            benignRatios.append(f1 / f2)
        f1 = frequenciesPerIndividual[individual]['pathogenic']['homo']
        f2 = frequenciesPerIndividual[individual]['pathogenic']['hetero']
        if f2 != 0:
            pathogenicRatios.append(f1 / f2)
        f1 = frequenciesPerIndividual[individual]['vus']['homo']
        f2 = frequenciesPerIndividual[individual]['vus']['hetero']
        if f2 != 0:
            vusRatios.append(f1 / f2)
    binPlot(benignRatios, 10, "homozygous benign homo:hetero ratio", "number of individuals", int, 0, None,
            outputDir, 'homoBenignRatios.png')
    binPlot(pathogenicRatios, 10, "homozygous pathogenic homo:hetero ratio", "number of individuals", int, 0, None,
            outputDir, 'homoPathogenicRatios.png')
    binPlot(vusRatios, 10, "homozygous vus homo:hetero ratio", "number of individuals", int, 0, None,
            outputDir, 'homoVUSRatios.png')


def plotFrequenciesPerIndividual(frequenciesPerIndividual, outputDir):
    # count the number of individuals per frequency
    homo_ben_counts = list()
    homo_path_counts = list()
    homo_vus_counts = list()
    hetero_ben_counts = list()
    hetero_path_counts = list()
    hetero_vus_counts = list()


    for individual in frequenciesPerIndividual:
        homo_ben_counts.append(frequenciesPerIndividual[individual]['benign']['homo'])
        #hetero_ben_counts.append(frequenciesPerIndividual[individual]['benign']['hetero'])
        #homo_path_counts.append(frequenciesPerIndividual[individual]['pathogenic']['homo'])
        #hetero_path_counts.append(frequenciesPerIndividual[individual]['pathogenic']['hetero'])
        #homo_vus_counts.append(frequenciesPerIndividual[individual]['vus']['homo'])
        #hetero_vus_counts.append(frequenciesPerIndividual[individual]['vus']['hetero'])

    if len(homo_ben_counts) == 0:
        print('empty list')
    else:
        binPlot(homo_ben_counts, 10, "homozygous benign variant count bins", "number of individuals", int, 0, None, outputDir, 'homoBenignFPI.png')
    #binPlot(hetero_ben_counts, 10, "heterozygous benign variant count bins", "number of individuals", int, 0, None)
    #binPlot(homo_path_counts, 10, "homozygous pathogenic variant count bins", "number of individuals", int, 0, None)
    #binPlot(hetero_path_counts, 10, "heterozygous pathogenic variant count bins", "number of individuals", int, 0, [0, 1])
    #binPlot(homo_vus_counts, 10, "homozygous vus variant count bins", "number of individuals", int, 0, None)
    #binPlot(hetero_vus_counts, 10, "heterozygous vus variant count bins", "number of individuals", int, 0, None)

def plotGenotypeCounts(genotypeCounts, rare, outputDir):
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

    if rare:
        labels = ['Homo VUS',  'Hetero VUS']
        colors = ['green', 'brown']
        sizes = [genotypeCounts['vus']['homo'], genotypeCounts['vus']['hetero']]
        explode = (0, 0)
        ax1.pie(sizes, explode=explode, colors=colors, labels=labels, shadow=False, startangle=90)

    else:
        labels = ['Homo benign', 'Homo path', 'Homo VUS', 'Hetero benign', 'Hetero path', 'Hetero VUS']
        sizes = [genotypeCounts['benign']['homo'], genotypeCounts['pathogenic']['homo'],
                    genotypeCounts['vus']['homo'], genotypeCounts['benign']['hetero'],
                    genotypeCounts['pathogenic']['hetero'], genotypeCounts['vus']['hetero']]
        colors = ['red', 'yellow', 'green', 'orange', 'blue', 'brown']
        explode = (0.1, 0.1, 0.1, 0, 0, 0)

        patches, texts = ax1.pie(sizes, explode=explode, colors=colors, shadow=False, startangle=90)
        plt.legend(patches, labels, fontsize=8)

        #ax1.pie(sizes, explode=explode, labels=labels, colors = colors, autopct='%1.3f%%', shadow=False, startangle=90)
        ax1.pie(sizes, explode=explode, colors = colors, shadow=False, startangle=90)

    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    #plt.show()
    plt.savefig(outputDir + '/genotypeCounts.png')

def divide(n, d):
   res = list()
   qu = int(n/d)
   rm = n%d
   for i in range(d):
       if i < rm:
           res.append(qu + 1)
       else:
           res.append(qu)
   return res

def getStartAndEnd(partitionSizes, threadID):
    start = 0
    for i in range(threadID):
        start += partitionSizes[i]

    end = start + partitionSizes[threadID]

    return start, end

def countTotalGenotypesForVariants(q1, q2, vpiDF, rareThreshold, brcaDF, ancestriesDict, hgVersion, rare, threadID, numProcesses):

    genotypeCounts = {'benign': {'homo':0, 'hetero': 0},
                     'pathogenic': {'homo': 0, 'hetero': 0},
                     'vus': {'homo': 0, 'hetero': 0}}
    frequenciesPerIndividual = dict()
    #notInGnomad = {'benign': [], 'pathogenic': [], 'vus': []}
    individuals = list()
    for individual in vpiDF:
        individuals.append(individual)
    n = len(individuals)
    logger.info('total number of individuals = ' + str(n))
    partitionSizes = divide(n, numProcesses)
    start, end = getStartAndEnd(partitionSizes, threadID)

    logger.info('threadID = ' + str(threadID) + ' processing from ' + str(start) + ' to ' + str(end))
    for i in range(start, end):
        individual = individuals[i]
        ethnicity = ancestriesDict[individual]['gnomadPop']
        logger.debug(individual)
        frequenciesPerIndividual[individual] = {'benign': {'homo': 0, 'hetero': 0},
                                            'pathogenic': {'homo': 0, 'hetero': 0},
                                              'vus': {'homo': 0, 'hetero': 0}}

        for b in vpiDF[individual]['benign']:
            if b:
                if rare and getGnomadData(brcaDF, tuple(b[0]), hgVersion, None)['max']['frequency'] > rareThreshold:
                    continue
                    '''gnomadData = getGnomadData(brcaDF, tuple(b[0]), hgVersion, ethnicity)
                    if gnomadData[ethnicity] > rareThreshold:
                        continue
                    elif gnomadData[ethnicity] == 0.0:
                        notInGnomad['benign'].append(b)'''
                if b[1] == '3':
                    genotypeCounts['benign']['homo'] += 1
                    frequenciesPerIndividual[individual]['benign']['homo'] += 1
                else:
                    genotypeCounts['benign']['hetero'] += 1
                    frequenciesPerIndividual[individual]['benign']['hetero'] += 1
        for p in vpiDF[individual]['pathogenic']:
            if p:
                if rare and getGnomadData(brcaDF, tuple(p[0]), hgVersion, ethnicity)['max']['frequency'] > rareThreshold:
                    '''gnomadData = getGnomadData(brcaDF, tuple(p[0]), hgVersion, ethnicity)
                    if gnomadData[ethnicity] > rareThreshold:
                        continue
                    elif gnomadData[ethnicity] == 0.0:
                        notInGnomad['pathogenic'].append(p)'''
                    continue
                if p[1] == '3':
                    genotypeCounts['pathogenic']['homo'] += 1
                    frequenciesPerIndividual[individual]['pathogenic']['homo'] += 1
                else:
                    genotypeCounts['pathogenic']['hetero'] += 1
                    frequenciesPerIndividual[individual]['pathogenic']['hetero'] += 1
        for v in vpiDF[individual]['vus']:
            if v:
                if rare and getGnomadData(brcaDF, tuple(v[0]), hgVersion, ethnicity)['max']['frequency'] > rareThreshold:
                    '''gnomadData = getGnomadData(brcaDF, tuple(v[0]), hgVersion, ethnicity)
                    if gnomadData[ethnicity] > rareThreshold:
                        continue
                    elif gnomadData[ethnicity] == 0.0:
                        notInGnomad['vus'].append(v)'''
                    continue
                if v[1] == '3':
                    genotypeCounts['vus']['homo'] += 1
                    frequenciesPerIndividual[individual]['vus']['homo'] += 1

                else:
                    genotypeCounts['vus']['hetero'] += 1
                    frequenciesPerIndividual[individual]['vus']['hetero'] += 1

        
    q1.put(genotypeCounts)
    q2.put(frequenciesPerIndividual)
    #q3.put(notInGnomad)
    logger.debug('finished putting results in queue')

if __name__ == "__main__":
    main()