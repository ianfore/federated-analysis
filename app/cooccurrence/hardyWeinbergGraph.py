import matplotlib.pyplot as plt
import logging
import numpy as np
import pandas as pd
import sys
from scipy.stats import pearsonr
import math


logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

def sumList(a, b):
    if len(a) != len(b):
        raise Exception('lists must be same length')
    mySum = list()
    for i in range(len(a)):
        mySum.append(a[i] + b[i])
    return mySum

def divideList(a, b):
    if len(a) != len(b):
        raise Exception('lists must be same length')
    myDivide = list()
    for i in range(len(a)):
        myDivide.append((1.0 * a[i]) / (1.0 * b[i]))
    return myDivide

def main():
    if len(sys.argv) != 6:
        print('brca1-report.tsv brca2-report.tsv output-dir all|homo brca1|brca2|both')
        sys.exit(1)

    brca1_report = sys.argv[1]
    brca2_report = sys.argv[2]
    outputDir = sys.argv[3]
    subset = sys.argv[4]
    gene = sys.argv[5]

    df_1 = pd.read_csv(brca1_report, header=0, sep='\t')
    df_2 = pd.read_csv(brca2_report, header=0, sep='\t')

    outputFileName = subset + '_' + gene

    if subset == 'homo':
        samples_1 = df_1[df_1['homozygousSample'] != 'None']
        samples_2 = df_2[df_2['homozygousSample'] != 'None']
        homoAltList_1 = list(samples_1['homo_alt'])
        homoAltList_2 = list(samples_2['homo_alt'])
        homoRefList_1 = list(samples_1['homo_ref'])
        homoRefList_2 = list(samples_2['homo_ref'])
        heteroList_1 = list(samples_1['hetero'])
        heteroList_2 = list(samples_2['hetero'])

        if gene == 'both':
            pList = list(samples_1['p']) + list(samples_2['p'])
            qList = list(samples_1['q']) + list(samples_2['q'])
            homoRefList = homoRefList_1 + homoRefList_2
            homoAltList = homoAltList_1 + homoAltList_2
            heteroList = heteroList_1 + heteroList_2
            totalList = sumList(sumList(homoRefList, homoAltList), heteroList)
            hrList = divideList(homoRefList, totalList)
            haList = divideList(homoAltList, totalList)
            hList = divideList(heteroList, totalList)
        elif gene == 'brca1':
            pList = list(samples_1['p'])
            qList = list(samples_1['q'])
            homoRefList = homoRefList_1
            homoAltList = homoAltList_1
            heteroList = heteroList_1
            totalList = sumList(sumList(homoRefList, homoAltList), heteroList)
            hrList = divideList(homoRefList, totalList)
            haList = divideList(homoAltList, totalList)
            hList = divideList(heteroList, totalList)
        elif gene == 'brca2':
            pList = list(samples_2['p'])
            qList = list(samples_2['q'])
            homoRefList = homoRefList_2
            homoAltList = homoAltList_2
            heteroList = heteroList_2
            totalList = sumList(sumList(homoRefList, homoAltList), heteroList)
            hrList = divideList(homoRefList, totalList)
            haList = divideList(homoAltList, totalList)
            hList = divideList(heteroList, totalList)
        else:
            print('brca1-report.tsv brca2-report.tsv output-dir all|homo brca1|brca2|both')
            sys.exit(1)
    elif subset == 'all':
        homoAltList_1 = list(df_1['homo_alt'])
        homoAltList_2 = list(df_2['homo_alt'])
        homoRefList_1 = list(df_1['homo_ref'])
        homoRefList_2 = list(df_2['homo_ref'])
        heteroList_1 = list(df_1['hetero'])
        heteroList_2 = list(df_2['hetero'])

        if gene == 'both':
            pList = list(df_1['p']) + list(df_2['p'])
            qList = list(df_1['q']) + list(df_2['q'])
            homoRefList = homoRefList_1 + homoRefList_2
            homoAltList = homoAltList_1 + homoAltList_2
            heteroList = heteroList_1 + heteroList_2
            totalList = sumList(sumList(homoRefList, homoAltList), heteroList)
            hrList = divideList(homoRefList, totalList)
            haList = divideList(homoAltList, totalList)
            hList = divideList(heteroList, totalList)
        elif gene == 'brca1':
            pList = list(df_1['p'])
            qList = list(df_1['q'])
            homoRefList = homoRefList_1
            homoAltList = homoAltList_1
            heteroList = heteroList_1
            totalList = sumList(sumList(homoRefList, homoAltList), heteroList)
            hrList = divideList(homoRefList, totalList)
            haList = divideList(homoAltList, totalList)
            hList = divideList(heteroList, totalList)
        elif gene == 'brca2':
            pList = list(df_2['p'])
            qList = list(df_2['q'])
            homoRefList = homoRefList_2
            homoAltList = homoAltList_2
            heteroList = heteroList_2
            totalList = sumList(sumList(homoRefList, homoAltList), heteroList)
            hrList = divideList(homoRefList, totalList)
            haList = divideList(homoAltList, totalList)
            hList = divideList(heteroList, totalList)
        else:
            print('brca1-report.tsv brca2-report.tsv output-dir all|homo brca1|brca2|both')
            sys.exit(1)
    else:
        print('brca1-report.tsv brca2-report.tsv output-dir all|homo brca1|brca2|both')
        sys.exit(1)

    pSq = list()
    qSq = list()
    twoPQ = list()
    x = list()
    for i in range(len(pList)):
        pSq.append(pList[i]**2)
        qSq.append(qList[i]**2)
        twoPQ.append(2. * pList[i] * qList[i])
        x.append(pList[i])

    plt.scatter(x=x, y=pSq, c='K', s=1)
    plt.scatter(x=x, y=qSq, c='K', s=1)
    plt.scatter(x=x, y=twoPQ, c='K', s=1)
    plt.scatter(x=x, y=haList, c='R', s=1)
    plt.scatter(x=x, y=hrList, c='B', s=1)
    plt.scatter(x=x, y=hList, c='G', s=1)

    plt.title('gene: ' + gene + ' subset: ' + subset)
    plt.savefig(outputDir + '/f8-observed-hw-dists_' + outputFileName + '.png')
    print('gene = ' + gene + ' subset = ' + subset + ' n= ' + str(len(pSq)))

    plt.show()


if __name__ == "__main__":
    main()