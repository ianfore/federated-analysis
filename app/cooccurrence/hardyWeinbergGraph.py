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
        if gene == 'both':
            pList = list(df_1[df_1['homozygousSample'] != 'None']['p']) + list(df_2[df_2['homozygousSample'] != 'None']['p'])
            qList = list(df_1[df_1['homozygousSample'] != 'None']['q']) + list(df_2[df_2['homozygousSample'] != 'None']['q'])
        elif gene == 'brca1':
            pList = list(df_1[df_1['homozygousSample'] != 'None']['p'])
            qList = list(df_1[df_1['homozygousSample'] != 'None']['q'])
        elif gene == 'brca2':
            pList = list(df_2[df_2['homozygousSample'] != 'None']['p'])
            qList = list(df_2[df_2['homozygousSample'] != 'None']['q'])
        else:
            print('brca1-report.tsv brca2-report.tsv output-dir all|homo brca1|brca2|both')
            sys.exit(1)
    elif subset == 'all':
        if gene == 'both':
            pList = list(df_1['p']) + list(df_2['p'])
            qList = list(df_1['q']) + list(df_2['q'])
        elif gene == 'brca1':
            pList = list(df_1['p'])
            qList = list(df_1['q'])
        elif gene == 'brca2':
            pList = list(df_2['p'])
            qList = list(df_2['q'])
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

    plt.scatter(x=x, y=pSq, c='B', s=1)
    plt.scatter(x=x, y=qSq, c='R', s=1)
    plt.scatter(x=x, y=twoPQ, c='G', s=1)
    plt.title('gene: ' + gene + ' subset: ' + subset)
    plt.savefig(outputDir + '/f8-observed-hw-dists_' + outputFileName + '.png')
    print('gene = ' + gene + ' subset = ' + subset + ' n= ' + str(len(pSq)))

    #plt.show()


if __name__ == "__main__":
    main()