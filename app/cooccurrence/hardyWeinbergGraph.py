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
    AA = list(df_1['homo_ref'])
    Aa = list(df_1['hetero'])
    aa = list(df_1['homo_alt'])
    pList = list(df_1['p'])
    qList = list(df_1['q'])

    pSq = list()
    pSqFreq = list()
    qSq = list()
    qSqFreq = list()
    twoPQ = list()
    twoPQFreq = list()
    n = len(df_1)
    for i in range(n):
        pSq.append(pList[i]**2)
        pSqFreq.append(1. * AA[i] / 1. * (aa[i] + Aa[i] + AA[i]))
        qSq.append((qList[i] - 1)**2)
        qSqFreq.append(1. * aa[i] / 1. * (aa[i] + Aa[i] + AA[i]))
        twoPQ.append(2 * pList[i] * qList[i])
        twoPQFreq.append(1. * Aa[i] / 1. * (aa[i] + Aa[i] + AA[i]))

    plt.scatter(x=pSq, y=pSqFreq)
    plt.scatter(x=qSq, y=qSqFreq)
    plt.scatter(x=twoPQ, y=twoPQFreq)
    plt.show()




if __name__ == "__main__":
    main()