import math
import pandas as pd


def main():
    brca1DF = pd.read_csv('/Users/jcasaletto/Desktop/RESEARCH/TOPMED/REPORT/F8/f8_chr17_brca1_gruhmb_report.tsv', header=0,
                          sep='\t')
    brca2DF = pd.read_csv('/Users/jcasaletto/Desktop/RESEARCH/TOPMED/REPORT/F8/f8_chr13_brca2_gruhmb_report.tsv', header=0,
                          sep='\t')

    b1_zeroList = list(brca1DF['hail_hweafp'])
    b1_oneList = list(brca1DF['q'])
    b1_psb = pearsonCorrelation(b1_zeroList, b1_oneList)
    print(b1_psb)

    b2_zeroList = list(brca2DF['hail_hweafp'])
    b2_oneList = list(brca2DF['q'])
    b2_psb = pearsonCorrelation(b2_zeroList, b2_oneList)
    print(b2_psb)

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

def pointBiserialCorrelation(zeroList, oneList):
    # rpb = ((M1 - M0) / sn ) * sqrt(n1 * n0 / n**2)
    M0 = mean(zeroList)
    M1 = mean(oneList)
    sn = standDev(zeroList + oneList)
    n0 = len(zeroList)
    n1 = len(oneList)
    n = n0 + n1
    # return rpb
    return ((M1 - M0) / sn) * math.sqrt(n1 * n0/ n**2)



if __name__ == "__main__":
    main()
