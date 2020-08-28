import math
import pandas as pd
from scipy.stats import pointbiserialr
import sys
from sklearn import preprocessing



def main():
    brca1DF = pd.read_csv('/Users/jcasaletto/Desktop/FINAL_REPORT/F8/f8_chr17_brca1_gruhmb_report.tsv', header=0,
                          sep='\t', dtype={'exonic': 'bool', 'inGnomad': 'bool'})
    brca2DF = pd.read_csv('/Users/jcasaletto/Desktop/FINAL_REPORT/F8/f8_chr13_brca2_gruhmb_report.tsv', header=0,
                          sep='\t', dtype={'exonic': 'bool', 'inGnomad': 'bool'})

    list1 = list(brca1DF['q'])
    list2 = list(brca1DF['HWEAF_P'])

    le = preprocessing.LabelEncoder()

    if isinstance(list1[0], bool) and not isinstance(list2[0], bool):
        corr = pointbiserialr(x=list1, y=list2)
    elif isinstance(list2[0], bool) and not isinstance(list1[0], bool):
        corr = pointbiserialr(x=list2, y=list1)
    elif isinstance(list1[0], bool) and isinstance(list2[0], bool):
        print('both lists cannot be boolean')
        sys.exit(1)
    elif isinstance(list1[0], str) and not isinstance(list2[0], str):
        list1 = le.fit_transform(list1)
        corr = pearsonCorrelation(list1, list2)
    elif isinstance(list2[0], str) and not isinstance(list1[0], str):
        list2 = le.fit_transform(list2)
        corr = pearsonCorrelation(list1, list2)
    else:
        corr = pearsonCorrelation(list1, list2)
    print(corr)

def convertStringToBool(someList):
    retList = list()
    if isinstance(someList[0], str):
        for x in someList:
            if x.lower() == 'true':
                retList.append(True)
            else:
                retList.append(False)
    return retList


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
