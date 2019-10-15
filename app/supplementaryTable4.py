import pandas
import sys
import statistics
import pandasql as psql
from tabulate import tabulate
import scipy.stats as stats
from math import exp, log, sqrt

def removeNA(myList):
    return [elt for elt in myList if str(elt) != 'NA']


def run(myFDA):

    results = dict()

    try:
        dfWithPath = myFDA.dataFile[myFDA.dataFile.CarrierGene != 'NonCarrier']
        dfWithoutPath = myFDA.dataFile[myFDA.dataFile.CarrierGene == 'NonCarrier']

        # No. of subjects
        results['No. of subjects'] = {
            'with': len(dfWithPath),
            'without': len(dfWithoutPath)}

        # Age at onset
        results['Age at onset'] = {
            'with': round(statistics.mean(removeNA(dfWithPath['Age at onset'].tolist())), 2),
            'without': round(statistics.mean(removeNA(dfWithoutPath['Age at onset'].tolist())), 2)}

        # Age at entry (?)

        # Age at diagnosis (?)

        # History of ovarian cancer
        results['History of ovarian cancer'] = {
            'with': {
                'Yes': psql.sqldf("select * from dfWithPath where `Ovarian cancer history` = 1",locals()).shape[0],
                'No': psql.sqldf("select * from dfWithPath where `Ovarian cancer history` = 0",locals()).shape[0]},
            'without': {
                'Yes': psql.sqldf("select * from dfWithoutPath where `Ovarian cancer history` = 1",locals()).shape[0],
                'No': psql.sqldf("select * from dfWithoutPath where `Ovarian cancer history` = 0",locals()).shape[0]}}

        getPctg(results, 'History of ovarian cancer', 'Yes', ['Yes', 'No'])
        getFisherExact(results, 'History of ovarian cancer', ['Yes', 'No'])

        # Location of cancer (?)

        # TNM clinical classification: N
        results['TNM clinical classification N'] = {
            'with': {
                '0': psql.sqldf("select * from dfWithPath where `TNM classification / N` = 0", locals()).shape[0],
                '1': psql.sqldf("select * from dfWithPath where `TNM classification / N` = 1", locals()).shape[0],
                '2': psql.sqldf("select * from dfWithPath where `TNM classification / N` = 2", locals()).shape[0],
                '3': psql.sqldf("select * from dfWithPath where `TNM classification / N` = 3", locals()).shape[0]},
            'without': {
                '0': psql.sqldf("select * from dfWithoutPath where `TNM classification / N` = 0", locals()).shape[0],
                '1': psql.sqldf("select * from dfWithoutPath where `TNM classification / N` = 1",locals()).shape[0],
                '2': psql.sqldf("select * from dfWithoutPath where `TNM classification / N` = 2",locals()).shape[0],
                '3': psql.sqldf("select * from dfWithoutPath where `TNM classification / N` = 3",locals()).shape[0]}}

        getPctg(results, 'TNM clinical classification N', '0', ['0', '1', '2', '3'])
        getPctg(results, 'TNM clinical classification N', '1', ['0', '1', '2', '3'])
        getPctg(results, 'TNM clinical classification N', '2', ['0', '1', '2', '3'])
        getPctg(results, 'TNM clinical classification N', '3', ['0', '1', '2', '3'])



        # TNM clinical classification: M
        results['TNM clinical classification M'] = {
            'with': {
                '0': psql.sqldf("select * from dfWithPath where `TNM classification / M` = 0", locals()).shape[0],
                '1': psql.sqldf("select * from dfWithPath where `TNM classification / M` = 1", locals()).shape[0]},
            'without': {
                '0': psql.sqldf("select * from dfWithoutPath where `TNM classification / M` = 0", locals()).shape[0],
                '1': psql.sqldf("select * from dfWithoutPath where `TNM classification / M` = 1",locals()).shape[0]}}

        getPctg(results, 'TNM clinical classification M', '0', ['0', '1'])
        getPctg(results, 'TNM clinical classification M', '1', ['0', '1'])

        # Estrogen-receptor status
        results['Estrogen-receptor status'] = {
            'with': {
                'Positive': psql.sqldf("select * from dfWithPath where `ER` = 'Positive' and `HER2` != 'NA' and \
            `PgR` != 'NA'", locals()).shape[0],
                'Negative': psql.sqldf("select * from dfWithPath where `ER` = 'Negative' and `HER2` != 'NA' and \
            `PgR` != 'NA'", locals()).shape[0]},
            'without': {
                'Positive': psql.sqldf("select * from dfWithoutPath where `ER` = 'Positive' and `HER2` != 'NA' and \
            `PgR` != 'NA'", locals()).shape[0],
                'Negative': psql.sqldf("select * from dfWithoutPath where `ER` = 'Negative' and `HER2` != 'NA' and \
            `PgR` != 'NA'", locals()).shape[0]}}

        getPctg(results, 'Estrogen-receptor status', 'Positive', ['Positive', 'Negative'])
        getFisherExact(results, 'Estrogen-receptor status', ['Positive', 'Negative'])


        # Progesterone-receptor status
        results['Progesterone-receptor status'] = {
            'with': {
                'Positive': psql.sqldf("select * from dfWithPath where `ER` != 'NA' and `HER2` != 'NA' and \
                `PgR` = 'Positive'", locals()).shape[0],
                'Negative': psql.sqldf("select * from dfWithPath where `ER` != 'NA' and `HER2` != 'NA' and \
                `PgR` = 'Negative'", locals()).shape[0]},
            'without': {
                'Positive': psql.sqldf("select * from dfWithoutPath where `ER` != 'NA' and `HER2` != 'NA' and \
                `PgR` = 'Positive'", locals()).shape[0],
                'Negative': psql.sqldf("select * from dfWithoutPath where `ER` != 'NA' and `HER2` != 'NA' and \
                `PgR` = 'Negative'", locals()).shape[0]}}

        getPctg(results, 'Progesterone-receptor status', 'Positive', ['Positive', 'Negative'])
        getFisherExact(results, 'Progesterone-receptor status', ['Positive', 'Negative'])


        # Triple negative breast cancer
        results['Triple negative breast cancer'] = {
            'with': {
                'Yes': psql.sqldf("select * from dfWithPath where (`PgR` = 'Negative' and `ER` = 'Negative') and \
                (`HER2` = '0'  or `HER2` = '1+')", locals()).shape[0],
                'No': psql.sqldf("select * from (select * from dfWithPath where  `PgR` != 'NA' and  `ER` != 'NA' \
                and  `HER2` != 'NA') T where  T.`PgR` != 'Negative' or  T.`ER` != 'Negative' or  (T.`HER2` != '0' and \
                T.`HER2` != '1+')", locals()).shape[0]},
            'without': {
                'Yes': psql.sqldf("select * from dfWithoutPath where (`PgR` = 'Negative' and `ER` = 'Negative') and \
                (`HER2` = '0'  or `HER2` = '1+')", locals()).shape[0],
                'No': psql.sqldf("select * from (select * from dfWithoutPath where  `PgR` != 'NA' and  `ER` != 'NA' \
                and  `HER2` != 'NA') T where  T.`PgR` != 'Negative' or  T.`ER` != 'Negative' or  (T.`HER2` != '0' and \
                T.`HER2` != '1+')", locals()).shape[0]}}

        getPctg(results, 'Triple negative breast cancer', 'Yes', ['Yes', 'No'])
        getFisherExact(results, 'Triple negative breast cancer', ['Yes', 'No'])


        # Family history of breast cancer
        results['Family history of breast cancer'] = {
            'with': {
                'Yes': psql.sqldf("select * from dfWithPath where `Family history / breast cancer` = 1",locals()).shape[0],
                'No': psql.sqldf("select * from dfWithPath where `Family history / breast cancer` = 0", locals()).shape[0]},
            'without': {
                'Yes': psql.sqldf("select * from dfWithoutPath where `Family history / breast cancer` = 1", locals()).shape[0],
                'No': psql.sqldf("select * from dfWithoutPath where `Family history / breast cancer` = 0",locals()).shape[0]}}

        getPctg(results, 'Family history of breast cancer', 'Yes', ['Yes', 'No'])
        getFisherExact(results, 'Family history of breast cancer', ['Yes', 'No'])


        # Family history of ovarian cancer
        results['Family history of ovarian cancer'] = {
            'with': {
                'Yes': psql.sqldf("select * from dfWithPath where `Family history / ovarian cancer` = 1", locals()).shape[0],
                'No': psql.sqldf("select * from dfWithPath where `Family history / ovarian cancer` = 0", locals()).shape[0]},
            'without': {
                'Yes': psql.sqldf("select * from dfWithoutPath where `Family history / ovarian cancer` = 1", locals()).shape[0],
                'No': psql.sqldf("select * from dfWithoutPath where `Family history / ovarian cancer` = 0", locals()).shape[0]}}

        getPctg(results, 'Family history of ovarian cancer', 'Yes', ['Yes', 'No'])
        getFisherExact(results, 'Family history of ovarian cancer', ['Yes', 'No'])


        # Family history of pancreas cancer
        results['Family history of pancreas cancer'] = {
            'with': {
                'Yes': psql.sqldf("select * from dfWithPath where `Family history / pancreatic cancer` = 1", locals()).shape[0],
                'No': psql.sqldf("select * from dfWithPath where `Family history / pancreatic cancer` = 0", locals()).shape[0]},
            'without': {
                'Yes': psql.sqldf("select * from dfWithoutPath where `Family history / pancreatic cancer` = 1", locals()).shape[0],
                'No': psql.sqldf("select * from dfWithoutPath where `Family history / pancreatic cancer` = 0", locals()).shape[0]}}

        getPctg(results, 'Family history of pancreas cancer', 'Yes', ['Yes', 'No'])
        getFisherExact(results, 'Family history of pancreas cancer', ['Yes', 'No'])


        # Family history of gastric (stomach?) cancer
        results['Family history of stomach cancer'] = {
            'with': {
                'Yes': psql.sqldf("select * from dfWithPath where `Family history / stomach cancer` = 1", locals()).shape[0],
                'No': psql.sqldf("select * from dfWithPath where `Family history / stomach cancer` = 0", locals()).shape[0]},
            'without': {
                'Yes': psql.sqldf("select * from dfWithoutPath where `Family history / stomach cancer` = 1", locals()).shape[0],
                'No': psql.sqldf("select * from dfWithoutPath where `Family history / stomach cancer` = 0", locals()).shape[0]}}

        getPctg(results, 'Family history of stomach cancer', 'Yes', ['Yes', 'No'])
        getFisherExact(results, 'Family history of stomach cancer', ['Yes', 'No'])


        # Family history of liver cancer
        results['Family history of liver cancer'] = {
            'with': {
                'Yes': psql.sqldf("select * from dfWithPath where `Family history / liver cancer` = 1", locals()).shape[0],
                'No': psql.sqldf("select * from dfWithPath where `Family history / liver cancer` = 0", locals()).shape[0]},
            'without': {
                'Yes': psql.sqldf("select * from dfWithoutPath where `Family history / liver cancer` = 1", locals()).shape[0],
                'No': psql.sqldf("select * from dfWithoutPath where `Family history / liver cancer` = 0", locals()).shape[0]}}

        getPctg(results, 'Family history of liver cancer', 'Yes', ['Yes', 'No'])
        getFisherExact(results, 'Family history of liver cancer', ['Yes', 'No'])


        # Family history of bone tumor
        results['Family history of bone tumor'] = {
            'with': {
                'Yes': psql.sqldf("select * from dfWithPath where `Family history / bone tumor` = 1", locals()).shape[0],
                'No': psql.sqldf("select * from dfWithPath where `Family history / bone tumor` = 0", locals()).shape[0]},
            'without': {
                'Yes': psql.sqldf("select * from dfWithoutPath where `Family history / bone tumor` = 1", locals()).shape[0],
                'No': psql.sqldf("select * from dfWithoutPath where `Family history / bone tumor` = 0", locals()).shape[0]}}

        getPctg(results, 'Family history of bone tumor', 'Yes', ['Yes', 'No'])
        getFisherExact(results, 'Family history of bone tumor', ['Yes', 'No'])


        # Family history of bladder cancer
        results['Family history of bladder cancer'] = {
            'with': {
                'Yes': psql.sqldf("select * from dfWithPath where `Family history / bladder cancer` = 1", locals()).shape[0],
                'No': psql.sqldf("select * from dfWithPath where `Family history / bladder cancer` = 0", locals()).shape[0]},
            'without': {
                'Yes': psql.sqldf("select * from dfWithoutPath where `Family history / bladder cancer` = 1", locals()).shape[0],
                'No': psql.sqldf("select * from dfWithoutPath where `Family history / bladder cancer` = 0", locals()).shape[0]}}

        getPctg(results, 'Family history of bladder cancer', 'Yes', ['Yes', 'No'])
        getFisherExact(results, 'Family history of bladder cancer', ['Yes', 'No'])


        '''results['his history of ovarian cancer'] = {
            'with': {
                'Yes': 7,
                'No': 397},
            'without': {
                'Yes': 40,
                'No': 6607}}
        getPctg(results, 'his history of ovarian cancer', 'Yes', ['Yes', 'No'])
        getFisherExact(results, 'his history of ovarian cancer', ['Yes', 'No'])'''

        # define fileObject based on config
        if myFDA.configFile.outputFile == "":
            fileObject = sys.stdout
        else:
            fileObject = open(myFDA.configFile.outputFile, mode='a')

        prettyPrint(results, fileObject)

        return True

    except Exception as e:
        print('exception in supplementaryTable4.run() method: ' + str(e))
        return False


def getPctg(results, key, value, allValues):
    for path in ['with', 'without']:
        denom = 0
        num = results[key][path][value]
        for v in allValues:
            denom += results[key][path][v]
        results[key][path]['% ' + value] = str(round(100 * num / denom, 2)) + '%'

def getFisherExact(results, key, allValues):
    # create 2x2 contingency table
    a = results[key]['with'][allValues[0]]
    b = results[key]['without'][allValues[0]]
    c = results[key]['with'][allValues[1]]
    d = results[key]['without'][allValues[1]]

    # run fisher exact test
    oddsRatio, pValue = stats.fisher_exact([[a, b], [c, d]])

    # get confidence interval for odds ratio: CI = e^(ln(OR) +/- [1.96 * sqrt(1/a + 1/b + 1/c + 1/d)])

    x = log(oddsRatio)
    y = 1.96 * sqrt(1/a + 1/b + 1/c + 1/d)
    lowerBound = exp(x - y)
    upperBound = exp(x + y)

    # insert OR and p-value into results dict
    results[key]['OR'] = round(oddsRatio, 3)
    results[key]['P value'] =  round(pValue, 3)
    results[key]['95% CI'] = (round(lowerBound, 2), round(upperBound, 2))


def prettyPrint(results, fileObject):
    # convert dict to data frame and reorder rows (which will become columns after transposing it)
    df = pandas.DataFrame(results).reindex(['with', 'without', 'P value', 'OR', '95% CI']).T

    # pretty print the data frame
    print(tabulate(df, headers='keys', tablefmt='psql'), file=fileObject)
