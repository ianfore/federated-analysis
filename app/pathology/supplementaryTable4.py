import pandas
import sys
import statistics
import pandasql as psql
from tabulate import tabulate
import scipy.stats as stats
from math import exp, log, sqrt
import subprocess

def removeNA(myList):
    return [elt for elt in myList if str(elt) != 'NA']


def run(myFDA):

    results = dict()
    RScriptPath = myFDA.configFile.RScriptPath
    fieldsOfInterest = dict()
    for fieldDict in myFDA.configFile.fieldFilters:
        fieldsOfInterest[fieldDict['fieldName']] = fieldDict['fieldValues']

    try:
        dfWithPath = myFDA.dataFile[myFDA.dataFile.CarrierGene != 'NonCarrier']
        dfWithoutPath = myFDA.dataFile[myFDA.dataFile.CarrierGene == 'NonCarrier']

        # No. of subjects
        results['No. of subjects'] = {
            'with': len(dfWithPath),
            'without': len(dfWithoutPath)}

        # Age at onset
        try:
            results['Age at onset'] = {
                'with': round(statistics.mean(removeNA(dfWithPath['Age at onset'].tolist())), 2),
                'without': round(statistics.mean(removeNA(dfWithoutPath['Age at onset'].tolist())), 2)}
        except Exception as e:
            pass
            #return False


        # TNM clinical classification: N
        if 'TNM classification / N' in fieldsOfInterest:
            try:
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

                getPercentage(results, 'TNM clinical classification N', '0', ['0', '1', '2', '3'])
                getPercentage(results, 'TNM clinical classification N', '1', ['0', '1', '2', '3'])
                getPercentage(results, 'TNM clinical classification N', '2', ['0', '1', '2', '3'])
                getPercentage(results, 'TNM clinical classification N', '3', ['0', '1', '2', '3'])
            except Exception as e:
                pass
                #return False


        # Estrogen-receptor status
        if 'ER' in fieldsOfInterest:
            try:
                results['Estrogen-receptor status'] = {
                    'with': {
                        'Positive': psql.sqldf("select * from dfWithPath where `ER` = 'Positive'", locals()).shape[0],
                        'Negative': psql.sqldf("select * from dfWithPath where `ER` = 'Negative'", locals()).shape[0]},
                    'without': {
                        'Positive': psql.sqldf("select * from dfWithoutPath where `ER` = 'Positive'", locals()).shape[0],
                        'Negative': psql.sqldf("select * from dfWithoutPath where `ER` = 'Negative'", locals()).shape[0]}}

                getPercentage(results, 'Estrogen-receptor status', 'Positive', ['Positive', 'Negative'])
                getFisherExact(RScriptPath, results, 'Estrogen-receptor status', ['Positive', 'Negative'])
            except Exception as e:
                pass
                #return False

        # Progesterone-receptor status
        if 'PgR' in fieldsOfInterest:
            try:
                results['Progesterone-receptor status'] = {
                    'with': {
                        'Positive': psql.sqldf("select * from dfWithPath where `PgR` = 'Positive'", locals()).shape[0],
                        'Negative': psql.sqldf("select * from dfWithPath where `PgR` = 'Negative'", locals()).shape[0]},
                    'without': {
                        'Positive': psql.sqldf("select * from dfWithoutPath where `PgR` = 'Positive'", locals()).shape[0],
                        'Negative': psql.sqldf("select * from dfWithoutPath where `PgR` = 'Negative'", locals()).shape[0]}}

                getPercentage(results, 'Progesterone-receptor status', 'Positive', ['Positive', 'Negative'])
                getFisherExact(RScriptPath, results, 'Progesterone-receptor status', ['Positive', 'Negative'])
            except Exception as e:
                pass
                #return False

        # Triple negative breast cancer
        if 'PgR' in fieldsOfInterest and 'ER' in fieldsOfInterest and 'HER2' in fieldsOfInterest:
            try:
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

                getPercentage(results, 'Triple negative breast cancer', 'Yes', ['Yes', 'No'])
                getFisherExact(RScriptPath, results, 'Triple negative breast cancer', ['Yes', 'No'])
            except Exception as e:
                pass
            #return False

        # TNM clinical classification: M
        callSQLOnBinary(results, dfWithPath, dfWithoutPath, 'TNM classification / M', fieldsOfInterest, RScriptPath)

        # Family history of breast cancer
        callSQLOnBinary(results, dfWithPath, dfWithoutPath, 'Family history / breast cancer', fieldsOfInterest, RScriptPath)

        # Family history of ovarian cancer
        callSQLOnBinary(results, dfWithPath, dfWithoutPath, 'Family history / ovarian cancer', fieldsOfInterest, RScriptPath)

        # Family history of pancreas cancer
        callSQLOnBinary(results, dfWithPath, dfWithoutPath, 'Family history / pancreatic cancer', fieldsOfInterest, RScriptPath)

        # Family history of gastric (stomach?) cancer
        callSQLOnBinary(results, dfWithPath, dfWithoutPath, 'Family history / stomach cancer', fieldsOfInterest, RScriptPath)

        # Family history of liver cancer
        callSQLOnBinary(results, dfWithPath, dfWithoutPath, 'Family history / liver cancer', fieldsOfInterest, RScriptPath)

        # Family history of bone tumor
        callSQLOnBinary(results, dfWithPath, dfWithoutPath, 'Family history / bone tumor', fieldsOfInterest, RScriptPath)

        # Family history of bladder cancer
        callSQLOnBinary(results, dfWithPath, dfWithoutPath, 'Family history / bladder cancer', fieldsOfInterest, RScriptPath)


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

def callSQLOnBinary(results, dfWithPath, dfWithoutPath, fieldName, fieldsOfInterest, RScriptPath):
    if not fieldName in fieldsOfInterest:
        return
    try:
        results[fieldName] = {
            'with': {
                'Yes':
                    psql.sqldf("select * from dfWithPath  where `" + fieldName + "` = 1", locals()).shape[
                        0],
                'No':
                    psql.sqldf("select * from dfWithPath where `" + fieldName + "` = 0", locals()).shape[
                        0]},
            'without': {
                'Yes': psql.sqldf("select * from dfWithoutPath where `" + fieldName + "` = 1",
                                  locals()).shape[0],
                'No': psql.sqldf("select * from dfWithoutPath where `" + fieldName + "` = 0",
                                 locals()).shape[0]}}
    except Exception as e:
        pass
    getPercentage(results, fieldName, 'Yes', ['Yes', 'No'])
    getFisherExact(RScriptPath, results, fieldName, ['Yes', 'No'])


def getPercentage(results, key, value, allValues):
    for path in ['with', 'without']:
        denom = 0
        num = results[key][path][value]
        for v in allValues:
            denom += results[key][path][v]
        results[key][path]['% ' + value] = str(round(100 * num / denom, 2)) + '%'

def getFisherExact(RScriptPath, results, key, allValues):
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
    CIcmd_1 = str(RScriptPath + ' -e "fisher.test(rbind(' + row1 + ',' +  row2 + '))"')
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


def prettyPrint(results, fileObject):
    # convert dict to data frame and reorder rows (which will become columns after transposing it)
    df = pandas.DataFrame(results).reindex(['with', 'without', 'P value', 'OR', '95% CI']).T

    # pretty print the data frame
    print(tabulate(df, headers='keys', tablefmt='psql'), file=fileObject)
