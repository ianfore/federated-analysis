import pandas
import sys
import statistics
import pandasql as psql
from tabulate import tabulate

def removeNA(myList):
    return [elt for elt in myList if str(elt) != 'NA']


def run(myFDA):

    results = dict()

    try:

        dfWithPath = myFDA.dataFile[myFDA.dataFile.CarrierGene != 'NonCarrier']
        dfWithoutPath = myFDA.dataFile[myFDA.dataFile.CarrierGene == 'NonCarrier']


        # No. of subjects
        results['No. of subjects'] = list()
        results['No. of subjects'].append(len(dfWithPath))
        results['No. of subjects'].append(len(dfWithoutPath))


        # Age at onset
        results['Age at onset'] = list()
        results['Age at onset'].append(statistics.mean(removeNA(dfWithPath['Age at onset'].tolist())))
        results['Age at onset'].append(statistics.mean(removeNA(dfWithoutPath['Age at onset'].tolist())))

        # Age at entry (?)

        # Age at diagnosis (?)

        # History of ovarian cancer
        results['History of ovarian cancer Yes'] = list()
        results['History of ovarian cancer No'] = list()

        results['History of ovarian cancer Yes'].append(
            psql.sqldf("select * from dfWithPath where `Ovarian cancer history` = 1",locals()).shape[0])
        results['History of ovarian cancer Yes'].append(
            psql.sqldf("select * from dfWithoutPath where `Ovarian cancer history` = 1",locals()).shape[0])
        results['History of ovarian cancer No'].append(
            psql.sqldf("select * from dfWithPath where `Ovarian cancer history` = 0",locals()).shape[0])
        results['History of ovarian cancer No'].append(
            psql.sqldf("select * from dfWithoutPath where `Ovarian cancer history` = 0",locals()).shape[0])


        # Location of cancer (?)

        # TNM clinical classification: N
        results['TNM clinical classification: N 0'] = list()
        results['TNM clinical classification: N 1'] = list()
        results['TNM clinical classification: N 2'] = list()
        results['TNM clinical classification: N 3'] = list()


        results['TNM clinical classification: N 0'].append(
            psql.sqldf("select * from dfWithPath where `TNM classification / N` = 0", locals()).shape[0])
        results['TNM clinical classification: N 0'].append(
            psql.sqldf("select * from dfWithoutPath where `TNM classification / N` = 0", locals()).shape[0])

        results['TNM clinical classification: N 1'].append(
            psql.sqldf("select * from dfWithPath where `TNM classification / N` = 1", locals()).shape[0])
        results['TNM clinical classification: N 1'].append(
            psql.sqldf("select * from dfWithoutPath where `TNM classification / N` = 1", locals()).shape[0])

        results['TNM clinical classification: N 2'].append(
            psql.sqldf("select * from dfWithPath where `TNM classification / N` = 2", locals()).shape[0])
        results['TNM clinical classification: N 2'].append(
            psql.sqldf("select * from dfWithoutPath where `TNM classification / N` = 2", locals()).shape[0])

        results['TNM clinical classification: N 3'].append(
            psql.sqldf("select * from dfWithPath where `TNM classification / N` = 3", locals()).shape[0])
        results['TNM clinical classification: N 3'].append(
            psql.sqldf("select * from dfWithoutPath where `TNM classification / N` = 3", locals()).shape[0])



        # TNM clinical classification: M
        results['TNM clinical classification: M 0'] = list()
        results['TNM clinical classification: M 1'] = list()

        results['TNM clinical classification: M 0'].append(
            psql.sqldf("select * from dfWithPath where `TNM classification / M` = 0", locals()).shape[0])
        results['TNM clinical classification: M 0'].append(
            psql.sqldf("select * from dfWithoutPath where `TNM classification / M` = 0", locals()).shape[0])
        results['TNM clinical classification: M 1'].append(
            psql.sqldf("select * from dfWithPath where `TNM classification / M` = 1", locals()).shape[0])
        results['TNM clinical classification: M 1'].append(
            psql.sqldf("select * from dfWithoutPath where `TNM classification / M` = 1", locals()).shape[0])



        # Estrogen-receptor status
        results['Estrogen-receptor status Positive'] = list()
        results['Estrogen-receptor status Negative'] = list()

        results['Estrogen-receptor status Positive'].append(
            psql.sqldf("select * from dfWithPath where `ER` = 'Positive'", locals()).shape[0])
        results['Estrogen-receptor status Positive'].append(
            psql.sqldf("select * from dfWithoutPath where `ER` = 'Positive'", locals()).shape[0])
        results['Estrogen-receptor status Negative'].append(
            psql.sqldf("select * from dfWithPath where `ER` = 'Negative'", locals()).shape[0])
        results['Estrogen-receptor status Negative'].append(
            psql.sqldf("select * from dfWithoutPath where `ER` = 'Negative'", locals()).shape[0])


        # Progesterone-receptor status
        results['Progesterone-receptor status Positive'] = list()
        results['Progesterone-receptor status Negative'] = list()

        results['Progesterone-receptor status Positive'].append(
            psql.sqldf("select * from dfWithPath where `PgR` = 'Positive'",locals()).shape[0])
        results['Progesterone-receptor status Positive'].append(
            psql.sqldf("select * from dfWithoutPath where `PgR` = 'Positive'", locals()).shape[0])
        results['Progesterone-receptor status Negative'].append(
            psql.sqldf("select * from dfWithPath where `PgR` = 'Negative'", locals()).shape[0])
        results['Progesterone-receptor status Negative'].append(
            psql.sqldf("select * from dfWithoutPath where `PgR` = 'Negative'", locals()).shape[0])



        # Triple negative breast cancer
        results['Triple negative breast cancer Yes'] = list()
        results['Triple negative breast cancer No'] = list()

        results['Triple negative breast cancer Yes'].append(
            psql.sqldf("select * from dfWithPath where `PgR` = 'Negative' and `ER` = 'Negative' and `HER2` = '0'",
                       locals()).shape[0])
        results['Triple negative breast cancer Yes'].append(
            psql.sqldf("select * from dfWithoutPath where `PgR` = 'Negative' and `ER` = 'Negative' and `HER2` = '0'",
                       locals()).shape[0])
        results['Triple negative breast cancer No'].append(
            psql.sqldf("select * from dfWithPath where not `PgR` = 'Negative' or not `ER` = 'Negative' or not `HER2` = '0'",
                       locals()).shape[0])
        results['Triple negative breast cancer No'].append(
            psql.sqldf("select * from dfWithoutPath where not `PgR` = 'Negative' or not `ER` = 'Negative' or not `HER2` = '0'",
                       locals()).shape[0])

        # Family history of breast cancer
        results['Family history of breast cancer Yes'] = list()
        results['Family history of breast cancer No'] = list()

        results['Family history of breast cancer Yes'].append(
            psql.sqldf("select * from dfWithPath where `Family history / breast cancer` = 1",locals()).shape[0])
        results['Family history of breast cancer Yes'].append(
            psql.sqldf("select * from dfWithoutPath where `Family history / breast cancer` = 1", locals()).shape[0])
        results['Family history of breast cancer No'].append(
            psql.sqldf("select * from dfWithPath where `Family history / breast cancer` = 0", locals()).shape[0])
        results['Family history of breast cancer No'].append(
            psql.sqldf("select * from dfWithoutPath where `Family history / breast cancer` = 0",locals()).shape[0])

        # Family history of ovarian cancer
        results['Family history of ovarian cancer Yes'] = list()
        results['Family history of ovarian cancer No'] = list()

        results['Family history of ovarian cancer Yes'].append(
            psql.sqldf("select * from dfWithPath where `Family history / ovarian cancer` = 1", locals()).shape[0])
        results['Family history of ovarian cancer Yes'].append(
            psql.sqldf("select * from dfWithoutPath where `Family history / ovarian cancer` = 1", locals()).shape[0])
        results['Family history of ovarian cancer No'].append(
            psql.sqldf("select * from dfWithPath where `Family history / ovarian cancer` = 0", locals()).shape[0])
        results['Family history of ovarian cancer No'].append(
            psql.sqldf("select * from dfWithoutPath where `Family history / ovarian cancer` = 0", locals()).shape[0])

        # Family history of pancreas cancer
        results['Family history of pancreas cancer Yes'] = list()
        results['Family history of pancreas cancer No'] = list()

        results['Family history of pancreas cancer Yes'].append(
            psql.sqldf("select * from dfWithPath where `Family history / pancreatic cancer` = 1", locals()).shape[0])
        results['Family history of pancreas cancer Yes'].append(
            psql.sqldf("select * from dfWithoutPath where `Family history / pancreatic cancer` = 1", locals()).shape[0])
        results['Family history of pancreas cancer No'].append(
            psql.sqldf("select * from dfWithPath where `Family history / pancreatic cancer` = 0", locals()).shape[0])
        results['Family history of pancreas cancer No'].append(
            psql.sqldf("select * from dfWithoutPath where `Family history / pancreatic cancer` = 0", locals()).shape[0])

        # Family history of gastric (stomach?) cancer
        results['Family history of stomach cancer Yes'] = list()
        results['Family history of stomach cancer No'] = list()

        results['Family history of stomach cancer Yes'].append(
            psql.sqldf("select * from dfWithPath where `Family history / stomach cancer` = 1", locals()).shape[0])
        results['Family history of stomach cancer Yes'].append(
            psql.sqldf("select * from dfWithoutPath where `Family history / stomach cancer` = 1", locals()).shape[0])
        results['Family history of stomach cancer No'].append(
            psql.sqldf("select * from dfWithPath where `Family history / stomach cancer` = 0", locals()).shape[0])
        results['Family history of stomach cancer No'].append(
            psql.sqldf("select * from dfWithoutPath where `Family history / stomach cancer` = 0", locals()).shape[0])

        # Family history of liver cancer
        results['Family history of liver cancer Yes'] = list()
        results['Family history of liver cancer No'] = list()

        results['Family history of liver cancer Yes'].append(
            psql.sqldf("select * from dfWithPath where `Family history / liver cancer` = 1", locals()).shape[0])
        results['Family history of liver cancer Yes'].append(
            psql.sqldf("select * from dfWithoutPath where `Family history / liver cancer` = 1", locals()).shape[0])
        results['Family history of liver cancer No'].append(
            psql.sqldf("select * from dfWithPath where `Family history / liver cancer` = 0", locals()).shape[0])
        results['Family history of liver cancer No'].append(
            psql.sqldf("select * from dfWithoutPath where `Family history / liver cancer` = 0", locals()).shape[0])

        # Family history of bone tumor
        results['Family history of bone tumor Yes'] = list()
        results['Family history of bone tumor No'] = list()

        results['Family history of bone tumor Yes'].append(
            psql.sqldf("select * from dfWithPath where `Family history / bone tumor` = 1", locals()).shape[0])
        results['Family history of bone tumor Yes'].append(
            psql.sqldf("select * from dfWithoutPath where `Family history / bone tumor` = 1", locals()).shape[0])
        results['Family history of bone tumor No'].append(
            psql.sqldf("select * from dfWithPath where `Family history / bone tumor` = 0", locals()).shape[0])
        results['Family history of bone tumor No'].append(
            psql.sqldf("select * from dfWithoutPath where `Family history / bone tumor` = 0", locals()).shape[0])

        # Family history of bladder cancer
        results['Family history of bladder cancer Yes'] = list()
        results['Family history of bladder cancer No'] = list()

        results['Family history of bladder cancer Yes'].append(
            psql.sqldf("select * from dfWithPath where `Family history / bladder cancer` = 1", locals()).shape[0])
        results['Family history of bladder cancer Yes'].append(
            psql.sqldf("select * from dfWithoutPath where `Family history / bladder cancer` = 1", locals()).shape[0])
        results['Family history of bladder cancer No'].append(
            psql.sqldf("select * from dfWithPath where `Family history / bladder cancer` = 0", locals()).shape[0])
        results['Family history of bladder cancer No'].append(
            psql.sqldf("select * from dfWithoutPath where `Family history / bladder cancer` = 0", locals()).shape[0])


        # define fileObject based on config
        if myFDA.configFile.outputFile == "":
            fileObject = sys.stdout
        else:
            fileObject = open(myFDA.configFile.outputFile, mode='a')

        # print results to fileObject
        '''for key in results.keys():
            print(key, results[key], sep='\t', file=fileObject)'''
        prettyPrint(results)



        return True

    except Exception as e:
        print('exception in supplementaryTable4.run() method: ' + str(e))
        return False

def prettyPrint(results):
    # convert dict to data frame
    df = pandas.DataFrame(results)

    # pretty print the data frame
    print(tabulate(df.T, headers='keys', tablefmt='psql'))
