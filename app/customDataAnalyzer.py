import sys
import json

class ReceptorCounts:
    """
    Purpose:
        Encapsulate count data for each of the receptors ER, PgR, and HER2.

    data members:
        counts: dictionary of counts for each "(receptor, status)" tuple string
    """

    def __init__(self):
        self.counts = dict()

    def increment(self, myRow):
        for receptor in [ 'ER', 'PgR', 'HER2']:
            # add to overall count
            if myRow[receptor] == 'NA':
                continue
            else:
                if str((receptor, myRow[receptor])) not in self.counts:
                    self.counts[str((receptor, myRow[receptor]))] = 1
                else:
                    self.counts[str((receptor, myRow[receptor]))] += 1
                # add to <50
                if myRow['Age at onset'] != 'NA' and int(myRow['Age at onset']) < 50:
                    if str((receptor, myRow[receptor], "age<50")) not in self.counts:
                        self.counts[str((receptor, myRow[receptor], "age<50"))] = 1
                    else:
                        self.counts[str((receptor, myRow[receptor], "age<50"))] += 1
                # add to >= 50
                else:
                    if str((receptor, myRow[receptor], "age>=50")) not in self.counts:
                        self.counts[str((receptor, myRow[receptor], "age>=50"))] = 1
                    else:
                        self.counts[str((receptor, myRow[receptor], "age>=50"))] += 1


    def print(self):
        print(json.dumps(self.__dict__, sort_keys=True))


def run(myFDA):
    try:
        # initialize counters
        countsPerGene = dict()
        tripleNegatives = list()
        triplePositives = list()
        tripleValues = list()
        nonTripleNegatives = list()

        skipTheseGenes = ['NonCarrier']

        # iterate through data frame to gather counts
        for myIndex, myRow in myFDA.dataFile.iterrows():
            # if all fields != NA, then add to that list
            if (isTripleValue(myRow)):
                tripleValues.append(myIndex)
                # if triple negative, add to that list
                if(isTripleNegative(myRow)):
                    tripleNegatives.append(myIndex)
                # if triple positive, add to that list
                elif(isTriplePositive(myRow)):
                    triplePositives.append(myIndex)
            else:
                nonTripleNegatives.append(myIndex)
            # don't count for genes in skipTheseGenes
            if myRow['CarrierGene'] in skipTheseGenes:
                continue
            else:
                if myRow['CarrierGene']  not in countsPerGene.keys():
                    countsPerGene[myRow['CarrierGene']] = ReceptorCounts()
                countsPerGene[myRow['CarrierGene']].increment(myRow)

        # define fileObject based on config
        if myFDA.configFile.outputFile == "":
            fileObject = sys.stdout
        else:
            fileObject = open(myFDA.configFile.outputFile, mode='a')

        # print results to fileObject
        print("number of triple negatives = " + str(len(tripleNegatives)), file=fileObject)
        print('============================================')
        print("number of triple positives = " + str(len(triplePositives)), file=fileObject)
        print('============================================')
        print("number of triple values = " + str(len(tripleValues)), file=fileObject)
        print('============================================')
        print("number of non-triple negatives = " + str(len(nonTripleNegatives)), file=fileObject)
        print('============================================')
        for gene in list(myFDA.dataFile['CarrierGene'].unique()):
            if gene == 'NonCarrier':
                continue
            else:
                print("gene: " + gene)
                countsPerGene[gene].print()
                print('============================================')

        return True

    except Exception as e:
        print('exception in customDataAnalylzer.run() method: ' + str(e))
        return False


def isTripleNegative(myRow):
    # HER2 = 0, 1, or 2 => 'negative'; =3 => 'positive'?
    if myRow['ER'] == 'Negative' and myRow['PgR'] == 'Negative' and (myRow['HER2'] == '0'):
        return True
    else:
        return False

def isTripleValue(myRow):
    if myRow['ER'] != 'NA' and myRow['PgR'] != 'NA' and myRow['HER2'] != 'NA':
        return True
    else:
        return False

def isTriplePositive(myRow):
    # HER2 = 0, 1, or 2 => 'negative'; =3 => 'positive'?
    if myRow['ER'] == 'Positive' and myRow['PgR'] == 'Positive' and (myRow['HER2'] != '0'):
        return True
    else:
        return False

