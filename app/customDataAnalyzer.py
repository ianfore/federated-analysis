import sys
import json

class ReceptorCounts:
    """
    Purpose:
        Encapsulate count data for each of the receptors ER, PgR, and HER2.

    data members:
        counts: dictionary of counts for each (receptor, status) tuple
    """

    def __init__(self):
        self.counts = dict()
        self.counts["('ER', 'Positive')"] = 0
        self.counts["('ER', 'Negative')"] = 0
        self.counts["('PgR', 'Positive')"] = 0
        self.counts["('PgR', 'Negative')"] = 0
        self.counts["('HER2', '0')"] = 0
        self.counts["('HER2', '1+')"] = 0
        self.counts["('HER2', '2+')"] = 0
        self.counts["('HER2', '3+')"] = 0

    def increment(self, myRow):
        for receptor in [ 'ER', 'PgR', 'HER2']:
            if myRow[receptor] == 'NA':
                continue
            else:
                self.counts[str((receptor, myRow[receptor]))] += 1

    def print(self):
        print(json.dumps(self.__dict__))
        '''for status in self.counts.keys():
            print(str(status) + ': ' + str(self.counts[status]))'''


def run(myFDA):
    try:
        countsPerGene = dict()
        tripleNegatives = list()
        triplePositives = list()
        tripleValues = list()
        nonTripleNegatives = list()
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
            if myRow['CarrierGene'] == 'NonCarrier':
                continue
            else:
                if myRow['CarrierGene']  not in countsPerGene.keys():
                    countsPerGene[myRow['CarrierGene']] = ReceptorCounts()
                countsPerGene[myRow['CarrierGene']].increment(myRow)

        if myFDA.configFile.outputFile == "":
            fileObject = sys.stdout
        else:
            fileObject = open(myFDA.configFile.outputFile, mode='a')
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

        # for each gene, pgr status overall, <50, and >=50

        # for each gene, her2 grade overall , <50, and >=50


        return True


    except Exception as e:
        print(str(e))
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

