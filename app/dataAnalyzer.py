import sys
import statistics
import json
import pandas
import numpy
import app.customDataAnalyzer
import app.supplementaryTable4
import os

class ConfigFile:
    """
    Purpose:
        Encapsulate filter information from config file.

    data members
        fileName:       string name of data file
        fileHeader:     boolean header present in file
        fieldDelimiter: char that separates field in data file
        fieldFilters:    json object encapsulating field filters
        outputFile:     name of file to write results
    """

    def __init__(self, fileName, fileHeader, fieldDelimiter, fieldFilters, outputFile):
        self.fileName = fileName
        self.fileHeader = fileHeader
        self.fieldDelimiter = fieldDelimiter
        self.fieldFilters = fieldFilters
        self.outputFile = outputFile

        if self.outputFile != '' and os.path.exists(self.outputFile):
            print('output file ' + self.outputFile + ' exists!', file=sys.stderr)
            sys.exit(1)

class FieldCounter:
    """
    Purpose:
        class to store counts of fields in data file

    Data members:
        fieldCount:     dict of {field, count}

    """
    def __init__(self, fieldName):
        self.fieldCount = dict()
        self.incrementCounter(fieldName)

    def incrementCounter(self, fieldName):
        # increment count of field value
        if (fieldName not in self.fieldCount.keys()):
            self.fieldCount[fieldName] = 1
        else:
            self.fieldCount[fieldName] += 1

    def print(self):
        # print count of each value in field
        print(json.dumps(self.__dict__, sort_keys=True))

class FederatedDataAnalyzer:
    """
    Purpose:
        class to store and analyze federated data

    Data members:
        configFileName:     string full path to configuration file
        configFile:         object of type ConfigFile containing configuration
        dataFile:           object of type pandas.DataFrame containing data
        fieldFilters:       dictionary object mapping each field to their filters {field: filter}
        valueFrequency:     dictionary object mapping field name to value count {fieldName: FieldCount}
        badValues:          dictionary object tracking bad values {rowIndex: [list of tuples (fieldName:fieldValue)]}
        missingValues:      dictionary object tracking missing values {rowIndex: [list of missing values in row]}

    """
    def __init__(self, configFileName):
        self.configFileName = configFileName
        self.configFile, self.fieldFilters = self.readConfigFile()
        self.dataFile = self.readDataFile()
        self.fieldFilter = dict()
        for field in self.fieldFilters:
            self.fieldFilter[self.fieldFilters[field]['fieldName']] = self.fieldFilters[field]
        self.valueFrequency = dict()
        self.badValues = dict() #
        self.missingValues = dict()


    def readConfigFile(self):
        """
        Purpose:
            Read data from config file
        Input:
            fileName:     name of file where filter config is located
        Returns:
            data from filter file and the field filter (object of type FieldFilter)
        """
        with open(self.configFileName, 'r') as myFile:
            jsonData = myFile.read()
        data = json.loads(jsonData)

        configFile = ConfigFile(data['fileName'], data['fileHeader'], data['fieldDelimiter'], data['fieldFilters'],
                                data['outputFile'])
        fieldFilters = dict()
        for f in configFile.fieldFilters:
            fieldFilters[f['fieldName']] = f
        return configFile, fieldFilters

    def readDataFile(self):
        """
        Purpose:
            Read data from data file using config in config file
        Input:
            myConfigFile:     object of type ConfigFile
        Returns:
            data frame containing filtered data from data file
        """
        fieldFilters = self.configFile.fieldFilters
        fieldNames = list()
        for f in fieldFilters:
            fieldNames.append(f['fieldName'])
        if self.configFile.fileHeader:
            header = 0
        else:
            header = 1

        # read data from file into data frame and replace 'NaN' (literal string "NA") with string 'NA'
        # but keep blank fields ('') as null strings so we can differentiate b/w "NA" and ""
        return pandas.read_csv(self.configFile.fileName, sep=self.configFile.fieldDelimiter, header=header,
                               keep_default_na=False, na_values='NA')[fieldNames].replace(numpy.nan, 'NA', regex=True)

    def validateField(self, fieldValue, fieldFilter):
        """
        Purpose:
            function to validate the value of a field using the user-provided filter criteria
        Input:
            fieldValue:     int, float, string, ... literal
            fieldFilter:    filter for this field , parsed from conf file
        Returns:
            True (valid) or False (not valid)
        """
        # if field is categorical and there are categories to match, test categories.
        if fieldFilter['fieldType'] == 'categorical':
            if len(fieldFilter['fieldValues']) == 0 or fieldValue in fieldFilter['fieldValues']:
                return True
            else:
                return False
        # else if field is numerical and data is int or float, then test if values match.
        elif fieldFilter['fieldType'] == 'numerical':
            if isinstance(fieldValue, (int)) or isinstance(fieldValue, (float)):
                if len(fieldFilter['fieldValues']) == 0 or fieldValue in fieldFilter['fieldValues']:
                    return True
                else:
                    print(fieldFilter['fieldType'] + ' is bad first')
                    return False
            elif isinstance(fieldValue, (str)):
                # try to cast it to float (floats and ints will cast to float)
                try:
                    fieldValue = float(fieldValue.strip("'"))
                except Exception as e:
                    return False
                return True
            else:
                return False
        # else if field is free-form and data is a string, then test character encodings of string.
        elif fieldFilter['fieldType'] == 'free-form':
            if isinstance(fieldValue, (str)):
                for charSet in fieldFilter['fieldValues']:
                    try:
                        fieldValue.encode(charSet)
                    except:
                        continue
                    else:
                        return True
                return False
            else: # shouldn't get here, but just in case :)
                return False
        # you'll reach here if it's not a valid field type as defined in filter file (TODO consider throwing exception here)
        else:
            return False


    def getStatistics(self, valueFrequency, myField):
        """
        Purpose:
            function to get basic stats for a given field
        Input:
            valueFrequency:     dictionary of {value: count}
            myField:              string name of field
        Returns:
            min, max, mean, median of values in myCol field
        """
        # create a list of values from the keys in the valueFrequency dict, and calculate stats from list
        allVals = list()
        for x in valueFrequency[myField].fieldCount.keys():
            for n in range(0, valueFrequency[myField].fieldCount[x]):
                if isinstance(x, (int)) or isinstance(x, (float)):
                    allVals.append(x)
                elif isinstance(x, (str)):
                    try:
                        x = float(x.strip("'"))
                        allVals.append(x)
                    except:
                        continue
        return min(allVals), max(allVals), statistics.mean(allVals), statistics.median(allVals)


    def printResults(self):
        """
        Purpose:
            Prints the values to sys.stdout
        Input:
            self
        Returns:
            void
        """

        if self.configFile.outputFile == "":
            fileObject = sys.stdout
        else:
            fileObject = open(self.configFile.outputFile, mode='a')

        print("============================================", file=fileObject)
        print('config file info: ' + str(json.dumps(self.configFile.__dict__)), file=fileObject)
        print("============================================", file=fileObject)
        print('total records read from data file: ' + str(len(self.dataFile)), file=fileObject)
        print("============================================", file=fileObject)
        for myField in self.valueFrequency.keys():
            print('column: ' + myField + ' / type: ' + self.fieldFilter[myField]['fieldType'], file=fileObject)
            print(json.dumps(self.valueFrequency[myField].__dict__, sort_keys=True), file=fileObject)
            if self.fieldFilter[myField]['fieldType'] == 'numerical':
                minVal, maxVal, meanVal, medianVal = self.getStatistics(self.valueFrequency, myField)
                print('min = ' + str(minVal) + ', max = ' + str(maxVal) + ', mean = ' + str(meanVal) + ' median = ' +
                      str(medianVal), file=fileObject)
            print("============================================", file=fileObject)
        print('bad values: ' + str(self.badValues), file=fileObject)
        print("============================================", file=fileObject)
        print('missing values: ' + str(self.missingValues), file=fileObject)
        print("============================================", file=fileObject)


    def run(self):
        """
        Purpose:
            This is the big loop that iterates through all rows and performs analysis on all fields of interest
        Input:
            self
        Returns:
            void
        """
        try:
            for myIndex, myRow in self.dataFile.iterrows():
                for myField in self.dataFile.columns:
                    # check if field has a value
                    if myRow[myField] == '':
                        if myIndex not in self.missingValues.keys():
                            self.missingValues[myIndex] = list()
                        self.missingValues[myIndex].append(myField)
                    # validate field value
                    elif(not self.validateField(myRow[myField], self.fieldFilters[myField])):
                        if myIndex not in self.badValues.keys():
                            self.badValues[myIndex] = list()
                        self.badValues[myIndex].append((myField, myRow[myField]))
                    else:
                        # add +1 to value counter
                        if myField not in self.valueFrequency.keys():
                            self.valueFrequency[myField] = FieldCounter(myRow[myField])
                        else:
                            self.valueFrequency[myField].incrementCounter(myRow[myField])
            # print data summary to stdout
            self.printResults()
            return True

        except Exception as e:
            print('exception in dataAnalyzer.run() method: ' + str(e))
            return False

def main():

    # get file name of config file from command line
    configFileName = sys.argv[1]

    # instantiate analyzer object
    myFederatedDataAnalyzer = FederatedDataAnalyzer(configFileName)

    # run analyzer and any custom code
    return myFederatedDataAnalyzer.run() and app.customDataAnalyzer.run(myFederatedDataAnalyzer) \
            and app.supplementaryTable4.run(myFederatedDataAnalyzer)

if __name__ == "__main__":
    main()
