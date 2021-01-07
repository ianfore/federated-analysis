import pandas as pd
import sys
import json


def main():
    if len(sys.argv) != 3:
        print('ancestry.tsv outputDir')
        sys.exit(1)

    ancestryFileName = sys.argv[1]
    outputDir = sys.argv[2]
    ancestry = pd.read_csv(ancestryFileName, sep='\t', header=0)
    topmed2gnomAD = dict()
    topmed2gnomAD['Sub_Saharan_Africa'] = 'AFR'
    topmed2gnomAD['Central_and_South_Asia'] = 'SAS'
    topmed2gnomAD['East_Asia'] = 'EAS'
    topmed2gnomAD['Europe'] = 'NFE'
    topmed2gnomAD['Native_America'] = 'AMR'
    topmed2gnomAD['Oceania'] = 'OTH'
    topmed2gnomAD['Middle_East'] = 'OTH'

    populationPerIndividual = dict()
    for individual in ancestry['individual']:
        row = ancestry[ancestry['individual'] == individual]
        # row = Index(['individual', 'Sub_Saharan_Africa', 'Central_and_South_Asia',
        #       'East_Asia', 'Europe', 'Native_America', 'Oceania', 'Middle_East'],
        #       dtype='object')
        tempPop = None
        tempMax = 0.0
        for pop in row.columns[1:7]:
            # pop = 107046    0.7787
            #       Name: Sub_Saharan_Africa, dtype: float64
            try:
                if float(row[pop].values[0]) > tempMax:
                    tempPop = row[pop].name
                    tempMax = row[pop].values[0]
            except Exception as e:
                continue
        if not individual in populationPerIndividual:
            populationPerIndividual[individual] = dict()
        populationPerIndividual[individual]['topmedPop'] = (tempPop, tempMax)
        try:
            populationPerIndividual[individual]['gnomadPop'] = topmed2gnomAD[tempPop]
        except:
            populationPerIndividual[individual]['gnomadPop'] = 'OTH'

    with open(outputDir + '/ancestries.json', 'w') as f:
        json.dump(populationPerIndividual, f)
    f.close()


if __name__ == "__main__":
    main()