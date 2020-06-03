import pandas as pd
import sys
import json


def main():
    if len(sys.argv) != 4:
        print('ancestry.tsv cohort.txt outputDir')
        sys.exit(1)

    ancestryFileName = sys.argv[1]
    cohortFileName = sys.argv[2]
    outputDir = sys.argv[3]

    ancestry = pd.read_csv(ancestryFileName, sep='\t')

    cohort = pd.read_csv(cohortFileName, header=None)

    populationPerIndividual = dict()
    for individual in cohort[0]:
        row = ancestry.loc[ancestry['individual'] == individual]
        # row = Index(['individual', 'Sub_Saharan_Africa', 'Central_and_South_Asia',
        #       'East_Asia', 'Europe', 'Native_America', 'Oceania', 'Middle_East'],
        #       dtype='object')
        tempPop = None
        tempMax = 0.0
        for pop in row.columns[1:]:
            # pop = 107046    0.7787
            #       Name: Sub_Saharan_Africa, dtype: float64
            try:
                if float(row[pop].values[0]) > tempMax:
                    tempPop = row[pop].name
                    tempMax = row[pop].values[0]
            except Exception as e:
                continue
        populationPerIndividual[individual] = (tempPop, tempMax)

    with open(outputDir + '/ancestries.json', 'w') as f:
        json.dump(populationPerIndividual)
    f.close()

if __name__ == "__main__":
    main()