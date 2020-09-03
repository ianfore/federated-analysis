import pandas as pd
import sys
import json


def main():
    if len(sys.argv) != 3:
        print('input-dir output-dir')
        sys.exit(1)

    inputDir =sys.argv[1]
    outputDir = sys.argv[2]
    df1 = pd.read_csv(inputDir + '/f8_chr17_brca1_gruhmb_report.tsv', sep='\t')
    df2 = pd.read_csv(inputDir + '/f8_chr13_brca2_gruhmb_report.tsv', sep='\t')

    # read .txt files
    onlyHomo_13_file = open(inputDir + "/13-only-homo.txt", "r")
    onlyHomo_13 = onlyHomo_13_file.readlines()
    onlyHomo_17_file = open(inputDir + "/17-only-homo.txt", "r")
    onlyHomo_17 = onlyHomo_17_file.readlines()

    onlyCooc_13_file = open(inputDir + "/13-only-coocs.txt", "r")
    onlyCoocs_13 = onlyCooc_13_file.readlines()
    onlyCooc_17_file = open(inputDir + "/17-only-coocs.txt", "r")
    onlyCoocs_17 = onlyCooc_17_file.readlines()

    both_13_file = open(inputDir + "/13-both.txt", "r")
    both_13 = both_13_file.readlines()
    both_17_file = open(inputDir + "/17-both.txt", "r")
    both_17 = both_17_file.readlines()

    cols = ['variant', 'popFreq', 'cohortFreq', 'homo_alt', 'hetero', 'homo_ref', 'p', 'q', 'exonic']
    eggReport_13 = pd.DataFrame(columns = cols)
    varTypeList = list()
    for i in range(len(df2)):
        j = 0
        if df2.iloc[i]['variant'] in onlyHomo_13:
            eggReport_13.iloc[j] = df2.iloc[i][cols]
            varTypeList.append('homozygous')
            j += 1
        elif df2.iloc[i]['variant'] in onlyCoocs_13:
            eggReport_13.iloc[j] = df2.iloc[i][cols]
            varTypeList.append('cooccurring')
            j += 1
        elif df2.iloc[i]['variant'] in both_13:
            eggReport_13.iloc[j] = df2.iloc[i][cols]
            varTypeList.append('both')
            j += 1

    eggReport_13.assign(varType = varTypeList)
    print(eggReport_13)



