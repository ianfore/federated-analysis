import pandas as pd
import sys
import json


def main():
    if len(sys.argv) != 4:
        print('input-dir brca-variants.tsv output-dir')
        sys.exit(1)

    inputDir =sys.argv[1]
    outputDir = sys.argv[2]
    brcaVariantsFile = sys.argv[3]
    df1 = pd.read_csv(inputDir + '/f8_chr17_brca1_gruhmb_report.tsv', sep='\t')
    df2 = pd.read_csv(inputDir + '/f8_chr13_brca2_gruhmb_report.tsv', sep='\t')
    brcaDF = pd.read_csv(brcaVariantsFile, sep='\t', header=0)

    # read .txt files
    onlyHomo_13_file = open(inputDir + "/13-only-homo.txt", "r")
    onlyHomo_13 = onlyHomo_13_file.read().splitlines()
    onlyHomo_17_file = open(inputDir + "/17-only-homo.txt", "r")
    onlyHomo_17 = onlyHomo_17_file.read().splitlines()

    onlyCooc_13_file = open(inputDir + "/13-only-coocs.txt", "r")
    onlyCoocs_13 = onlyCooc_13_file.read().splitlines()
    onlyCooc_17_file = open(inputDir + "/17-only-coocs.txt", "r")
    onlyCoocs_17 = onlyCooc_17_file.read().splitlines()

    both_13_file = open(inputDir + "/13-both.txt", "r")
    both_13 = both_13_file.read().splitlines()
    both_17_file = open(inputDir + "/17-both.txt", "r")
    both_17 = both_17_file.read().splitlines()

    cols = ['variant', 'popFreq', 'cohortFreq', 'homo_alt', 'hetero', 'homo_ref', 'p', 'q', 'exonic', 'inGnomad']
    eggReport_13 = pd.DataFrame(columns = cols)
    varTypeList = list()
    otherClassList = list()
    for i in range(len(df2)):
        otherClass = "NA"
        x = df2.iloc[i]['variant']
        if x in onlyHomo_13:
            eggReport_13 = eggReport_13.append(df2.iloc[i][cols])
            varTypeList.append('homozygous')
        elif x in onlyCoocs_13:
            eggReport_13 = eggReport_13.append(df2.iloc[i][cols])
            varTypeList.append('cooccurring')
        if x in both_13:
            eggReport_13 = eggReport_13.append(df2.iloc[i][cols])
            varTypeList.append('both')
        pos = x[1]
        try:
            otherClass = brcaDF.loc[brcaDF['Pos'] == pos]['Pathogenicity_all']
        except Exception as e:
            pass
        otherClassList.append(otherClass)

    eggReport_13['varType'] = varTypeList
    eggReport_13.to_csv(outputDir + '/f8-brca2-forEGG.tsv', index=False, header=True, sep='\t')
    eggReport_13['otherClass'] = otherClassList

    eggReport_17 = pd.DataFrame(columns = cols)
    varTypeList = list()
    otherClassList = list()
    for i in range(len(df1)):
        otherClass = "NA"
        x = df1.iloc[i]['variant']
        if x in onlyHomo_17:
            eggReport_17 = eggReport_17.append(df1.iloc[i][cols])
            varTypeList.append('homozygous')
        elif x in onlyCoocs_17:
            eggReport_17 = eggReport_17.append(df1.iloc[i][cols])
            varTypeList.append('cooccurring')
        if x in both_17:
            eggReport_17 = eggReport_17.append(df1.iloc[i][cols])
            varTypeList.append('both')
        pos = x[1]
        try:
            otherClass = brcaDF.loc[brcaDF['Pos'] == pos]['Pathogenicity_all']
        except Exception as e:
            pass
        otherClassList.append(otherClass)

    eggReport_17['varType'] = varTypeList
    eggReport_17.to_csv(outputDir + '/f8-brca1-forEGG.tsv', index=False, sep='\t', header=True)
    eggReport_17['otherClass'] = otherClassList

if __name__ == "__main__":
    main()



