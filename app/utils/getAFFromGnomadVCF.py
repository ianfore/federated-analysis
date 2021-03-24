import pandas as pd
import logging
import argparse
import json


logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', help='vcf input file')
    parser.add_argument('-o', '--output', help='json output file')
    options = parser.parse_args()
    return options


def main():
    vcfFileName = parse_args().vcf
    outputFileName = parse_args().output

    logger.info('finding variants from ' + vcfFileName)
    vcfDF = pd.read_csv(vcfFileName, delimiter='\t', header=0, dtype=str)
    variantsDict = dict()

    keys = ['AF-non_topmed-afr', 'AF-afr', 'AF-non_topmed-amr', 'AF-amr',
            'AF-non_topmed-nfe', 'AF-nfe']
    for i in range(len(vcfDF)):
        fields = vcfDF.iloc[i]['INFO'].split(';')
        chrom = vcfDF.iloc[i]['#CHROM']
        pos = vcfDF.iloc[i]['POS']
        ref = vcfDF.iloc[i]['REF']
        alt = vcfDF.iloc[i]['ALT']
        mykey = "(" + chrom + ", " + pos + ", " + ref + ", " + alt + ")"
        variantsDict[mykey] = dict()
        for field in fields:
            if '=' in field:
                key = field.split('=')[0]
                if key in keys:
                    value = field.split('=')[1]
                    variantsDict[mykey][str(i) + str(key)] = str(value)

    with open(outputFileName, 'w') as f:
        json.dump(variantsDict, f)

    f.close()






if __name__ == "__main__":
    main()