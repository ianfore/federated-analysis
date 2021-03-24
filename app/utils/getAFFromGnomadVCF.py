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

    fields = vcfDF['INFO'][21].split(';')

    variantsDict = dict()

    for field in fields:
        if '=' in field:
            key = field.split('=')[0]
            value = field.split('=')[1]
            variantsDict[key] = value

    with open(outputFileName, 'w') as f:
        json.dump(variantsDict)

    f.close()






if __name__ == "__main__":
    main()