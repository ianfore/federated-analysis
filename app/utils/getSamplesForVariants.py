import json
import argparse


def main():
    parser = argparse.ArgumentParser(usage="getSamplesForReads")
    parser.add_argument("--v", dest="variantFileName", help="variant file name", default=None)
    parser.add_argument("--i", dest="ipvFileName", help="ipv file name", default=None)
    parser.add_argument("--o", dest="outputFileName", help="output file name", default=None)
    options = parser.parse_args()
    variantFileName = options.variantFileName
    ipvFileName = options.ipvFileName
    outputFileName = options.outputFileName

    variantFile = open(variantFileName)
    variants = variantFile.readlines()
    variantFile.close()

    with open(ipvFileName, 'r') as f:
        ipv = json.load(f)
    f.close()

    for variant in variants:
        variant = variant.replace('(', '').replace(')', '')
        variantArray = variant.split(',')
        chrom = int(variantArray[0])
        pos = int(variantArray[1])
        ref = str(variantArray[2])
        alt = str(variantArray[3].replace('\n', ''))
        variantString = str((chrom, pos, ref, alt))
        print(variantString)




if __name__ == "__main__":
    main()