import pandas as pd
import logging
import argparse
import json
import matplotlib.pyplot as plt

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
    topmedKeys = ['AF-afr', 'AF-amr', 'AF-nfe']
    nontopmedKeys = ['AF-non_topmed-afr', 'AF-non_topmed-amr', 'AF-non_topmed-nfe']
    keys = topmedKeys + nontopmedKeys

    logger.info('getting allele freqs')
    for i in range(len(vcfDF)):
        fields = vcfDF.iloc[i]['INFO'].split(';')
        chrom = vcfDF.iloc[i]['#CHROM'].split('chr')[1]
        pos = vcfDF.iloc[i]['POS']
        ref = vcfDF.iloc[i]['REF']
        alt = vcfDF.iloc[i]['ALT']
        mykey = "(" + chrom + ", " + pos + ", '" + ref + "', '" + alt + "')"
        variantsDict[mykey] = dict()
        for field in fields:
            if '=' in field:
                key = field.split('=')[0]
                if key in keys:
                    value = field.split('=')[1]
                    variantsDict[mykey][key] = value

    logger.info('saving allele freqs')
    with open(outputFileName, 'w') as f:
        json.dump(variantsDict, f)
    f.close()

    logger.info('plotting allele freqs')
    plotProbability(variantsDict, topmedKeys, nontopmedKeys)

def plotProbability(variantsDict, topmedKeys, nontopmedKeys):

    topmedKeys.sort()
    nontopmedKeys.sort()
    topmedList = list()
    nontopmedList = list()
    for variant in variantsDict:
        topmedSum = 0
        nontopmedSum = 0
        for key in topmedKeys:
            topmedSum += float(variantsDict[variant][key])
        topmedList.append(topmedSum)
        for key in nontopmedKeys:
            nontopmedSum += float(variantsDict[variant][key])
        nontopmedList.append(nontopmedSum)

    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.plot(nontopmedList, topmedList, marker='.', color='red')
    plt.ylabel('topmed AF', fontsize=18)
    plt.xlabel('nontopmed AF', fontsize=18)
    plt.title('nontopmed vs topmed AF')

    plt.savefig('qq.png')
    #plt.show()
    plt.close()


if __name__ == "__main__":
    main()