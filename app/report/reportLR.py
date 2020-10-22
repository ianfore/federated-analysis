import pandas as pd
import json
import sys
import logging
import pyhgvs as hgvs

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

def main():
    if(len(sys.argv) != 5):
        print('usage: bayesdel.vcf intersection.json out.json to-translate.txt')
        sys.exit(1)

    variantsFileName = sys.argv[1]
    intersectionFileName = sys.argv[2]
    outFileName = sys.argv[3]

    logger.info('reading ' + variantsFileName)
    variantsDF = pd.read_csv(variantsFileName, sep='\t', header=0)

    logger.info('reading ' + intersectionFileName)
    with open(intersectionFileName, 'r') as f:
        intersectionDict = json.load(f)
    f.close()

    logger.info('reading ' + outFileName)
    with open(outFileName, 'r') as f:
        outDict = json.load(f)
    f.close()

    for vus in outDict['cooccurring vus']:
        vus = eval(vus)
        chrom = vus[0]
        pos = vus[1]
        ref = vus[2]
        alt = vus[3]
        hgvsCoord = vus[4]
        infoField = tuple(variantsDF[variantsDF['POS'] == pos]['INFO'])
        try:
            bayesdelScore = infoField[0].split('BayesDel')[1].split('=')[1]
        except Exception as e:
            continue
        print(bayesdelScore)

    for vus in outDict['homozygous vus']:
        chrom = vus[0]
        pos = vus[1]
        ref = vus[2]
        alt = vus[3]
        hgvsCoord = vus[4]
        infoField = variantsDF[variantsDF['POS']==pos]['INFO']
        try:
            bayesdelScore = float(infoField.iloc[1].split('=')[3])
        except Exception as e:
            continue
        print(bayesdelScore)



if __name__ == "__main__":
    main()