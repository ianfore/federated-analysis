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

    reportDict = dict()

    for vus in outDict['homozygous vus']:
        vus = eval(vus)
        chrom = vus[0]
        pos = vus[1]
        ref = vus[2]
        alt = vus[3]
        hg38Coord = 'chr' + str(chrom) + ':' + str(pos) + ':' + ref + '>' + alt
        hsgvsCoord = variantsDF[variantsDF['Genomic_Coordinate_hg38'] == hg38Coord]['pyhgvs_cDNA']
        bayesDel = variantsDF[variantsDF['Genomic_Coordinate_hg38'] == hg38Coord]['BayesDel_nsfp33a_noAF']
        print(hsgvsCoord)
        print(bayesDel)




if __name__ == "__main__":
    main()