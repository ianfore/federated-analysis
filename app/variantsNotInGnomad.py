import allel
import sys
import pandas as pd
import logging
import time
import json

coordinateColumnBase = 'Genomic_Coordinate_hg'
hgVersion = 38

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)



def main():
    if len(sys.argv) != 4:
        print('vcf-input brca-input output-file')
        sys.exit(1)

    vcfFileName = sys.argv[1]
    brcaFileName = sys.argv[2]
    outputFile = sys.argv[3]

    logger.info('finding variants from ' + brcaFileName)
    brcaDF = findVariantsInBRCA(brcaFileName)

    logger.info('reading VCF file ' + vcfFileName)
    t = time.time()
    vcf = readVCFFile(vcfFileName)
    logger.info('elapsed time in readVCFFile() ' + str(time.time() - t))

    varsNotInGnomad = list()

    for variant in range(len(vcf['calldata/GT'])):
        if 1 in vcf['calldata/GT'][variant][0]:
            c = int(vcf['variants/CHROM'][variant].replace('chr', ''))
            p = int(vcf['variants/POS'][variant])
            r = str(vcf['variants/REF'][variant])
            a = str(vcf['variants/ALT'][variant][0])
            if not checkGnomad(brcaDF, (c,p,r,a), 38):
                varsNotInGnomad.append((c,p,r,a))

    with open(outputFile, 'r') as f:
        json.dump(varsNotInGnomad, f)
    f.close()

def readVCFFile(vcfFileName):
    return allel.read_vcf(vcfFileName, fields=['samples', 'calldata/GT', 'variants/ALT', 'variants/CHROM',
                'variants/FILTER_PASS', 'variants/ID', 'variants/POS', 'variants/QUAL','variants/REF', 'variants/INFO'])


def findVariantsInBRCA(fileName):
    return pd.read_csv(fileName, sep='\t', header=0, dtype=str)


def checkGnomad(brcaDF, vus, hgVersion):
    # TODO write a unit test
    # 13:g.32393468:C>CT
    #hgString = 'chr' + str(vus[0][0]) + ':g.' + str(vus[0][1]) + ':' + str(vus[0][2]) + '>' + str(vus[0][3])
    if type(vus) is str:
        vus = eval(vus)
    hgString = 'chr' + str(vus[0]) + ':g.' + str(vus[1]) + ':' + str(vus[2]) + '>' + str(vus[3])

    # first, get list of columns for GnomAD allleles
    #alleleFrequencyPrefixes = ['Allele_frequency_genome_', 'Allele_frequency_']
    # 'Allele_count_exome_AMR_GnomAD', 'Allele_count_hemi_exome_AMR_GnomAD', 'Allele_count_hom_exome_AMR_GnomAD',
    # 'Allele_number_exome_AMR_GnomAD', 'Allele_frequency_exome_AMR_GnomAD' ,'Allele_frequency_AMR_GnomAD'

    gnomad = [v for v in list(brcaDF.columns) if 'GnomAD' in v]

    # second, get frequencies across exomes and genomes to determine max
    # return population, frequency, count, and number
    # replace "frequency" with "count" and "number" in Allele_frequency_genome_AFR_GnomAD

    allDict = dict()

    alleleFrequencies = [v for v in gnomad if 'Allele_frequency' in v]
    allDict['max'] = getMaxGnomad(brcaDF, hgString, hgVersion, alleleFrequencies)
    if allDict['max']:
        return True
    else:
        return False


def getMaxGnomad(brcaDF, hgString, hgVersion, alleleFrequencies):
    maxData = {'frequency': 0.0, 'population': None}
    for af in alleleFrequencies:
        freq = 0.0
        alleleFreqList = brcaDF[brcaDF[coordinateColumnBase + str(hgVersion)] == hgString][af].tolist()
        if alleleFreqList:
            try:
                if len(alleleFreqList) > 1:
                    print('more than one')
                # yes, it always returns a list of length 0
                freq = float(alleleFreqList[0])
            except ValueError:
                continue
            if freq > maxData['frequency']:
                # TODO get homo and abs counts as well
                maxData['frequency'] = freq
                maxData['population'] = af
    return maxData


if __name__ == "__main__":
    main()