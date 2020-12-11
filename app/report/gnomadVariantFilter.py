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


'''TODO  genome_ac_hom_delta and exome_ac_hom_delta, which describe the number of 
homozygous observations in gnomAD that came from topmed for the genome and exome data respectively. Take
their sum.'''


def main():
    if len(sys.argv) != 6:
        print('vcf-input brca-input mc-file notinsubset-file inubset-file')
        sys.exit(1)

    vcfFileName = sys.argv[1]
    brcaFileName = sys.argv[2]
    mcFileName = sys.argv[3]
    notInSubsetFileName = sys.argv[4]
    inSubsetFileName = sys.argv[5]

    logger.info('finding variants from ' + brcaFileName)
    brcaDF = findVariantsInBRCA(brcaFileName)

    logger.info('finding variants from ' + mcFileName)
    mcDF = pd.read_csv(mcFileName, delimiter='\t', header=0, dtype=str)

    logger.info('reading VCF file ' + vcfFileName)
    t = time.time()
    vcf = readVCFFile(vcfFileName)
    logger.info('elapsed time in readVCFFile() ' + str(time.time() - t))

    effectivelyZeroValues = ['0', '0.0', '-', None]

    notInSubset = open(notInSubsetFileName, 'w')
    inSubset = open(inSubsetFileName, 'w')
    for variant in range(len(vcf['calldata/GT'])):
        c = int(vcf['variants/CHROM'][variant].replace('chr', ''))
        p = int(vcf['variants/POS'][variant])
        r = str(vcf['variants/REF'][variant])
        a = str(vcf['variants/ALT'][variant][0])
        logger.debug((c,p,r,a))
        #if not checkGnomad(brcaDF, (c,p,r,a), 38, keyStrings):
        #    f.write('chr' + str(c) + '\t' + str(p) + '\t' + str(r) + '\t' + str(a) + '\n')
        exomeDelta, genomeDelta, hg37 = checkMelissaTable((c,p,r,a), mcDF, brcaDF)
        if exomeDelta in effectivelyZeroValues and genomeDelta in effectivelyZeroValues:
            notInSubset.write('(' + str(c) + ',' + str(p) + ',' + str(r) + ',' + str(a)  + ')' + '\n')
        else:
            ed = 0.0
            gd = 0.0
            try:
                ed = float(exomeDelta)
            except:
                pass
            try:
                gd = float(genomeDelta)
            except:
                pass
            deltaSum = ed + gd
            print(deltaSum)
            inSubset.write('(' + str(c) + ',' + str(p) + ',' + str(r) + ',' + str(a)  + ')' + '\n')

    inSubset.close()
    notInSubset.close()

def checkMelissaTable(variant, mcDF, brcaDF):
    chrom = variant[0]
    pos = variant[1]
    ref = variant[2]
    alt = variant[3]
    hgString = 'chr' + str(chrom) + ':g.' + str(pos) + ':' + str(ref) + '>' + str(alt)
    row = brcaDF[brcaDF[coordinateColumnBase + str(hgVersion)] == hgString]
    if len(row) == 0:
        return None, None, 0
    coord = str(row['Genomic_Coordinate_hg37'])
    pos_37 = int(coord.split(':')[1].split('.')[1])
    row = mcDF[(mcDF['chrom'] == str(chrom)) & (mcDF['pos'] == str(pos_37)) & (mcDF['ref'] == ref) & (mcDF['alt'] == alt)]
    if len(row) == 0:
        return None, None, pos_37
    exomeDelta = row['exome_ac_hom_delta'].iloc[0]
    genomeDelta = row['genome_ac_hom_delta'].iloc[0]

    return exomeDelta, genomeDelta, pos_37



def generateKeyStrings():
    # Allele_count_hom_genome_AFR_GnomAD
    ethnicities = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS']
    prefix = 'Allele_count_hom'
    suffix = '_GnomAD'
    omes = ['_genome_', '_exome_']
    keyStrings = list()
    for e in ethnicities:
        for o in omes:
            keyStrings.append(prefix + o + e + suffix)

    return keyStrings

def readVCFFile(vcfFileName):
    return allel.read_vcf(vcfFileName, fields=['samples', 'calldata/GT', 'variants/ALT', 'variants/CHROM',
                'variants/FILTER_PASS', 'variants/ID', 'variants/POS', 'variants/QUAL','variants/REF', 'variants/INFO'])


def findVariantsInBRCA(fileName):
    return pd.read_csv(fileName, sep='\t', header=0, dtype=str)


def checkGnomad(brcaDF, vus, hgVersion, keyStrings):
    # TODO write a unit test
    # 13:g.32393468:C>CT
    #hgString = 'chr' + str(vus[0][0]) + ':g.' + str(vus[0][1]) + ':' + str(vus[0][2]) + '>' + str(vus[0][3])
    if type(vus) is str:
        vus = eval(vus)
    hgString = 'chr' + str(vus[0]) + ':g.' + str(vus[1]) + ':' + str(vus[2]) + '>' + str(vus[3])

    retValue = getInGnomad(brcaDF, hgString, hgVersion, keyStrings)
    if retValue:
        return True
    else:
        return False


def getInGnomad(brcaDF, hgString, hgVersion, alleleFrequencies):
    for af in alleleFrequencies:
        freq = None
        alleleFreqList = brcaDF[brcaDF[coordinateColumnBase + str(hgVersion)] == hgString][af].tolist()
        if alleleFreqList:
            try:
                if len(alleleFreqList) > 1:
                    print('more than one')
                # yes, it always returns the first elt of list
                freq = float(alleleFreqList[0])
                return freq
            except ValueError:
                continue
    return freq


if __name__ == "__main__":
    main()