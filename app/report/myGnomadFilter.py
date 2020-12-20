import pandas as pd
import logging
import json
from pyliftover import LiftOver
import argparse

coordinateColumnBase = 'Genomic_Coordinate_hg'

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)


'''TODO  genome_ac_hom_delta and exome_ac_hom_delta, which describe the number of 
homozygous observations in gnomAD that came from topmed for the genome and exome data respectively. Take
their sum.'''

# vcf-hgversion and gnomad-hgversion must of the form hg38, hg19, ...
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--homozygous', help='only homozygous variants input file')
    parser.add_argument('-e', '--heterozygous', help='only heterozygous variants input file')
    parser.add_argument('-b', '--both', help='both homo- and hetero-zygous variants input file')
    parser.add_argument('-m', '--mcfile', help='melissas gnomad input file')
    parser.add_argument('-n', '--notinfile', help='not in gnomad output file')
    parser.add_argument('-i', '--infile', help='in gnomad output file')
    parser.add_argument('-v', '--version', help='hg version of input files')
    parser.add_argument('-g', '--gnomadversion', help='hg version of gnomad')
    options = parser.parse_args()
    return options


def main():
    onlyHomozygousVariantsFile = parse_args().homozygous
    onlyHeterozygousVariantsFile = parse_args().heterozygous
    bothVariantsFile = parse_args().both
    mcFile = parse_args().mcfile
    notInSubsetFile = parse_args().notinfile
    inSubsetFile = parse_args().infile
    inputHG = parse_args().version
    gnomadHG = parse_args().gnomadversion

    logger.info('finding variants from ' + mcFile)
    mcDF = pd.read_csv(mcFile, delimiter='\t', header=0, dtype=str)

    try:
        lo = LiftOver(inputHG, gnomadHG)
    except Exception as e:
        lo = None

    processVariants(onlyHomozygousVariantsFile, inSubsetFile +'_homo.txt', notInSubsetFile + '_homo.txt', mcDF, 'homo', lo)
    processVariants(onlyHeterozygousVariantsFile, inSubsetFile + '_hetero.txt', notInSubsetFile + '_hetero.txt', mcDF, 'hetero', lo)
    processVariants(bothVariantsFile, inSubsetFile + '_both.txt', notInSubsetFile + '_both.txt', mcDF, 'both', lo)


def processVariants(variantsFile, inSubsetFile, notInSubsetFile, mcDF, fileType, lo):
    effectivelyZeroValues = ['0', '0.0', '-', None]

    variants = open(variantsFile, 'r')
    notInSubset = open(notInSubsetFile, 'w')
    inSubset = open(inSubsetFile, 'w')

    for line in variants.readlines():
        variant = eval(line)
        c = int(variant[0])
        p = int(variant[1])
        r = str(variant[2])
        a = str(variant[3])

        exomeDelta, genomeDelta, hg = checkMelissaTable((c,p,r,a), mcDF, lo, fileType)
        if exomeDelta in effectivelyZeroValues and genomeDelta in effectivelyZeroValues:
            notInSubset.write('(' + str(c) + ',' + str(p) + ',' + str(r) + ',' + str(a)  + ')' + '\n')
        else:
            ed = 0.0
            gd = 0.0
            try:
                ed = float(exomeDelta)
            except:
                logger.error('couldnt convert exomedelta to float')
                pass
            try:
                gd = float(genomeDelta)
            except:
                logger.error('couldnt convert genomedelta to float')
                pass
            deltaSum = ed + gd
            logger.debug(deltaSum)
            if deltaSum in effectivelyZeroValues:
                notInSubset.write('(' + str(c) + ',' + str(p) + ',' + str(r) + ',' + str(a) + ')' + '\n')
            else:
                inSubset.write('(' + str(c) + ',' + str(p) + ',' + str(r) + ',' + str(a)  + ')' + '\n')
    inSubset.close()
    notInSubset.close()

def checkMelissaTable(variant, mcDF, lo, fileType):
    chrom = variant[0]
    pos = variant[1]
    ref = variant[2]
    alt = variant[3]

    if not lo is None:
        coord = lo.convert_coordinate('chr' + str(chrom), int(pos))
        if coord is None or len(coord) != 1:
            return None, None, None
        try:
            pos = int(coord[0][1])
        except Exception as e:
            return None, None, pos
    row = mcDF[(mcDF['chrom'] == str(chrom)) & (mcDF['pos'] == str(pos)) & (mcDF['ref'] == ref) & (mcDF['alt'] == alt)]
    if len(row) == 0:
        return None, None, pos

    if fileType == 'homo':
        exomeDelta = row['exome_ac_hom_delta'].iloc[0]
        genomeDelta = row['genome_ac_hom_delta'].iloc[0]
    elif fileType == 'hetero':
        exomeDelta = row['exome_ac_delta'].iloc[0]
        genomeDelta = row['genome_ac_delta'].iloc[0]
    elif fileType == 'both':
        exomeDelta = row['exome_ac_delta'].iloc[0]
        genomeDelta = row['genome_ac_delta'].iloc[0]
    else:
        return None, None, None
    return exomeDelta, genomeDelta, pos


if __name__ == "__main__":
    main()