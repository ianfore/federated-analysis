import json
import argparse


def main():
    parser = argparse.ArgumentParser(usage="getSamplesForReads")
    parser.add_argument("--v", dest="variantFileName", help="variant file name", default=None)
    parser.add_argument("--i", dest="ipvFileName", help="ipv file name", default=None)
    parser.add_argument("--g", dest="gen3FileName", help="gen3 file name", default=None)
    parser.add_argument("--o", dest="outputFileName", help="output file name", default=None)
    options = parser.parse_args()
    variantFileName = options.variantFileName
    ipvFileName = options.ipvFileName
    gen3FileName = options.gen3FileName
    outputFileName = options.outputFileName

    variantFile = open(variantFileName)
    variants = variantFile.readlines()
    variantFile.close()

    with open(ipvFileName, 'r') as f:
        ipv = json.load(f)
    f.close()

    gen3File = open(gen3FileName)
    gen3 = set(gen3File.readlines())
    gen3File.close()

    variant2samplesDict = dict()

    for variant in variants:
        variant = variant.replace('(', '').replace(')', '')
        variantArray = variant.split(',')
        chrom = int(variantArray[0])
        pos = int(variantArray[1])
        ref = str(variantArray[2])
        alt = str(variantArray[3].replace('\n', ''))
        variantString = str((chrom, pos, ref, alt))

        variant2samplesDict[variant] = dict()
        # now get the samples for the variant from ipv
        print(variantString)
        homoSamples = set(ipv[variantString]['homozygous individuals'])
        if len(homoSamples) == 0:
            variant2samplesDict[variant]['homo'] = None
        elif len(gen3.intersection(homoSamples)) == 0:
            variant2samplesDict[variant]['homo'] = None
        else:
            variant2samplesDict[variant]['homo'] = gen3.intersection(homoSamples)[0]

    print(variant2samplesDict)





if __name__ == "__main__":
    main()