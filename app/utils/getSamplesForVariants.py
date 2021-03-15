import json
import argparse


def main():
    parser = argparse.ArgumentParser(usage="getSamplesForReads")
    parser.add_argument("--v", dest="variantFileName", help="variant file name", default=None)
    parser.add_argument("--i", dest="ipvFileName", help="ipv file name", default=None)
    parser.add_argument("--g", dest="gen3FileName", help="gen3 file name", default=None)
    parser.add_argument("--o", dest="outputFileName", help="output tsv file name", default=None)
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
    gen3 = gen3File.readlines()
    gen3File.close()
    for i in range(len(gen3)):
        gen3[i] = gen3[i].replace('\n', '')
    gen3 = set(gen3)

    outputFile = open(outputFileName, 'w')
    # NWD825720       chr13:32329240:G:A      chr13:32329140-32329340

    for variant in variants:
        variant = variant.replace('(', '').replace(')', '')
        variantArray = variant.split(',')
        chrom = int(variantArray[0])
        pos = int(variantArray[1])
        ref = str(variantArray[2])
        alt = str(variantArray[3].replace('\n', ''))
        variantString = str((chrom, pos, ref, alt))

        # now get the samples for the variant from ipv
        homoSamples = set(ipv[variantString]['homozygous individuals'])
        if len(homoSamples) == 0:
            pass
        elif len(gen3.intersection(homoSamples)) == 0:
            pass
        else:
            sample = list(gen3.intersection(homoSamples))[0]
            v = 'chr' + str(chrom) + ':' + str(pos) + ':' + ref + ':' + alt
            r = 'chr' + str(chrom) + ':' + str(pos - 100) + '-' + str(pos+100)
            outputFile.write(sample + '\t' + v + '\t' + r + '\n')
            continue

        heteroSamples = set(ipv[variantString]['heterozygous individuals'])
        if len(heteroSamples) == 0:
            continue
        elif len(gen3.intersection(heteroSamples)) == 0:
            continue
        else:
            sample = list(gen3.intersection(heteroSamples))[0]
            v = 'chr' + str(chrom) + ':' + str(pos) + ':' + ref + ':' + alt
            r = 'chr' + str(chrom) + ':' + str(pos - 100) + '-' + str(pos + 100)
            outputFile.write(sample + '\t' + v + '\t' + r + '\n')

    outputFile.close()
    #print(variant2samplesDict)







if __name__ == "__main__":
    main()