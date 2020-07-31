import allel
import sys
import logging
import json
import pandas as pd

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

def main():

    if len(sys.argv) != 4:
        print('cnv.vcf sample.txt vpi.json')
        sys.exit(1)
    vcfFileName = sys.argv[1]
    sampleFileName = sys.argv[2]
    vpiFileName = sys.argv[3]

    # read cnv.vcf from file into df
    vcf = pd.read_csv(vcfFileName, sep='\t', header = 0)

    # read sample from file into list
    with open(sampleFileName) as f:
        samples = [line.rstrip() for line in f]


    # read vpi.json into dict
    logger.info('reading data from ' + vpiFileName)
    with open(vpiFileName, 'r') as f:
        vpiDict = json.load(f)
    f.close()

    # subset dict samples from cnv list
    vpi = dict()
    for key in samples:
        if key in vpiDict:
            vpi[key] = vpiDict[key]
    #vpi = {key: vpiDict[key] for key in samples}

    potentiallyHemi = list()
    # iterate through samples of vcf
    for i in range(len(vcf)):
        # ref = CNV_chr13_32313762_32316562
        ref = vcf.iloc[i]['ID']
        start = int(ref.split('_')[2])
        stop = int(ref.split('_')[3])
        for sample in vpi:
            cn = int(vcf.iloc[i][sample].split(':')[1])
            for b in vpi[sample]['vus']:
                gt = int(b[1])
                chr = int(b[0][0])
                pos = int(b[0][1])
                if gt == 3 and start <= pos and pos <= stop:
                    potentiallyHemi.append((sample,chr, pos, cn, ref))

        # if cnv is a deletion, then see if sample has that homozygous VUS in that deletion range
    print(potentiallyHemi)
        # add hemizygous vus to list

    # save list of hemizygous to file


if __name__ == "__main__":
    main()