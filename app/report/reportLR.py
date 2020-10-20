import json
import sys
import logging

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

def writeHeaders(f):
    header = '#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO'
    print('##fileformat=VCFv4.0', file=f)
    print('##source=brcaexchange', file=f)
    print('##reference=GRCh38', file=f)
    for word in header:
        f.write(word + '\t')
    f.write('\n')

def main():
    if(len(sys.argv) != 3):
        print('usage: input-out.json output-out.json')
        sys.exit(1)

    inputFileName = sys.argv[1]
    outputFileName = sys.argv[2]

    logger.info('reading ' + inputFileName)
    with open(inputFileName, 'r') as f:
        inputDict = json.load(f)
    f.close()

    f = open(outputFileName, 'w')
    writeHeaders(f)

    for vus in inputDict['cooccurring vus']:
        row = [vus[0], vus[1], '.', vus[2], vus[3], '.', '.', '.']
        for word in row:
            f.write(word + '\t')
        f.write('\n')

    f.close()

if __name__ == "__main__":
    main()