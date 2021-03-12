import json
import logging
import argparse

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

def writeHeader(f):
    header = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',  'INFO']
    for word in header:
        f.write(word + '\t')
    f.write('\n')

def main():
    parser = argparse.ArgumentParser(usage="prep4bayesDel --i input-json --o output-vcf")
    parser.add_argument("--i", dest="input", help="input json file", default=None)
    parser.add_argument("--o", dest="output", help="output vcf file", default=None)
    options = parser.parse_args()
    inputFileName = options.input
    outputFileName = options.output
    
    logger.info('reading ' + inputFileName)
    with open(inputFileName, 'r') as f:
        inputDict = json.load(f)
    f.close()

    # sort input on chr then pos
    coocList = list(inputDict['cooccurring vus'].keys())
    homoList = list(inputDict['homozygous vus'].keys())
    inputList = coocList + homoList
    inputList = list(set(inputList))
    inputList.sort()
    f = open(outputFileName, 'w')
    writeHeader(f)

    for vus in inputList:
        vus = eval(vus)
        row = [vus[0], vus[1], '.', vus[2], vus[3], '.', '.', '.']
        for word in row:
            f.write(str(word) + '\t')
        f.write('\n')

    f.close()

if __name__ == "__main__":
    main()