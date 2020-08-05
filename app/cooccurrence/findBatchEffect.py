import pandas as pd
import json
import logging
import sys
from collections import defaultdict

logger = logging.getLogger()
defaultLogLevel = "DEBUG"

def main():
    # read in vpi
    if len(sys.argv) != 2:
        print('vpi.json')
        sys.exit(1)
    vpiFileName = sys.argv[1]
    logger.info('reading data from ' + vpiFileName)
    with open(vpiFileName, 'r') as f:
        vpiDict = json.load(f)
    f.close()

    centersPerHomoVus = defaultdict(set)
    for individual in vpiDict:
        for vus in vpiDict[individual]['vus']:
            variant = vus[0]
            genotype = vus[1]
            seqCenter = vus[2]
            if genotype == '3':
                centersPerHomoVus[str(variant)].add(seqCenter)
    print(centersPerHomoVus)


if __name__ == "__main__":
    main()

'''
{
    "NWD977487": {
        "benign": [
            [
                [
                    13,
                    32315655,
                    "A",
                    "G"
                ],
                "2",
                "BAYLOR"
            ],
            [
                [
                    13,
                    32315831,
                    "G",
                    "A"
                ],
                "1",
                "BAYLOR"
            ],

'''