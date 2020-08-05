import pandas as pd
import json
import logging
import sys
from collections import defaultdict
import numpy as np

logger = logging.getLogger()
defaultLogLevel = "DEBUG"

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, set):
            return list(obj)
        elif isinstance(obj, tuple):
            return list(obj)
        else:
            return super(NpEncoder, self).default(obj)


def main():
    # read in vpi
    if len(sys.argv) != 3:
        print('usage: findBatchEffect vpi.json output.json')
        sys.exit(1)
    vpiFileName = sys.argv[1]
    outputFileName = sys.argv[2]

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
                centersPerHomoVus[str(tuple(variant))].add(seqCenter)

    with open(outputFileName, 'w') as f:
        json.dump(centersPerHomoVus, f, cls=NpEncoder)
    f.close()

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