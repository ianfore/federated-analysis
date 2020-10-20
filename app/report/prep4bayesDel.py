import json
import logging
import sys
from collections import defaultdict
import numpy as np

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

centersPerHomo = defaultdict(set)
studiesPerHomo = defaultdict(set)
countsPerCenter = dict()
countsPerStudy = dict()

def main():
    # read in vpi
    if len(sys.argv) != 6:
        print('usage: findBatchEffect 13-vpi.json 13-out.json 17-vpi.json 17-out.json output-dir')
        sys.exit(1)
    vpi13FileName = sys.argv[1]
    out13FileName = sys.argv[2]
    vpi17FileName = sys.argv[3]
    out17FileName = sys.argv[4]
    centersPerHomoOutputFileName = sys.argv[5] + '/centersPerHomo.json'
    studiesPerHomoOutputFileName = sys.argv[5] + '/studiesPerHomo.json'
    countsPerCenterOutputFileName = sys.argv[5] + '/countsPerCenter.json'
    countsPerStudyOutputFileName = sys.argv[5] + '/countsPerStudy.json'

    logger.info('reading data from ' + vpi13FileName)
    with open(vpi13FileName, 'r') as f:
        vpi13Dict = json.load(f)
    f.close()

    logger.info('reading data from ' + out13FileName)
    with open(out13FileName, 'r') as f:
        out13Dict = json.load(f)
    f.close()

    logger.info('reading data from ' + vpi17FileName)
    with open(vpi17FileName, 'r') as f:
        vpi17Dict = json.load(f)
    f.close()

    logger.info('reading data from ' + out17FileName)
    with open(out17FileName, 'r') as f:
        out17Dict = json.load(f)
    f.close()

    findBatch(vpi13Dict, out13Dict)
    findBatch(vpi17Dict, out17Dict)

    with open(centersPerHomoOutputFileName, 'w') as f:
        json.dump(centersPerHomo, f, cls=NpEncoder)
    f.close()

    with open(countsPerCenterOutputFileName, 'w') as f:
        json.dump(countsPerCenter, f, cls=NpEncoder)
    f.close()

    with open(studiesPerHomoOutputFileName, 'w') as f:
        json.dump(studiesPerHomo, f, cls=NpEncoder)
    f.close()

    with open(countsPerStudyOutputFileName, 'w') as f:
        json.dump(countsPerStudy, f, cls=NpEncoder)
    f.close()

if __name__ == "__main__":
    main()