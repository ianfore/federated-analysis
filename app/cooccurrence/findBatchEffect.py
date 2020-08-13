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
    if len(sys.argv) != 4:
        print('usage: findBatchEffect 13-vpi.json 17-vpi.json output-dir')
        sys.exit(1)
    vpi13FileName = sys.argv[1]
    vpi17FileName = sys.argv[2]
    centersPerHomoOutputFileName = sys.argv[3] + '/centersPerHomo.json'
    countsPerCenterOutputFileName = sys.argv[3] + '/countsPerCenter.json'

    logger.info('reading data from ' + vpi13FileName)
    with open(vpi13FileName, 'r') as f:
        vpi13Dict = json.load(f)
    f.close()

    logger.info('reading data from ' + vpi17FileName)
    with open(vpi17FileName, 'r') as f:
        vpi17Dict = json.load(f)
    f.close()

    centersPerHomo = defaultdict(set)
    countsPerCenter = dict()

    for individual in vpi13Dict:
        for vus in vpi13Dict[individual]['vus']:
            variant = vus[0]
            genotype = vus[1]
            seqCenter = vus[2]
            if not seqCenter in countsPerCenter:
                countsPerCenter[seqCenter] = {'homoVUS': 0, 'heteroVUS': 0}
            if genotype == '3':
                centersPerHomo[str(tuple(variant))].add(seqCenter)
                countsPerCenter[seqCenter]['homoVUS'] += 1
            else:
                countsPerCenter[seqCenter]['heteroVUS'] += 1

        for ben in vpi13Dict[individual]['benign']:
            variant = ben[0]
            genotype = ben[1]
            seqCenter = ben[2]
            if not seqCenter in countsPerCenter:
                countsPerCenter[seqCenter] = {'homoBen': 0, 'heteroBen': 0}
            if genotype == '3':
                centersPerHomo[str(tuple(variant))].add(seqCenter)
                countsPerCenter[seqCenter]['homoBen'] += 1
            else:
                countsPerCenter[seqCenter]['heteroBen'] += 1


    for individual in vpi17Dict:
        for vus in vpi17Dict[individual]['vus']:
            variant = vus[0]
            genotype = vus[1]
            seqCenter = vus[2]
            if not seqCenter in countsPerCenter:
                countsPerCenter[seqCenter] = {'homoVUS': 0, 'heteroVUS': 0}
            if genotype == '3':
                centersPerHomo[str(tuple(variant))].add(seqCenter)
                countsPerCenter[seqCenter]['homoVUS'] += 1
            else:
                countsPerCenter[seqCenter]['heteroVUS'] += 1

        for ben in vpi17Dict[individual]['benign']:
            variant = ben[0]
            genotype = ben[1]
            seqCenter = ben[2]
            if not seqCenter in countsPerCenter:
                countsPerCenter[seqCenter] = {'homoBen': 0, 'heteroBen': 0}
            if genotype == '3':
                centersPerHomo[str(tuple(variant))].add(seqCenter)
                countsPerCenter[seqCenter]['homoBen'] += 1
            else:
                countsPerCenter[seqCenter]['heteroBen'] += 1

    with open(centersPerHomoOutputFileName, 'w') as f:
        json.dump(centersPerHomo, f, cls=NpEncoder)
    f.close()

    with open(countsPerCenterOutputFileName, 'w') as f:
        json.dump(countsPerCenter, f, cls=NpEncoder)
    f.close()

if __name__ == "__main__":
    main()