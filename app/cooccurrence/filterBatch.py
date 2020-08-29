import json
import logging
import sys
import numpy as np

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

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
    if len(sys.argv) != 5:
        print('vpi.json batchPerHomoVUS.json batchName filtered-vpi.json')
        sys.exit(1)

    vpiFileName = sys.argv[1]
    perBatchFileName = sys.argv[2]
    batchName = sys.argv[3]
    filteredVPIFileName = sys.argv[4]

    logger.info('reading data from ' + vpiFileName)
    with open(vpiFileName, 'r') as f:
        vpi = json.load(f)
    f.close()

    logger.info('reading data from ' + perBatchFileName)
    with open(perBatchFileName, 'r') as f:
        cphv = json.load(f)
    f.close()

    excludeList = list()
    for vus in cphv:
        if cphv[vus] == [batchName]:
            excludeList.append(vus)

    #for key in excludeList:
    #    del out['homozygous vus'][key]
    for key in excludeList:
        for individual in vpi:
            for v in individual['vus']:
                if v[2] == key:
                    del vpi[individual][v]

    with open(filteredVPIFileName, 'w') as f:
        json.dump(vpi, f, cls=NpEncoder)
    f.close()


if __name__ == "__main__":
    main()