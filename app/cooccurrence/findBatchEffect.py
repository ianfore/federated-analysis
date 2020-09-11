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

def findBatch(vpiDict, outDict):

    for individual in vpiDict:
        for vus in vpiDict[individual]['vus']:
            variant = vus[0]
            genotype = vus[1]
            seqCenter = vus[2]
            study = vus[3]
            varStr = str(tuple(variant))
            if not seqCenter in countsPerCenter:
                countsPerCenter[seqCenter] = {'homoVUS': 0, 'heteroVUS': 0,
                                              'homoBen': 0, 'heteroBen': 0,
                                              'totalHomo': 0, 'totalHetero': 0,
                                              'homoVUS_0': 0,
                                              'homoVUS_0.1': 0, 'homoVUS_0.01': 0,
                                              'homoVUS_0.001': 0, 'homoVUS_0.0001': 0,
                                              'homoVUS_0.00001':0
                                              }
            if not study in countsPerStudy:
                countsPerStudy[study] = {'homoVUS': 0, 'heteroVUS': 0,
                                         'homoBen': 0, 'heteroBen': 0,
                                         'totalHomo': 0, 'totalHetero': 0,
                                         'homoVUS_0': 0,
                                         'homoVUS_0.1': 0, 'homoVUS_0.01': 0,
                                         'homoVUS_0.001': 0, 'homoVUS_0.0001': 0,
                                         'homoVUS_0.00001': 0
                                         }
            if genotype == '3':
                centersPerHomo[varStr].add(seqCenter)
                studiesPerHomo[varStr].add(study)
                countsPerCenter[seqCenter]['homoVUS'] += 1
                countsPerCenter[seqCenter]['totalHomo'] += 1
                countsPerStudy[study]['homoVUS'] += 1
                countsPerStudy[study]['totalHomo'] += 1
                if varStr in outDict['homozygous vus']:
                    freq = 0.5 * (outDict['homozygous vus'][varStr]['maxPopFreq'] + outDict['homozygous vus'][varStr]['cohortFreq'])
                    if freq <= 0.00001:
                        countsPerCenter[seqCenter]['homoVUS_0.00001'] += 1
                    elif freq <= 0.0001:
                        countsPerCenter[seqCenter]['homoVUS_0.0001'] += 1
                    elif freq <= 0.001:
                        countsPerCenter[seqCenter]['homoVUS_0.001'] += 1
                    elif freq < 0.01:
                        countsPerCenter[seqCenter]['homoVUS_0.01'] += 1
                    elif freq < 0.1:
                        countsPerCenter[seqCenter]['homoVUS_0.1'] += 1
                    else:
                        countsPerCenter[seqCenter]['homoVUS_0'] += 1
            else:
                countsPerCenter[seqCenter]['heteroVUS'] += 1
                countsPerCenter[seqCenter]['totalHetero'] += 1
                countsPerStudy[study]['heteroVUS'] += 1
                countsPerStudy[study]['totalHetero'] += 1


        for ben in vpiDict[individual]['benign']:
            variant = ben[0]
            genotype = ben[1]
            seqCenter = ben[2]
            study = ben[3]
            varStr = str(tuple(variant))
            if not seqCenter in countsPerCenter:
                countsPerCenter[seqCenter] = {'homoVUS': 0, 'heteroVUS': 0,
                                              'homoBen': 0, 'heteroBen': 0,
                                              'totalHomo': 0, 'totalHetero': 0,
                                              'homoVUS_0': 0,
                                              'homoVUS_0.1': 0, 'homoVUS_0.01': 0,
                                              'homoVUS_0.001': 0, 'homoVUS_0.0001': 0
                                              }
            if not study in countsPerStudy:
                countsPerStudy[study] = {'homoVUS': 0, 'heteroVUS': 0,
                                         'homoBen': 0, 'heteroBen': 0,
                                         'homoVUS_0': 0,
                                         'totalHomo': 0, 'totalHetero': 0,
                                         'homoVUS_0.1': 0, 'homoVUS_0.01': 0,
                                         'homoVUS_0.001': 0, 'homoVUS_0.0001': 0
                                         }
            if genotype == '3':
                centersPerHomo[varStr].add(seqCenter)
                studiesPerHomo[varStr].add(study)
                countsPerCenter[seqCenter]['homoBen'] += 1
                countsPerCenter[seqCenter]['totalHomo'] += 1
                countsPerStudy[study]['homoBen'] += 1
                countsPerStudy[study]['totalHomo'] += 1
            else:
                countsPerCenter[seqCenter]['heteroBen'] += 1
                countsPerCenter[seqCenter]['totalHetero'] += 1
                countsPerStudy[study]['heteroBen'] += 1
                countsPerStudy[study]['totalHetero'] += 1

    statList = ['homoVUS', 'heteroVUS', 'homoBen', 'heteroBen', 'homoVUS_0', 'totalHomo', 'totalHetero',
                'homoVUS_0.1', 'homoVUS_0.01', 'homoVUS_0.001', 'homoVUS_0.0001']
    for center in countsPerCenter:
        total = countsPerCenter[center]['totalHomo'] + countsPerCenter[center]['totalHetero']
        for stat in statList:
            countsPerCenter[center][stat + '_ratio'] = countsPerCenter[center][stat] / total
    for study in countsPerStudy:
        total = countsPerStudy[study]['totalHomo'] + countsPerStudy[study]['totalHetero']
        for stat in statList:
            countsPerStudy[study][stat + '_ratio'] = countsPerStudy[study][stat] / total


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