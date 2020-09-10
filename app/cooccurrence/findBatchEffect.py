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

def findBatch(vpiDict):

    for individual in vpiDict:
        for vus in vpiDict[individual]['vus']:
            variant = vus[0]
            genotype = vus[1]
            seqCenter = vus[2]
            study = vus[3]
            if not seqCenter in countsPerCenter:
                countsPerCenter[seqCenter] = {'homoVUS': 0, 'heteroVUS': 0, 'homoBen': 0, 'heteroBen': 0,
                                              'totalHomo': 0, 'totalHetero': 0}
            if not study in countsPerStudy:
                countsPerStudy[study] = {'homoVUS': 0, 'heteroVUS': 0, 'homoBen': 0, 'heteroBen': 0,
                                              'totalHomo': 0, 'totalHetero': 0}
            if genotype == '3':
                centersPerHomo[str(tuple(variant))].add(seqCenter)
                studiesPerHomo[str(tuple(variant))].add(study)
                countsPerCenter[seqCenter]['homoVUS'] += 1
                countsPerCenter[seqCenter]['totalHomo'] += 1
                countsPerStudy[study]['homoVUS'] += 1
                countsPerStudy[study]['totalHomo'] += 1
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
            if not seqCenter in countsPerCenter:
                countsPerCenter[seqCenter] = {'homoVUS': 0, 'heteroVUS': 0, 'homoBen': 0, 'heteroBen': 0,
                                              'totalHomo': 0, 'totalHetero': 0}
            if not study in countsPerStudy:
                countsPerStudy[study] = {'homoVUS': 0, 'heteroVUS': 0, 'homoBen': 0, 'heteroBen': 0,
                                              'totalHomo': 0, 'totalHetero': 0}
            if genotype == '3':
                centersPerHomo[str(tuple(variant))].add(seqCenter)
                studiesPerHomo[str(tuple(variant))].add(study)
                countsPerCenter[seqCenter]['homoBen'] += 1
                countsPerCenter[seqCenter]['totalHomo'] += 1
                countsPerStudy[study]['homoBen'] += 1
                countsPerStudy[study]['totalHomo'] += 1
            else:
                countsPerCenter[seqCenter]['heteroBen'] += 1
                countsPerCenter[seqCenter]['totalHetero'] += 1
                countsPerStudy[study]['heteroBen'] += 1
                countsPerStudy[study]['totalHetero'] += 1

    #return countsPerCenter, countsPerStudy, centersPerHomo, studiesPerHomo

def main():
    # read in vpi
    if len(sys.argv) != 6:
        print('usage: findBatchEffect 13-vpi.json 13-ipv-f.json 17-vpi.json 17-ipv-f.json output-dir')
        sys.exit(1)
    vpi13FileName = sys.argv[1]
    ipv13FileName = sys.argv[2]
    vpi17FileName = sys.argv[3]
    ipv17FileName = sys.argv[4]
    centersPerHomoOutputFileName = sys.argv[5] + '/centersPerHomo.json'
    studiesPerHomoOutputFileName = sys.argv[5] + '/studiesPerHomo.json'
    countsPerCenterOutputFileName = sys.argv[5] + '/countsPerCenter.json'
    countsPerStudyOutputFileName = sys.argv[5] + '/countsPerStudy.json'

    logger.info('reading data from ' + vpi13FileName)
    with open(vpi13FileName, 'r') as f:
        vpi13Dict = json.load(f)
    f.close()

    logger.info('reading data from ' + ipv13FileName)
    with open(ipv13FileName, 'r') as f:
        ipv13Dict = json.load(f)
    f.close()

    logger.info('reading data from ' + vpi17FileName)
    with open(vpi17FileName, 'r') as f:
        vpi17Dict = json.load(f)
    f.close()

    logger.info('reading data from ' + ipv17FileName)
    with open(ipv17FileName, 'r') as f:
        ipv17Dict = json.load(f)
    f.close()

    '''centersPerHomo = defaultdict(set)
    studiesPerHomo = defaultdict(set)
    countsPerCenter = dict()
    countsPerStudy = dict()

    for individual in vpi13Dict:
        for vus in vpi13Dict[individual]['vus']:
            variant = vus[0]
            genotype = vus[1]
            seqCenter = vus[2]
            study = vus[3]
            if not seqCenter in countsPerCenter:
                countsPerCenter[seqCenter] = {'homoVUS': 0, 'heteroVUS': 0, 'homoBen': 0, 'heteroBen': 0,
                                              'totalHomo': 0, 'totalHetero': 0}
            if not study in countsPerStudy:
                countsPerStudy[study] = {'homoVUS': 0, 'heteroVUS': 0, 'homoBen': 0, 'heteroBen': 0,
                                              'totalHomo': 0, 'totalHetero': 0}
            if genotype == '3':
                centersPerHomo[str(tuple(variant))].add(seqCenter)
                studiesPerHomo[str(tuple(variant))].add(study)
                countsPerCenter[seqCenter]['homoVUS'] += 1
                countsPerCenter[seqCenter]['totalHomo'] += 1
                countsPerStudy[study]['homoVUS'] += 1
                countsPerStudy[study]['totalHomo'] += 1
            else:
                countsPerCenter[seqCenter]['heteroVUS'] += 1
                countsPerCenter[seqCenter]['totalHetero'] += 1
                countsPerStudy[study]['heteroVUS'] += 1
                countsPerStudy[study]['totalHetero'] += 1

        for ben in vpi13Dict[individual]['benign']:
            variant = ben[0]
            genotype = ben[1]
            seqCenter = ben[2]
            study = ben[3]
            if not seqCenter in countsPerCenter:
                countsPerCenter[seqCenter] = {'homoVUS': 0, 'heteroVUS': 0, 'homoBen': 0, 'heteroBen': 0,
                                              'totalHomo': 0, 'totalHetero': 0}
            if not study in countsPerStudy:
                countsPerStudy[study] = {'homoVUS': 0, 'heteroVUS': 0, 'homoBen': 0, 'heteroBen': 0,
                                              'totalHomo': 0, 'totalHetero': 0}
            if genotype == '3':
                centersPerHomo[str(tuple(variant))].add(seqCenter)
                studiesPerHomo[str(tuple(variant))].add(study)
                countsPerCenter[seqCenter]['homoBen'] += 1
                countsPerCenter[seqCenter]['totalHomo'] += 1
                countsPerStudy[study]['homoBen'] += 1
                countsPerStudy[study]['totalHomo'] += 1
            else:
                countsPerCenter[seqCenter]['heteroBen'] += 1
                countsPerCenter[seqCenter]['totalHetero'] += 1
                countsPerStudy[study]['heteroBen'] += 1
                countsPerStudy[study]['totalHetero'] += 1


    for individual in vpi17Dict:
        for vus in vpi17Dict[individual]['vus']:
            variant = vus[0]
            genotype = vus[1]
            seqCenter = vus[2]
            study = vus[3]
            if not seqCenter in countsPerCenter:
                countsPerCenter[seqCenter] = {'homoVUS': 0, 'heteroVUS': 0, 'homoBen': 0, 'heteroBen': 0,
                                              'totalHomo': 0, 'totalHetero': 0}
            if not study in countsPerStudy:
                countsPerStudy[study] = {'homoVUS': 0, 'heteroVUS': 0, 'homoBen': 0, 'heteroBen': 0,
                                              'totalHomo': 0, 'totalHetero': 0}
            if genotype == '3':
                centersPerHomo[str(tuple(variant))].add(seqCenter)
                studiesPerHomo[str(tuple(variant))].add(study)
                countsPerCenter[seqCenter]['homoVUS'] += 1
                countsPerCenter[seqCenter]['totalHomo'] += 1
                countsPerStudy[study]['homoVUS'] += 1
                countsPerStudy[study]['totalHomo'] += 1
            else:
                countsPerCenter[seqCenter]['heteroVUS'] += 1
                countsPerCenter[seqCenter]['totalHetero'] += 1
                countsPerStudy[study]['heteroVUS'] += 1
                countsPerStudy[study]['totalHetero'] += 1

        for ben in vpi17Dict[individual]['benign']:
            variant = ben[0]
            genotype = ben[1]
            seqCenter = ben[2]
            study = ben[3]
            if not seqCenter in countsPerCenter:
                countsPerCenter[seqCenter] = {'homoVUS': 0, 'heteroVUS': 0, 'homoBen': 0, 'heteroBen': 0,
                                              'totalHomo': 0, 'totalHetero': 0}
            if not study in countsPerStudy:
                countsPerStudy[study] = {'homoVUS': 0, 'heteroVUS': 0, 'homoBen': 0, 'heteroBen': 0,
                                              'totalHomo': 0, 'totalHetero': 0}
            if genotype == '3':
                centersPerHomo[str(tuple(variant))].add(seqCenter)
                studiesPerHomo[str(tuple(variant))].add(study)
                countsPerCenter[seqCenter]['homoBen'] += 1
                countsPerCenter[seqCenter]['totalHomo'] += 1
                countsPerStudy[study]['homoBen'] += 1
                countsPerStudy[study]['totalHomo'] += 1
            else:
                countsPerCenter[seqCenter]['heteroBen'] += 1
                countsPerCenter[seqCenter]['totalHetero'] += 1
                countsPerStudy[study]['heteroBen'] += 1
                countsPerStudy[study]['totalHetero'] += 1'''

    findBatch(vpi13Dict)
    findBatch(vpi17Dict)


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