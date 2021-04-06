import json
import sys
from hardyWeinbergTest import hardyWeinberg

def calculatePandQ(vDict, cohortSize):
    for key in vDict:
        aa = vDict[key]['aa']
        Aa = vDict[key]['Aa']
        AA = cohortSize - (aa + Aa)
        vDict[key]['AA'] = AA
        p = (2 * int(AA) + int(Aa)) / (2 * (int(AA) + int(Aa) + int(aa)))
        q = 1 - p
        vDict[key]['p'] = p
        vDict[key]['q'] = q

def initializeGenotypes(vDict, vtype, key):
    if not key in vDict:
        vDict[key] = dict()
        vDict[key]['aa'] = 0
        vDict[key]['Aa'] = 0
        vDict[key]['AA'] = 0
    vDict[key]['class'] = vtype

def doHW(vDict, fileName):
    hw = hardyWeinberg(vDict)
    hw_res = hw.hwTest()
    f = open(fileName, 'w')
    f.write('variant\tclass\taa\tAa\tAA\tp\tq\tchisquare\thweafp\n')
    for key in hw_res:
        f.write(key + '\t' + hw_res[key]['class'] + '\t' +
                str(hw_res[key]['aa']) + '\t' +
                str(hw_res[key]['Aa']) + '\t' +
                str(hw_res[key]['AA']) + '\t' +
                str(hw_res[key]['p']) + '\t' +
                str(hw_res[key]['q']) + '\t' +
                str(hw_res[key]['chisquare']) + '\t' +
                str(hw_res[key]['hail_hweafp']) + '\n')
    f.close()

casesOnlyFile = sys.argv[1]
numCases = int(sys.argv[2])
controlsOnlyFile = sys.argv[3]
numControls = int(sys.argv[4])

with open(casesOnlyFile, 'r') as f:
    cases = json.load(f)
f.close()

with open(controlsOnlyFile, 'r') as f:
    controls = json.load(f)
f.close()

casesAndControlsDict = dict()
casesDict = dict()
controlsDict = dict()

for vtype in ['benign', 'pathogenic', 'vus']:
    # do the cases
    for i in range(len(cases[vtype])):
        key = str(cases[vtype][i][0])
        initializeGenotypes(casesAndControlsDict, vtype, key)
        initializeGenotypes(casesDict, vtype, key)
        if cases[vtype][i][1] == '1':
            casesAndControlsDict[key]['Aa'] += 1
            casesDict[key]['Aa'] += 1
        elif cases[vtype][i][1] == '3':
            casesAndControlsDict[key]['aa'] += 1
            casesDict[key]['aa'] += 1
        else:
            print('shouldnt find homo ref in output')
            sys.exit(1)

    # now do controls
    for j in range(len(controls[vtype])):
        key = str(controls[vtype][j][0])
        initializeGenotypes(casesAndControlsDict, vtype, key)
        initializeGenotypes(controlsDict, vtype, key)

        if controls[vtype][j][1] == '1':
            casesAndControlsDict[key]['Aa'] += 1
            controlsDict[key]['Aa'] += 1
        elif controls[vtype][j][1] == '3':
            casesAndControlsDict[key]['aa'] += 1
            controlsDict[key]['aa'] += 1
        else:
            print('shouldnt find homo ref in output')
            sys.exit(1)

calculatePandQ(casesAndControlsDict, numCases + numControls)
calculatePandQ(casesDict, numCases)
calculatePandQ(controlsDict, numControls)

# do HW for cases and controls
doHW(casesAndControlsDict, 'casesAndControls-hwe.tsv')

# do HW for cases
doHW(casesDict, 'cases-hwe.tsv')

# do HW
doHW(controlsDict, 'controls-hwe.tsv')

