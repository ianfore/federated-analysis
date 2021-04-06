import json
import sys

casesOnlyFile = sys.argv[1]
controlsOnlyFile = sys.argv[2]


with open(casesOnlyFile, 'r') as f:
    cases = json.load(f)
f.close()

with open(controlsOnlyFile, 'r') as f:
    controls = json.load(f)
f.close()

numCases = len(cases)
numControls = len(controls)
cohortSize = numCases + numControls

casesAndControlsDict = dict()
casesDict = dict()
controlsDict = dict()

for vtype in ['benign', 'pathogenic', 'vus']:
    # do the cases
    for i in range(len(cases[vtype])):
        key = str(cases[vtype][i][0])
        if not key in casesAndControlsDict:
            casesAndControlsDict[key] = dict()
            casesAndControlsDict[key]['aa'] = 0
            casesAndControlsDict[key]['Aa'] = 0
            casesAndControlsDict[key]['AA'] = 0
        casesAndControlsDict[key]['class'] = vtype

        if not key in casesDict:
            casesDict[key] = dict()
            casesDict[key]['aa'] = 0
            casesDict[key]['Aa'] = 0
            casesDict[key]['AA'] = 0
        casesDict[key]['class'] = vtype

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
        if not key in casesAndControlsDict:
            casesAndControlsDict[key] = dict()
            casesAndControlsDict[key]['aa'] = 0
            casesAndControlsDict[key]['Aa'] = 0
            casesAndControlsDict[key]['AA'] = 0
        casesAndControlsDict[key]['class'] = vtype

        if not key in controlsDict:
            controlsDict[key] = dict()
            controlsDict[key]['aa'] = 0
            controlsDict[key]['Aa'] = 0
            controlsDict[key]['AA'] = 0
        controlsDict[key]['class'] = vtype

        if controls[vtype][j][1] == '1':
            casesAndControlsDict[key]['Aa'] += 1
            controlsDict[key]['Aa'] += 1
        elif controls[vtype][j][1] == '3':
            casesAndControlsDict[key]['aa'] += 1
            controlsDict[key]['aa'] += 1
        else:
            print('shouldnt find homo ref in output')
            sys.exit(1)

# write casesAndControls to file
casesAndControlsFile = open('casesAndControls-hwe.tsv', 'w')
# print header
casesAndControlsFile.write('variant\tclass\thomo_alt\thetero\thomo_ref\tp\tq\n')
for key in casesAndControlsDict:
    aa = casesAndControlsDict[key]['aa']
    Aa = casesAndControlsDict[key]['Aa']
    AA = len(casesAndControlsDict) - (aa + Aa)
    casesAndControlsDict[key]['AA'] = AA
    p = (2 * int(AA) + int(Aa)) / (2 * (int(AA) + int(Aa) + int(aa)))
    q = 1 - p
    casesAndControlsDict[key]['p'] = p
    casesAndControlsDict[key]['q'] = q
    casesAndControlsFile.write(key + '\t' + casesAndControlsDict[key]['class'] + '\t' +
          str(casesAndControlsDict[key]['aa']) + '\t' +
          str(casesAndControlsDict[key]['Aa']) + '\t' +
          str(casesAndControlsDict[key]['AA']) + '\t' +
          str(casesAndControlsDict[key]['p']) + '\t' +
          str(casesAndControlsDict[key]['q']) + '\n')
casesAndControlsFile.close()

# write cases to file
# print header
casesFile = open('cases-hwe.tsv', 'w')
casesFile.write('variant\tclass\thomo_alt\thetero\thomo_ref\tp\tq\n')
for key in casesDict:
    aa = casesDict[key]['aa']
    Aa = casesDict[key]['Aa']
    AA = len(casesDict) - (aa + Aa)
    casesDict[key]['AA'] = AA
    p = (2 * int(AA) + int(Aa)) / (2 * (int(AA) + int(Aa) + int(aa)))
    q = 1 - p
    casesDict[key]['p'] = p
    casesDict[key]['q'] = q
    casesFile.write(key + '\t' + casesDict[key]['class'] + '\t' +
          str(casesDict[key]['aa']) + '\t' +
          str(casesDict[key]['Aa']) + '\t' +
          str(casesDict[key]['AA']) + '\t' +
          str(casesDict[key]['p']) + '\t' +
          str(casesDict[key]['q']) + '\n')
casesFile.close()

# write controls to file
controlsFile = open('controls-hwe.tsv', 'w')
# print header
controlsFile.write('variant\tclass\thomo_alt\thetero\thomo_ref\tp\tq\n')
for key in controlsDict:
    aa = controlsDict[key]['aa']
    Aa = controlsDict[key]['Aa']
    AA = len(controlsDict) - (aa + Aa)
    controlsDict[key]['AA'] = AA
    p = (2 * int(AA) + int(Aa)) / (2 * (int(AA) + int(Aa) + int(aa)))
    q = 1 - p
    controlsDict[key]['p'] = p
    controlsDict[key]['q'] = q
    controlsFile.write(key + '\t' + controlsDict[key]['class'] + '\t' +
          str(controlsDict[key]['aa']) + '\t' +
          str(controlsDict[key]['Aa']) + '\t' +
          str(controlsDict[key]['AA']) + '\t' +
          str(controlsDict[key]['p']) + '\t' +
          str(controlsDict[key]['q']) + '\n')
controlsFile.close()