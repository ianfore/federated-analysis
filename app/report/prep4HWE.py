import json
import sys

casesOnlyFile = sys.argv[1]
numCases = int(sys.argv[2])
controlsOnlyFile = sys.argv[3]
numControls = int(sys.argv[4])

cohortSize = numCases + numControls

with open(casesOnlyFile, 'r') as f:
    cases = json.load(f)
f.close()

with open(controlsOnlyFile, 'r') as f:
    controls = json.load(f)
f.close()

casesAndControls = dict()

for vtype in ['benign', 'pathogenic', 'vus']:
    for i in range(len(cases[vtype])):
        key = str(cases[vtype][i][0])
        if not key in casesAndControls:
            casesAndControls[key] = dict()
            casesAndControls[key]['homo_alt'] = 0
            casesAndControls[key]['hetero'] = 0
            casesAndControls[key]['homo_ref'] = 0
        casesAndControls[key]['class'] = vtype
        if cases[vtype][i][1] == '1':
            casesAndControls[key]['hetero'] += 1
        elif cases[vtype][i][1] == '3':
            casesAndControls[key]['homo_alt'] += 1
        else:
            print('shouldnt find homo ref in output')
            sys.exit(1)
    for j in range(len(controls[vtype])):
        key = str(controls[vtype][j][0])
        if not key in casesAndControls:
            casesAndControls[key] = dict()
            casesAndControls[key]['homo_alt'] = 0
            casesAndControls[key]['hetero'] = 0
            casesAndControls[key]['homo_ref'] = 0
        casesAndControls[key]['class'] = vtype
        if controls[vtype][j][1] == '1':
            casesAndControls[key]['hetero'] += 1
        elif controls[vtype][j][1] == '3':
            casesAndControls[key]['homo_alt'] += 1
        else:
            print('shouldnt find homo ref in output')
            sys.exit(1)

# print header
print('variant\tclass\thomo_alt\thetero\thomo_ref\tp\tq')
for key in casesAndControls:
    aa = casesAndControls[key]['homo_alt']
    Aa = casesAndControls[key]['hetero']
    AA = cohortSize - (aa + Aa)
    casesAndControls[key]['homo_ref'] = AA
    p = (2 * int(AA) + int(Aa)) / (2 * (int(AA) + int(Aa) + int(aa)))
    q = 1 - p
    casesAndControls[key]['p'] = p
    casesAndControls[key]['q'] = q
    print(key + '\t' + casesAndControls[key]['class'] + '\t' +
          str(casesAndControls[key]['homo_alt']) + '\t' +
          str(casesAndControls[key]['hetero']) + '\t' +
          str(casesAndControls[key]['homo_ref']) + '\t' +
          str(casesAndControls[key]['p']) + '\t' +
          str(casesAndControls[key]['q']))