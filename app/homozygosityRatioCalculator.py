import sys
import json

def main():
    if len(sys.argv) != 3:
        print('13-vpi.json 13-hrpi.json')
        sys.exit(1)

    vpiFileName = sys.argv[1]
    hrpiFileName = sys.argv[2]

    with open(vpiFileName, 'r') as f:
        vpiDict = json.load(f)
    f.close()

    ratios = getRatios(vpiDict)
    with open(hrpiFileName, 'w') as f:
        json.dump(ratios, f)
    f.close()


def getRatios(vpiDict):
    ratioDict = dict()
    for i in vpiDict:
        homoSum = 0
        hetSum = 0
        for b in vpiDict[i]['benign']:
            if b[1] == "3":
                homoSum += 1
            else:
                hetSum += 1
        for p in vpiDict[i]['pathogenic']:
            if p[1] == "3":
                homoSum += 1
            else:
                hetSum += 1
        for v in vpiDict[i]['vus']:
            if v[1] == "3":
                homoSum += 1
            else:
                hetSum += 1
        if homoSum + hetSum == 0:
            ratioDict[i] = 0
        else:
            ratioDict[i] = float(homoSum) / float(homoSum + hetSum)
    return ratioDict

if __name__ == "__main__":
    main()