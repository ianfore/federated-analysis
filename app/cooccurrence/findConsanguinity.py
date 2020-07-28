import sys
import json

def main():
    if len(sys.argv) != 2:
        print('vpi.json')
        sys.exit(1)
    vpiFileName = sys.argv[1]

    with open(vpiFileName, 'r') as f:
        vpiDict = json.load(f)
    f.close()

    for individual in vpiDict:
        homozygousCount = 0
        heterozygousCount = 0
        for b in vpiDict[individual]['benign']:
            # b = [[13, 32325741, 'C', 'T'], '3']
            if list(b)[1] == '3':
                homozygousCount += 1
            else:
                heterozygousCount += 1

        for v in vpiDict[individual]['vus']:
            if list(v)[1] == '3':
                homozygousCount += 1
            else:
                heterozygousCount += 1

        for p in vpiDict[individual]['pathogenic']:
            if list(p)[1] == '3':
                homozygousCount += 1
            else:
                heterozygousCount += 1

        print('individual: ' + str(individual) + ',homo: ' + str(homozygousCount) + ',hetero:' + str(heterozygousCount) +
              ',fraction:' + str(float(homozygousCount)/float(homozygousCount + heterozygousCount)))
if __name__ == "__main__":
    main()