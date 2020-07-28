import sys
import json

def main():
    if len(sys.argv) != 3:
        print('vpi.json output-dir')
        sys.exit(1)
    vpiFileName = sys.argv[1]
    outputDir = sys.argv[2]

    with open(vpiFileName, 'r') as f:
        vpiDict = json.load(f)
    f.close()

    homozygousCount = 0
    heterozygousCount = 0
    for individual in vpiDict:
        for b in vpiDict[individual]['benign']:
            # b = [[13, 32325741, 'C', 'T'], '3']
            print('list(b) = ' + str(list(b)))
            print('list(b)[0] = ' + str(list(b)[0]))
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

        print('individual: ' + str(individual) + ',homo: ' + str(homozygousCount) + ',hetero:' + str(heterozygousCount))
if __name__ == "__main__":
    main()