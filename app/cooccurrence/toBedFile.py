import sys

def main():
    if len(sys.argv) != 3:
        print('inputDir chrom')
        sys.exit(1)


    coocs, homos, inGnomad, notGnomad = readFiles(sys.argv[1], sys.argv[2])

    bedDict = generateBedDict(coocs, homos, inGnomad, notGnomad)

    print(bedDict)


def generateBedDict(coocs, homos, inGnomad, notGnomad):
    bedDict = {'homo-in':[],
               'homo-not': [],
               'cooc-in':[],
               'cooc-not':[],
               'both-in':[],
               'both-not':[]}
    for c in coocs:
        chrom = c[0]
        pos = c[1]
        length = max(len(c[2]), len(c[3]))
        if c in homos:
            if c in inGnomad:
                bedDict['both-in'].append((chrom, pos, length))
            else:
                bedDict['both-not'].append((chrom, pos, length))
        else:
            if c in inGnomad:
                bedDict['cooc-in'].append((chrom, pos, length))
            else:
                bedDict['cooc-not'].append((chrom, pos, length))

    for h in homos:
        chrom = h[0]
        pos = h[1]
        length = max(len(h[2]), len(h[3]))
        if h in coocs:
            if h in inGnomad:
                bedDict['both-in'].append((chrom, pos, length))
            else:
                bedDict['both-not'].append((chrom, pos, length))
        else:
            if h in inGnomad:
                bedDict['homo-in'].append((chrom, pos, length))
            else:
                bedDict['homo-not'].append((chrom, pos, length))

    return bedDict


def readFiles(inputDir, chrom):
    inputDir = sys.argv[1]
    chrom = sys.argv[2]

    coocsFileName = inputDir + '/' + chrom + '-coocs.txt'
    coocs = list()
    with open(coocsFileName, 'r') as f:
        coocs = [line.rstrip() for line in f]
    f.close()

    homosFileName = inputDir + '/' + chrom + '-homos.txt'
    homos = list()
    with open(homosFileName, 'r') as f:
        homos = [line.rstrip() for line in f]
    f.close()

    inFileName = inputDir + '/' + chrom + '-in.txt'
    inGnomad = list()
    with open(inFileName, 'r') as f:
        inGnomad = [line.rstrip() for line in f]
    f.close()

    notFileName = inputDir + '/' + chrom + '-not.txt'
    notGnomad = list()
    with open(inFileName, 'r') as f:
         notGnomad = [line.rstrip() for line in f]
    f.close()

    return coocs, homos, inGnomad, notGnomad



if __name__ == "__main__":
    main()