import sys


homozygousInGnomad = 'homo-in'
homozygousNotInGnomad = 'homo-not'
cooccurringInGnomad = 'cooc-in'
cooccurringNotInGnomad = 'cooc-not'
bothInGnomad = 'both-in'
bothNotInGnomad = 'both-not'

def main():
    if len(sys.argv) != 3:
        print('inputDir chrom')
        sys.exit(1)


    coocs, homos, inGnomad, notGnomad = readFiles(sys.argv[1], sys.argv[2])

    bedDict = generateBedDict(coocs, homos, inGnomad, notGnomad)

    writeBedFiles(bedDict, sys.argv[1], sys.argv[2])


def writeBedFiles(bedDict, inputDir, chrom):
    homoInFile = inputDir + '/' + chrom + '-' + homozygousInGnomad
    homoNotInFile = inputDir + '/' + chrom + '-' + homozygousNotInGnomad
    coocInFile = inputDir + '/' + chrom + '-' + cooccurringInGnomad
    coocNotInFile = inputDir + '/' + chrom + '-' + cooccurringNotInGnomad
    bothInFile = inputDir + '/' + chrom + '-' + bothInGnomad
    bothNotInFile = inputDir + '/' + chrom + '-' + bothNotInGnomad

    with open(homoInFile, 'w') as f:
        for item in bedDict[homozygousInGnomad]:
            f.write(str(item))
    f.close()

    with open(homoNotInFile, 'w') as f:
        for item in bedDict[homozygousNotInGnomad]:
            f.write(str(item))
    f.close()

    with open(coocInFile, 'w') as f:
        for item in bedDict[cooccurringInGnomad]:
            f.write(str(item))
    f.close()

    with open(coocNotInFile, 'w') as f:
        for item in bedDict[cooccurringNotInGnomad]:
            f.write(str(item))
    f.close()

    with open(bothInFile, 'w') as f:
        for item in bedDict[bothInGnomad]:
            f.write(str(item))
    f.close()

    with open(bothNotInFile, 'w') as f:
        for item in bedDict[bothNotInGnomad]:
            f.write(str(item))
    f.close()



def generateBedDict(coocs, homos, inGnomad, notGnomad):
    bedDict = {homozygousInGnomad:[],
               homozygousNotInGnomad: [],
               cooccurringInGnomad:[],
               cooccurringNotInGnomad:[],
               bothInGnomad:[],
               bothNotInGnomad:[]}
    for c in coocs:
        chrom = eval(c)[0]
        pos = eval(c)[1]
        length = max(len(eval(c)[2]), len(eval(c)[3]))
        if c in homos:
            if c in inGnomad:
                bedDict[bothInGnomad].append((chrom, pos, length))
            elif c in notGnomad:
                bedDict[bothNotInGnomad].append((chrom, pos, length))
            else:
                print('error: ' + str(c) + ' not in gnomad or not gnomad files!')
                sys.exit(1)
        else:
            if c in inGnomad:
                bedDict[cooccurringInGnomad].append((chrom, pos, length))
            elif c in notGnomad:
                bedDict[cooccurringNotInGnomad].append((chrom, pos, length))
            else:
                print('error: ' + str(c) + ' not in gnomad or not gnomad files!')
                sys.exit(1)

    for h in homos:
        chrom = eval(h)[0]
        pos = eval(h)[1]
        length = max(len(eval(h)[2]), len(eval(h)[3]))
        if h in coocs:
            if h in inGnomad:
                bedDict[bothInGnomad].append((chrom, pos, length))
            elif h in notGnomad:
                bedDict[bothNotInGnomad].append((chrom, pos, length))
            else:
                print('error: ' + str(h) + ' not in gnomad or not gnomad files!')
                sys.exit(1)
        else:
            if h in inGnomad:
                bedDict[homozygousInGnomad].append((chrom, pos, length))
            elif h in notGnomad:
                bedDict[homozygousNotInGnomad].append((chrom, pos, length))
            else:
                print('error: ' + str(h) + ' not in gnomad or not gnomad files!')
                sys.exit(1)

    return bedDict


def readFiles(inputDir, chrom):

    coocsFileName = inputDir + '/' + chrom + '-coocs.txt'
    with open(coocsFileName, 'r') as f:
        coocs = [line.rstrip() for line in f]
    f.close()

    homosFileName = inputDir + '/' + chrom + '-homos.txt'
    with open(homosFileName, 'r') as f:
        homos = [line.rstrip() for line in f]
    f.close()

    inFileName = inputDir + '/' + chrom + '-in.txt'
    with open(inFileName, 'r') as f:
        inGnomad = [line.rstrip() for line in f]
    f.close()

    notFileName = inputDir + '/' + chrom + '-not.txt'
    with open(notFileName, 'r') as f:
         notGnomad = [line.rstrip() for line in f]
    f.close()

    return coocs, homos, inGnomad, notGnomad



if __name__ == "__main__":
    main()