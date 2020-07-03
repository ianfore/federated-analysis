import sys


homozygousInGnomad = 'homo-in'
homozygousNotInGnomad = 'homo-not'
cooccurringInGnomad = 'cooc-in'
cooccurringNotInGnomad = 'cooc-not'
bothInGnomad = 'both-in'
bothNotInGnomad = 'both-not'
inGnomad = 'in'
notInGnomad = 'not'

def main():
    if len(sys.argv) != 3:
        print('inputDir chrom')
        sys.exit(1)


    coocs, homos, inGnomad, notGnomad = readFiles(sys.argv[1], sys.argv[2])

    bedDict = generateBedDict(coocs, homos, inGnomad, notGnomad)

    writeBedFiles(bedDict, sys.argv[1], sys.argv[2])


def writeBedFiles(bedDict, inputDir, chrom):
    homoInFile = inputDir + '/' + chrom + '-' + homozygousInGnomad + '.bed'
    homoNotInFile = inputDir + '/' + chrom + '-' + homozygousNotInGnomad + '.bed'
    coocInFile = inputDir + '/' + chrom + '-' + cooccurringInGnomad + '.bed'
    coocNotInFile = inputDir + '/' + chrom + '-' + cooccurringNotInGnomad + '.bed'
    bothInFile = inputDir + '/' + chrom + '-' + bothInGnomad + '.bed'
    bothNotInFile = inputDir + '/' + chrom + '-' + bothNotInGnomad + '.bed'
    inFile = inputDir + '/' + chrom + '-' + inGnomad + '.bed'
    notFile = inputDir + '/' + chrom + '-' + notInGnomad + '.bed'

    with open(homoInFile, 'w') as f:
        for item in bedDict[homozygousInGnomad]:
            f.write(str(item))
    f.close()

    with open(homoNotInFile, 'w') as f:
        for item in bedDict[homozygousNotInGnomad]:
            f.write(str(item[0]) + '\t' + str(item[1]) + '\t' + str(item[2]) + '\n')
    f.close()

    with open(coocInFile, 'w') as f:
        for item in bedDict[cooccurringInGnomad]:
            f.write(str(item[0]) + '\t' + str(item[1]) + '\t' + str(item[2]) + '\n')
    f.close()

    with open(coocNotInFile, 'w') as f:
        for item in bedDict[cooccurringNotInGnomad]:
            f.write(str(item[0]) + '\t' + str(item[1]) + '\t' + str(item[2]) + '\n')
    f.close()

    with open(bothInFile, 'w') as f:
        for item in bedDict[bothInGnomad]:
            f.write(str(item[0]) + '\t' + str(item[1]) + '\t' + str(item[2]) + '\n')
    f.close()

    with open(bothNotInFile, 'w') as f:
        for item in bedDict[bothNotInGnomad]:
            f.write(str(item[0]) + '\t' + str(item[1]) + '\t' + str(item[2]) + '\n')
    f.close()

    with open(inFile, 'w') as f:
        for item in bedDict[inGnomad]:
            f.write(str(item[0]) + '\t' + str(item[1]) + '\t' + str(item[2]) + '\n')
    f.close()

    with open(notFile, 'w') as f:
        for item in bedDict[notInGnomad]:
            f.write(str(item[0]) + '\t' + str(item[1]) + '\t' + str(item[2]) + '\n')
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
        if not type(chrom) is str or not chrom.startswith('chr'):
            chrom = 'chr' + str(chrom)
        start = eval(c)[1]
        length = max(len(eval(c)[2]), len(eval(c)[3]))
        end = start + length
        if c in homos:
            if c in inGnomad:
                bedDict[bothInGnomad].append((chrom, start, end))
            elif c in notGnomad:
                bedDict[bothNotInGnomad].append((chrom, start, end))
            else:
                print('error: ' + str(c) + ' not in gnomad or not gnomad files!')
                sys.exit(1)
        else:
            if c in inGnomad:
                bedDict[cooccurringInGnomad].append((chrom, start, end))
            elif c in notGnomad:
                bedDict[cooccurringNotInGnomad].append((chrom, start, end))
            else:
                print('error: ' + str(c) + ' not in gnomad or not gnomad files!')
                sys.exit(1)

    for h in homos:
        chrom = eval(h)[0]
        if not type(chrom) is str or not chrom.startswith('chr'):
            chrom = 'chr' + str(chrom)
        start = eval(h)[1]
        length = max(len(eval(h)[2]), len(eval(h)[3]))
        end = start + length
        if h in coocs:
            if h in inGnomad:
                bedDict[bothInGnomad].append((chrom, start, end))
            elif h in notGnomad:
                bedDict[bothNotInGnomad].append((chrom, start, end))
            else:
                print('error: ' + str(h) + ' not in gnomad or not gnomad files!')
                sys.exit(1)
        else:
            if h in inGnomad:
                bedDict[homozygousInGnomad].append((chrom, start, end))
            elif h in notGnomad:
                bedDict[homozygousNotInGnomad].append((chrom, start, end))
            else:
                print('error: ' + str(h) + ' not in gnomad or not gnomad files!')
                sys.exit(1)

    for i in inGnomad:
        chrom = eval(i)[0]
        if not type(chrom) is str or not chrom.startswith('chr'):
            chrom = 'chr' + str(chrom)
        start = eval(i)[1]
        length = max(len(eval(i)[2]), len(eval(i)[3]))
        end = start + length
        bedDict[inGnomad].append((chrom, start, end))

    for n in notInGnomad:
        chrom = eval(n)[0]
        if not type(chrom) is str or not chrom.startswith('chr'):
            chrom = 'chr' + str(chrom)
        start = eval(n)[1]
        length = max(len(eval(n)[2]), len(eval(n)[3]))
        end = start + length
        bedDict[notInGnomad].append((chrom, start, end))


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