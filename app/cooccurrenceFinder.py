import pandas
import itertools
import ast


mutationsFile = '/data/mutationsPerIndividual.txt'
cooccurrencesFile = '/data/cooccurringMutations.txt'
'''fileObject = open(mutationsFile, mode='a')
numMetaDataLines = 8
numColsBeforeSamples = 9
df = pandas.read_csv('/data/BreastCancer.shuffle.vcf', sep='\t', skiprows=numMetaDataLines)


# #CHROM  POS             ID      REF     ALT     QUAL    FILTER  INFO            FORMAT  0000057940      0000057950
# 10      89624243        .       A       G       .       .       AF=1.622e-05    GT      0/0             0/0


# 0/0 => does not have variant on either strand (homozygous negative)
# 0/1  => has variant on 1 strand (heterozygous positive)
# 1/1 =>  has variant on both strands (homozygous positive)


# find mutations
mutationsPerIndvidual = dict()
for index, row in df.iterrows():
    for col in range(numColsBeforeSamples, len(row)):
        if row[col] == '0/1' or row[col] == '1/1':
            if col not in mutationsPerIndvidual:
                mutationsPerIndvidual[col] = list()
            varTuple = (row['#CHROM'], row['POS'], row['REF'], row['ALT'], row['QUAL'])
            mutationsPerIndvidual[col].append(varTuple)

# write out mutations to file
for i in mutationsPerIndvidual:
    print(str(i) + ' = ' + str(mutationsPerIndvidual[i]), file=fileObject)
fileObject.close()'''

# now find co-occurrences

# userID = [(chrom, pos, ref, alt, qual), ... ]

# 10748 = [(10, 89624243, 'A', 'G', '.'), (13, 32910721, 'T', 'C', '.'), (13, 32910842, 'A', 'G', '.')]
# 27089 = [(10, 89624245, 'GA', 'G', '.'), (13, 32906729, 'A', 'C', '.'), (13, 32910721, 'T', 'C', '.')]
# => both have (13, 32910721, 'T', 'C', '.')

# create list of sets, one set for each individual
individuals = list()

# open output file
try:
    fp = open (mutationsFile, 'r')
    record = fp.readline()
    print('read in mutations')
    while record:
        userId, mutations = record.split("=")
        # mutations is a string now, so convert to list
        mutations = ast.literal_eval(mutations.strip())
        user = set()
        # add mutations to user set
        for m in mutations:
            user.add(m)
        individuals.append(user)
        record = fp.readline()

    # find intersections of sets
    print('find intersections')
    cooccurrences = list()
    for a, b in itertools.combinations(individuals, 2):
        intersection = a.intersection(b)
        cooccurrences.append(intersection)

# close file
finally:
    fp.close()

# print list of intersections to stdout
print('write cooccurrences to file')
fileObject = open(mutationsFile, mode='a')
for c in cooccurrences:
    print(str(c), file = fileObject)

fileObject.close()




'''mutationsPerIndividual = dict() # {id: [var1, var2, ...]}

i = 0
vcf_reader = vcf.Reader(open('/data/BreastCancer.shuffle.vcf', 'r'))
for record in vcf_reader:
    for sample in record.samples:
        if sample.is_variant:
            if sample.sample not in mutationsPerIndividual:
                mutationsPerIndividual[sample.sample] = list()
            #mutationsPerIndividual[sample.sample].append(sample.site)
            mutationsPerIndividual[sample.sample].append('1')
            print('length of this list is ' + str(len(mutationsPerIndividual[sample.sample])))
            print('length of map is ' + str(len(mutationsPerIndividual)))

print('okay now we start printing out ')
for individual in mutationsPerIndividual:
    print('mutations for individual ' + individual + ' is ' + str(mutationsPerIndividual[individual]))'''

