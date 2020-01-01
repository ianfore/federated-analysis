import vcf

vcf_reader = vcf.Reader(open('data/BreastCancer.shuffle-test.vcf.gz', 'r'))
mutationsPerIndividual = dict()
for record in vcf_reader:
    for sample in record.samples:
        if sample.is_variant:
            if sample.sample not in mutationsPerIndividual:
                mutationsPerIndividual[sample.sample] = list()
            mutationsPerIndividual[sample.sample].append(sample.site)

print('okay now we start printing out ')
for individual in mutationsPerIndividual:
    print('mutations for individual ' + individual + ' is ' + str(mutationsPerIndividual[individual]))


