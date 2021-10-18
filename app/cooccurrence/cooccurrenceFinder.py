import pandas
import itertools
import time
import sys
import json
from collections import defaultdict
import numpy as np
import os
import argparse
import logging
import allel
from multiprocessing import Process, Queue, cpu_count


logger = logging.getLogger()
defaultLogLevel = "DEBUG"

# you must define the PYENSEMBL_CACHE_DIR before importing the pyensembl module
logger.info('setting pyensembl dir to /var/tmp/pyensembl-cache')
os.environ['PYENSEMBL_CACHE_DIR'] = '/var/tmp/pyensembl-cache'
import pyensembl

# p2 = P(VUS is pathogenic and patient carries a pathogenic variant in trans) (arbitrarily set by goldgar et al)
# Integrated Evaluation of DNA Sequence Variants of Unknown Clinical Significance: Application to BRCA1 and BRCA2
#brca1_p2 = 0.0001
#brca2_p2 = 0.001

# TODO: make these an external file
classStrings = { 'Pathogenic':[ 'Pathogenic',
                                'Likely pathogenic',
                                'Likely_pathogenic',
                                'Pathogenic/Likely_pathogenic'],
                 'Benign':[ 'Benign',
                            'Likely benign',
                            'Likely_benign',
                            'Benign/Likely_benign'],
                 'Unknown': [ 'Uncertain significance',
                              'Uncertain_significance',
                              'Conflicting_interpretations_of_pathogenicity',
                              '-']}
sigColName = 'Clinical_significance'
coordinateColumnBase = 'Genomic_Coordinate_hg'


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, set):
            return list(obj)
        elif isinstance(obj, tuple):
            return list(obj)
        else:
            return super(NpEncoder, self).default(obj)

class NpDecoder(json.JSONDecoder):
    def default(self, obj):
        if isinstance(obj, list):
            return set(obj)
        else:
            return super(NpDecoder, self).default(obj)

def str2bool(string):
    # or just return bool(eval(string))
    '''if string == 'True' or string == 'true':
        return True
    else:
        return False'''
    return bool(eval(string))

def parseArgs():
    parser = argparse.ArgumentParser(usage="cooccurrenceFinder args [options]")
    parser.add_argument("--anno", dest="anno", help="annotation file name, default=None", default=None)
    parser.add_argument("--vcf", dest="vcf", help="vcf file name, default=None", default=None)
    parser.add_argument("--data", dest="data", help="data directory name, default=None", default=None)
    parser.add_argument("--save", dest="save", help="save intermediate files boolean, default=None", default=None)
    parser.add_argument("--h", dest="h", help="Human genome version (37 or 38). Default=None", default=None)
    parser.add_argument("--e", dest="e", help="Ensembl version - 75 (for 37) or 99 (for 38). Default=None", default=None)
    parser.add_argument("--c", dest="c", help="Chromosome of interest. Default=None", default=None)
    parser.add_argument("--g", dest="g", help="Gene of interest. Default=None", default=None)
    parser.add_argument("--p", dest="p", help="Phased (boolean). Default=True", default='True')
    parser.add_argument("--p2", dest="p2", help="Probability VUS is pathogenic and carries pathogenic variant in trans", default=0.001)
    parser.add_argument("--n", dest="n", help="Number of processes. Default=cpu_count()", default=cpu_count())
    parser.add_argument("--vpf", dest="vpf", help="variant pathogenicity  file. Default=None", default=None)
    parser.add_argument("--d", dest="d", help="directory containing pyensembl-cache. Default=/var/tmp/pyensembl-cache",
                        default='/var/tmp/pyensembl-cache')
    parser.add_argument("--spf", dest="spf", help="sample pathology input file. Default=None", default=None)
    parser.add_argument("--gf", dest="gf", help="gnomad sites file name", default=None)
    parser.add_argument("--log", dest="logLevel", help="Logging level. Default=%s" % defaultLogLevel, default=defaultLogLevel)
    return parser.parse_args()

def configureLogging(options):
    # Parse the log level
    numeric_level = getattr(logging, options.logLevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.logLevel)
    # Setup a logger
    logger.setLevel(numeric_level)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(numeric_level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.debug("Established logger")


def main():
    # main just parses the CLI and then calls the run() method with appropriate args
    options = parseArgs()

    configureLogging(options)

    print(options)

    dataDir = options.data
    pathologyFileName = None
    intersectionFile = None
    if dataDir != None:
        outFileName = dataDir + "/" + str(options.g) + "-cooccurrences.json"
        ipvFileName = dataDir + "/" + str(options.g) + "-ipv.json"
        vpiFileName = dataDir + "/" + str(options.g) + "-vpi.json"
        allFileName = dataDir + "/" + str(options.g) + "-all.json"
        toutFileName = dataDir + "/" + str(options.g) + "-tout.json"
        vcfFileName = dataDir + "/" + options.vcf
        pathogenicityFileName = dataDir + "/" + options.vpf
        gnomadFileName = dataDir + "/" + options.gf
        if not options.spf is "":
            pathologyFileName = dataDir + "/" + options.spf
            intersectionFile = dataDir + "/" + str(options.g) + '-intersection.json'
    else:
        outFileName = str(options.g) + "-cooccurrences.json"
        ipvFileName = str(options.g) + "-ipv.json"
        vpiFileName = str(options.g) + "-vpi.json"
        allFileName = str(options.g) + "-all.json"
        toutFileName = str(options.g) + "-tout.json"
        vcfFileName = options.vcf
        pathogenicityFileName =  options.vpf
        gnomadFileName = options.gf
        if not options.spf is None:
            pathologyFileName = options.spf
            intersectionFile = str(options.g) + '-intersection.json'
    saveFiles = str2bool(options.save)
    phased = str2bool(options.p)
    p2 = float(options.p2)


    run(int(options.h), int(options.e), options.c, options.g, phased, p2, vcfFileName,
        int(options.n), pathogenicityFileName, options.d, ipvFileName, vpiFileName, allFileName, options.anno,
        outFileName, toutFileName, saveFiles, pathologyFileName, intersectionFile, gnomadFileName)

def run(hgVersion, ensemblRelease, chromosome, gene, phased, p2, vcfFileName, numProcs,
        pathogenicityFileName, pyensemblDir, ipvFileName, vpiFileName, allVariantsFileName, annoFileName,
        outputFileName, toutFileName, saveFiles, pathologyFileName, intersectionFile, gnomadFileName):


    logger.info('setting pyensembl dir to ' + pyensemblDir)
    os.environ['PYENSEMBL_CACHE_DIR'] = '/var/tmp/pyensembl-cache'

    if not annoFileName is None and annoFileName != '':
        logger.info('reading annotation data from ' + annoFileName)
        with open(annoFileName, 'r') as f:
            annoDF = pandas.read_csv(annoFileName, header=0, sep='\t')
        f.close()
    else:
        annoDF = None

    logger.info('reading data from ' + pathogenicityFileName)
    t = time.time()
    df, pathogenicVariants, benignVariants, unknownVariants = findVariants(pathogenicityFileName, classStrings, hgVersion)
    logger.info('elapsed time in findVariants() ' + str(time.time() -t))

    logger.info('number of pathogenic variants is ' + str(len(pathogenicVariants)))
    logger.info('number of benign variants is ' + str(len(benignVariants)))
    logger.info('number of vus variants is ' + str(len(unknownVariants)))

    if saveFiles:
        myTout = {'benign': benignVariants, 'pathogenic': pathogenicVariants, 'vus': unknownVariants}
        logger.info('saving all variants to ' + toutFileName)
        with open(toutFileName, 'w') as f:
            json.dump(myTout, f, cls=NpEncoder)
        f.close()



    logger.info('reading VCF file ' + vcfFileName)
    t = time.time()
    vcf = readVCFFile(vcfFileName)
    logger.info('elapsed time in readVCFFile() ' + str(time.time() -t))

    logger.info('finding variants per individual in ' + vcfFileName)
    t = time.time()
    q = Queue()
    processList = list()
    for i in range(numProcs):
        p = Process(target=findVarsPerIndividual, args=(q, vcf, benignVariants, pathogenicVariants, chromosome, gene,
                                                        ensemblRelease, annoDF, i, numProcs, ))
        p.start()
        processList.append(p)
    logger.info('joining results from forked threads')
    variantsPerIndividual = dict()
    for i in range(numProcs):
        variantsPerIndividual.update(q.get())
    for i in range(numProcs):
        processList[i].join()
    logger.info('elapsed time in findVariantsPerIndividual() ' + str(time.time() -t))

    t = time.time()
    cohortSize = len(variantsPerIndividual)
    logger.info('number of samples is ' + str(cohortSize))
    individualsPerVariant = findIndividualsPerVariant(variantsPerIndividual, vcf, chromosome,df, hgVersion,
                                                        ensemblRelease, cohortSize, gnomadFileName)
    logger.info('number of records is ' + str(len(individualsPerVariant)))
    logger.info('elapsed time in updateIndividualsPerVariant() ' + str(time.time() -t))

    if saveFiles:
        logger.info('saving vpi to ' + vpiFileName)
        with open(vpiFileName, 'w') as f:
            json.dump(variantsPerIndividual, f, cls=NpEncoder)
        f.close()

        logger.info('saving ipv to ' + ipvFileName)
        with open(ipvFileName, 'w') as f:
            json.dump(individualsPerVariant, f, cls=NpEncoder)
        f.close()

    logger.info('finding homozygous individuals per vus')
    t = time.time()
    homozygousPerVus = countHomozygousPerVus(variantsPerIndividual, df, hgVersion, ensemblRelease, gene, gnomadFileName)
    logger.info('elapsed time in countHomozygousPerVus() ' + str(time.time() -t))

    logger.info('finding individuals per cooc')
    t = time.time()
    individualsPerPathogenicCooccurrence, n, k = findIndividualsPerCooccurrence(variantsPerIndividual, ensemblRelease,
                                                                                phased, gene)
    logger.info('elapsed time in findIndividualsPerCooccurrence() ' + str(time.time() -t))


    # TODO check this math!
    logger.info('calculating p1')
    # p1 = P(VUS is benign and patient carries a path variant in trans) = 0.5 * overall freq of path muts in cohort
    # calculate total number of benign, pathogenic, and vus variants in cohort
    logger.info('getting all variants for cohort')
    allVariants = getAllVariantsPerClass(variantsPerIndividual)
    numPathogenic = len(set(allVariants['pathogenic']))
    p1 =  0.5 * numPathogenic / cohortSize

    if saveFiles:
        logger.info('saving all variants to ' + allVariantsFileName)
        json_dump = json.dumps(allVariants, cls=NpEncoder)
        with open(allVariantsFileName, 'w') as f:
            f.write(json_dump)
        f.close()

    logger.info('putting all the data together per vus')
    dataPerVus = calculateLikelihood(individualsPerPathogenicCooccurrence, p1, p2, n, k, df, hgVersion, cohortSize, gnomadFileName)

    data_set = {"cooccurring vus": dataPerVus, "homozygous vus": homozygousPerVus}
    json_dump = json.dumps(data_set, cls=NpEncoder)

    logger.info('saving final VUS data  to ' + outputFileName)
    with open(outputFileName, 'w') as f:
        f.write(json_dump)
    f.close()

    if not pathologyFileName is None:
        logger.info('intersecting variants with pathology data in file ' + str(pathologyFileName))
        intersectPathology(pathologyFileName, data_set, individualsPerVariant, intersectionFile )

def intersectPathology(pathologyFile, data_set, ipvDF, intersectFile):
    logger.info('reading data from ' + pathologyFile)
    pathologyDF = pandas.read_csv(pathologyFile, sep='\t', header=0)
    fields = pathologyDF.columns

    # determine total number of cases and controls for cohort frequency calculations
    '''numCases = 0
    numSpecialCases = 0
    for i in range(len(pathologyDF)):
        if pandas.isna(pathologyDF.iloc[i]['Age at onset']):
            numSpecialCases += 1
        else:
            numCases += 1
    numTotalCases = numCases + numSpecialCases
    print('numCases = ' + str(numCases))
    print('numSpecialCases = ' + str(numSpecialCases))'''

    numMissing = 0

    pathologyPerCoocIndividual = dict()
    for variant in data_set['cooccurring vus']:
        #for pathogenicVariant in data_set['cooccurring vus'][variant]['pathogenic variants']:
            #pv = str(tuple(pathogenicVariant))
        heterozygousIndividuals = ipvDF[variant]['heterozygous individuals']
        pathologyPerCoocIndividual[variant] = dict()
        pathologyPerCoocIndividual[variant]['phenotype'] = list()
        #pathologyPerCoocIndividual[pv]['numCases'] = 0
        #pathologyPerCoocIndividual[pv]['numSpecialCases'] = 0
        for hi in heterozygousIndividuals:
            pathologies = dict()
            hiInt = int(hi)
            row = pathologyDF.loc[pathologyDF['ID'] == hiInt]
            if len(row) == 0:
                logger.warning('no pathology record for sample ' + hi)
                numMissing += 1

            for field in fields:
                try:
                    pathologies[field] = row[field].tolist()
                except Exception as e:
                    pass

            try:
                pathologyPerCoocIndividual[variant]['phenotype'].append(pathologies)
            except Exception as e:
                pass

    pathologyPerHomoIndividual = dict()
    for variant in data_set['homozygous vus']:
        homozygousIndividuals = ipvDF[variant]['homozygous individuals']
        pathologyPerHomoIndividual[variant] = dict()
        pathologyPerHomoIndividual[variant]['phenotype'] = list()
        #pathologyPerHomoIndividual[variant]['numCases'] = 0
        #pathologyPerHomoIndividual[variant]['numSpecialCases'] = 0
        for hi in homozygousIndividuals:
            pathologies = dict()
            hiInt = int(hi)
            row = pathologyDF.loc[pathologyDF['ID'] == hiInt]
            if len(row) == 0:
                logger.warning('no pathology record for sample ' + hi)
                numMissing += 1
            for field in fields:
                try:
                    pathologies[field] = row[field].tolist()
                except Exception as e:
                    pass

            try:
                pathologyPerHomoIndividual[variant]['phenotype'].append(pathologies)
            except Exception as e:
                pass
                #return False

    pathologyPerAllIndividuals = dict()
    pathologyPerAllIndividuals['homozygous'] = pathologyPerHomoIndividual
    pathologyPerAllIndividuals['cooccurring'] = pathologyPerCoocIndividual
    #pathologyPerAllIndividuals['numCases'] = numCases
    #pathologyPerAllIndividuals['numSpecialCases'] = numSpecialCases
    pathologyPerAllIndividuals['numMissing'] = numMissing

    logger.info('saving pathology intersection data  to ' + intersectFile)
    json_dump = json.dumps(pathologyPerAllIndividuals, indent=4, sort_keys=True)
    with open(intersectFile, 'w') as f:
        f.write(json_dump)
    f.close()

def findIndividualsPerVariant(variantsPerIndividual, vcf, chromosome, df, hgVersion, ensemblRelease, cohortSize, gnomadFileName):
    individualsPerVariant = dict()
    for individual in variantsPerIndividual:
        for b in variantsPerIndividual[individual]['benign']:
            v = str((b[0][0], b[0][1], b[0][2], b[0][3]))
            if not v in individualsPerVariant:
                individualsPerVariant[v] = {'heterozygous individuals': set(),
                                            'homozygous individuals': set()}
            if b[1] == '1' or b[1] == '2':
                individualsPerVariant[v]['heterozygous individuals'].add(individual)
            elif b[1] == '3':
                individualsPerVariant[v]['homozygous individuals'].add(individual)
            else:
                logger.warning('hmm - didnt add this ben ' + v)
        for p in variantsPerIndividual[individual]['pathogenic']:
            v = str((p[0][0], p[0][1], p[0][2], p[0][3]))
            if not v in individualsPerVariant:
                individualsPerVariant[v] = {'heterozygous individuals': set(),
                                            'homozygous individuals': set()}
            if p[1] == '1' or p[1] == '2':
                individualsPerVariant[v]['heterozygous individuals'].add(individual)
            elif p[1] == '3':
                individualsPerVariant[v]['homozygous individuals'].add(individual)
            else:
                logger.warning('hmm - didnt add this path ' + v)
        for vus in variantsPerIndividual[individual]['vus']:
            v = str((vus[0][0], vus[0][1], vus[0][2], vus[0][3]))
            if not v in individualsPerVariant:
                individualsPerVariant[v] = {'heterozygous individuals': set(),
                                            'homozygous individuals': set()}
            if vus[1] == '1' or vus[1] == '2':
                individualsPerVariant[v]['heterozygous individuals'].add(individual)
            elif vus[1] == '3':
                individualsPerVariant[v]['homozygous individuals'].add(individual)
            else:
                logger.warning('hmm - didnt add this vus ' + v)
    individualsPerVariant = addVariantInfo(individualsPerVariant, vcf, chromosome, ['FIBC_I', 'FIBC_P'], df,
                                           hgVersion, cohortSize, ensemblRelease, gnomadFileName)

    return individualsPerVariant


def isExonic(ensemblRelease, chrom, pos):
    ensembl = pyensembl.EnsemblRelease(release=ensemblRelease)
    try:
        exons = ensembl.exons_at_locus(contig=int(chrom), position=int(pos))
    except Exception as e:
        logger.error('exception: ' + str(e))
        return None
    return len(exons) > 0

def addVariantInfo(individualsPerVariant, vcf, chromosome, infoList, df, hgVersion, cohortSize, ensemblRelease, gnomadFileName):
    # add infoList stuff from INFO field
    # do getAFFromGnomadSites inline
    gnomad = pandas.read_csv(gnomadFileName, header=None, sep='\t', comment='#')
    # CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    gnomad.columns = columns
    for variant in range(len(vcf['calldata/GT'])):
        c = int(vcf['variants/CHROM'][variant].replace('chr', ''))
        if c != int(chromosome):
            continue
        p = int(vcf['variants/POS'][variant])
        r = str(vcf['variants/REF'][variant])
        a = str(vcf['variants/ALT'][variant][0])
        v = str((c,p,r,a))
        #if v in individualsPerVariant:
        '''for info in infoList:
            individualsPerVariant[v][info] = vcf['variants/' + info][variant]'''
        #maxPop, maxFreq, minPop, minFreq, allPopFreq = getGnomadData(df, eval(v), hgVersion)
        try:
            individualsPerVariant[v]['cohortFreq'] = float(len(individualsPerVariant[v]['homozygous individuals']) + \
                                                           len(individualsPerVariant[v][
                                                                   'heterozygous individuals'])) / float(cohortSize)
            chrom = 'chr' + str(c)
            info = gnomad[gnomad['#CHROM'] == chrom][gnomad['POS'] == p][gnomad['REF'] == r][gnomad["ALT"] == a][ 'INFO']

            #maxPop, maxFreq, allPopFreq = getAFFromGnomadSites(gnomadFileName, eval(v))
            if len(info) == 1:
                infoArray = info.iloc[0].split(';')
                varValDict = dict()
                for varVal in infoArray:
                    varValArray = varVal.split('=')
                    if len(varValArray) == 1:
                        continue
                    _var = varVal.split('=')[0]
                    _val = varVal.split('=')[1]
                    varValDict[_var] = _val
                popmax = None
                faf95 = None
                af = None
                if 'popmax' in varValDict:
                    popmax = varValDict['popmax']
                if 'faf95_popmax' in varValDict:
                    faf95 = varValDict['faf95_popmax']
                if 'AF' in varValDict:
                    af = varValDict['AF']
            individualsPerVariant[v]['maxPop'] = popmax
            individualsPerVariant[v]['maxFreq'] = af
            individualsPerVariant[v]['faf95'] = faf95
            individualsPerVariant[v]['exonic'] = isExonic(ensemblRelease, c, p)
        except Exception:
            logger.warning('variant ' + str(v) + ' not in ipv dict?')
        #else:
            #logger.warning('variant ' + str(v) + ' not in ipv dict?')
    return individualsPerVariant

def getAllVariantsPerClass(vpi):
    allVariants = dict()
    allVariants['benign'] = list()
    allVariants['pathogenic'] = list()
    allVariants['vus'] = list()

    for i in vpi:
        for b in vpi[i]['benign']:
            allVariants['benign'].append(b)
        for p in vpi[i]['pathogenic']:
            allVariants['pathogenic'].append(p)
        for v in vpi[i]['vus']:
            allVariants['vus'].append(v)

    return allVariants

def findVariants(fileName, classStrings, hgVersion):
    df = pandas.read_csv(fileName, sep='\t', header=0, dtype=str)
    # Genomic_Coordinate_hg37
    # chr13:g.32972575:G>T
    pathVars = set()
    benignVars = set()
    vusVars = set()
    for i in range(len(df)):
        # TODO use HGVS? problem is indel representation
        # TODO VCF has standard of using left-most pos where HGVS has standard of using right-most pos for indel
        # if cDNA (3' side or 5'?) => standard is using 5' strand
        logger.debug(df.loc[i, coordinateColumnBase + str(hgVersion)])
        coord = df.loc[i, coordinateColumnBase + str(hgVersion)].split(':')
        if 'chr' in coord[0]:
            chrom = int(coord[0].split('chr')[1])
        else:
            chrom = int(coord[0])
        if 'g.' in coord[1]:
            pos = int(coord[1].split('g.')[1])
        else:
            pos = int(coord[1])
        ref, alt = coord[2].split('>')
        var = (chrom, pos, ref, alt)
        if str(df.loc[i, sigColName]) in classStrings['Pathogenic']:
            pathVars.add(var)
        elif str(df.loc[i, sigColName]) in classStrings['Benign']:
            benignVars.add(var)
        elif str(df.loc[i, sigColName]) in classStrings['Unknown']:
            vusVars.add(var)
        else:
            continue

    return df, pathVars, benignVars, vusVars

def getAFFromGnomadSites(fileName, vus):
    if not str(vus[0]).startswith("chr"):
        _chrom="chr" + str(vus[0])
    else:
        _chrom=vus[0]
    _pos = vus[1]
    _ref = str(vus[2])
    _alt = str(vus[3])
    gnomad = pandas.read_csv(fileName, header=None, sep='\t', comment='#')
    # CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    gnomad.columns = columns
    info = gnomad[gnomad['#CHROM'] == _chrom][gnomad['POS'] == _pos][gnomad['REF'] == _ref][gnomad["ALT"] == _alt]['INFO']
    '''chromQuery = '#CHROM == ' + "\"" + _chrom + "\""
    posQuery = 'POS == ' + _pos
    refQuery = 'REF == ' + "\"" + _ref + "\""
    altQuery = 'ALT == ' + "\"" + _alt + "\""
    query = "\'" + chromQuery + "\'  and  \'" + posQuery + "\'  and  \'" + refQuery + "\'  and  \'" + altQuery + "\'"
    logger.info(query)
    info = gnomad.query(query)['INFO']'''

    if len(info) == 1:
        infoArray = info.iloc[0].split(';')
        varValDict = dict()
        for varVal in infoArray:
            varValArray = varVal.split('=')
            if len(varValArray) == 1:
                continue
            _var = varVal.split('=')[0]
            _val = varVal.split('=')[1]
            varValDict[_var] = _val
        popmax = None
        faf95 = None
        af = None
        if 'popmax' in varValDict:
            popmax = varValDict['popmax']
        if 'faf95_popmax' in varValDict:
            faf95 = varValDict['faf95_popmax']
        if 'AF' in varValDict:
            af = varValDict['AF']
    else:
        logger.error('malformed VCF file has more than one record for same (chr, pos, ref, alt) tuple')
        return None, None, None
    return (popmax, faf95, af)


def getGnomadData(df, vus, hgVersion):
    # TODO write a unit test
    # 13:g.32393468:C>CT
    #hgString = 'chr' + str(vus[0][0]) + ':g.' + str(vus[0][1]) + ':' + str(vus[0][2]) + '>' + str(vus[0][3])
    hgString = 'chr' + str(vus[0]) + ':g.' + str(vus[1]) + ':' + str(vus[2]) + '>' + str(vus[3])

    # first, get list of columns for GnomAD allleles
    gnomad = [v for v in list(df.columns) if 'GnomAD' in v]
    alleleFrequencies = [v for v in gnomad if 'Allele_frequency' in v]

    # second, get frequencies across exomes and genomes to determine max
    # return population, frequency, count, and number
    # replace "frequency" with "count" and "number" in Allele_frequency_genome_AFR_GnomAD

    maxFrequency = 0.0
    maxPopulation = None
    allPopFreq = dict()
    for af in alleleFrequencies:
        freq=0.0
        alleleFreqList = df[df[coordinateColumnBase + str(hgVersion)] == hgString][af].tolist()
        if alleleFreqList:
            try:
                freq = float(alleleFreqList[0])
            except ValueError:
                continue
            if freq > maxFrequency:
                maxFrequency = freq
                maxPopulation = af
            allPopFreq[af] = freq

    return (maxPopulation, maxFrequency, allPopFreq)

def countHomozygousPerBenign(variantsPerIndividual, df, hgVersion, ensemblRelease, geneOfInterest, gnomadFileName):
    homozygousPerBenign = dict()

    for individual in variantsPerIndividual:
        for ben in variantsPerIndividual[individual]['benign']:
            if (ben[1] == '3') and (getGenesForVariant(ben[0], ensemblRelease, geneOfInterest)):
                if str(ben[0]) not in homozygousPerBenign:
                    homozygousPerBenign[str(ben[0])] = dict()
                    homozygousPerBenign[str(ben[0])]['count'] = 0
                    #maxPop, maxPopFreq, minPop, minPopFreq, allPopFreq = getGnomadData(df, ben[0], hgVersion)
                    maxPop, maxPopFreq, allPopFreq = getAFFromGnomadSites(gnomadFileName, ben[0])
                    homozygousPerBenign[str(ben[0])]['maxPop'] = maxPop
                    homozygousPerBenign[str(ben[0])]['maxPopFreq'] = maxPopFreq
                homozygousPerBenign[str(ben[0])]['count'] += 1

    cohortSize = len(variantsPerIndividual)
    for ben in homozygousPerBenign:
        maxPopFreq = homozygousPerBenign[ben]['maxPopFreq']
        cohortFreq = float(homozygousPerBenign[ben]['count'])/ float(cohortSize)
        homozygousPerBenign[ben]['cohortFreq'] = float(cohortFreq)

    return homozygousPerBenign


def countHomozygousPerVus(variantsPerIndividual, df, hgVersion, ensemblRelease, geneOfInterest, gnomadFileName):
    homozygousPerVus = dict()

    for individual in variantsPerIndividual:
        for vus in variantsPerIndividual[individual]['vus']:
            if (vus[1] == '3') and (getGenesForVariant(vus[0], ensemblRelease, geneOfInterest)):
                if str(vus[0]) not in homozygousPerVus:
                    homozygousPerVus[str(vus[0])] = dict()
                    homozygousPerVus[str(vus[0])]['count'] = 0
                    #maxPop, maxPopFreq, minPop, minPopFreq, allPopFreq = getGnomadData(df, vus[0], hgVersion)
                    maxPop, maxPopFreq, allPopFreq = getAFFromGnomadSites(gnomadFileName, vus[0])
                    homozygousPerVus[str(vus[0])]['maxPop'] = maxPop
                    homozygousPerVus[str(vus[0])]['maxPopFreq'] = maxPopFreq
                homozygousPerVus[str(vus[0])]['count'] += 1

    cohortSize = len(variantsPerIndividual)
    for vus in homozygousPerVus:
        #maxPopFreq = homoZygousPerVus[vus][1][1]
        maxPopFreq = homozygousPerVus[vus]['maxPopFreq']
        cohortFreq = float(homozygousPerVus[vus]['count'])/ float(cohortSize)
        homozygousPerVus[vus]['cohortFreq'] = float(cohortFreq)


    return homozygousPerVus

def calculateLikelihood(pathCoocs, p1, p2, n, k, df, hgVersion, cohortSize, gnomadFileName):

    # vus coocs data: {(vus1, vus2):[individuals]}
    # "([10, 89624243, 'A', 'G'], [10, 89624304, 'C', 'T')]": ["0000057940", "0000057950"],
    # "([10, 89624304, 'C', 'T'], [10, 89624243, 'A', 'G')]": ["0000057940", "0000057950"],
    # "([10, 89624243, 'A', 'G'], [10, 89624298, 'C', 'T')]": ["0000057950"],
    # "([10, 89624304, 'C', 'T'], [10, 89624298, 'C', 'T')]": ["0000057950"],
    # "([10, 89624298, 'C', 'T'], [10, 89624243, 'A', 'G')]": ["0000057950"],
    # "([10, 89624298, 'C', 'T'], [10, 89624304, 'C', 'T')]": ["0000057950"]}

    # new version: {(vus, path):[individuals]}
    # {
    # ((13, 32906579, 'A', 'C'), (13, 32913597, 'CAGAA', 'C')): ['0000061858', '000061999'],
    # ((13, 32906579, 'A', 'C'), (13, 32907420, 'GA', 'G')): ['0999937461']
    # }

    # now calculate log likelihood ratios!
    likelihoodRatios = dict()
    for vus in n:
        n_ = n[vus]
        k_ = k[vus]
        denom = ((p1 ** k_) * (1 - p1) ** (n_ - k_))
        if k_ == 0:
            continue
        elif denom == 0:
            likelihoodRatios[vus] = sys.float_info.min
        else:
            likelihoodRatios[vus] = ((p2 ** k_) * (1 - p2) ** (n_ - k_)) / ((p1 ** k_) * (1 - p1) ** (n_ - k_))


    # find all the pathogenic variants this vus co-occurred with
    pathVarsPerVus = defaultdict(list)
    for cooc in pathCoocs:
        vus = cooc[0]
        pathVarsPerVus[vus].append(cooc[1])

    # put it all together in a single dict
    dataPerVus = dict()
    for vus in likelihoodRatios:
        #maxPop, maxPopFreq, minPop, minPopFreq, allPopFreq = getGnomadData(df, vus, hgVersion)
        maxPop, maxPopFreq, allPopFreq = getAFFromGnomadSites(gnomadFileName, vus)
        cohortFreq = float(n[vus]) / float(cohortSize)
        data = {'likelihood data': {'p1':p1, 'p2':p2, 'n':n[vus], 'k':k[vus], 'likelihood':likelihoodRatios[vus]},
                    'allele frequencies':{'maxPop':maxPop, 'maxPopFreq':maxPopFreq,
                                          'cohortFreq':cohortFreq}, 'pathogenic variants': pathVarsPerVus[vus]}


        dataPerVus[str(vus)] = data

    return dataPerVus

def readVCFFile(vcfFileName):
    return allel.read_vcf(vcfFileName, fields=['samples', 'calldata/GT', 'variants/ALT', 'variants/CHROM',
                                               'variants/FILTER_PASS', 'variants/ID', 'variants/POS', 'variants/QUAL',
                                               'variants/REF', 'variants/INFO'])

def divide(n, d):
   res = list()
   qu = int(n/d)
   rm = n%d
   for i in range(d):
       if i < rm:
           res.append(qu + 1)
       else:
           res.append(qu)
   return res

def getStartAndEnd(partitionSizes, threadID):
    start = 0
    for i in range(threadID):
        start += partitionSizes[i]
    end = start + partitionSizes[threadID]

    return start, end

def findVarsPerIndividual(q, vcf, benignVariants, pathogenicVariants, chromosome, gene, ensemblRelease, annoDF,
                          threadID, numProcesses):
    '''infoFields = ['variants/ABE', 'variants/ABZ', 'variants/AC', 'variants/AF',
     'variants/AN', 'variants/ANN', 'variants/AVGDP', 'variants/BETA_IF', 'variants/BQZ',
     'variants/CYZ', 'variants/FIBC_I', 'variants/FIBC_P', 'variants/FILTER_PASS', 'variants/FLT20', 'variants/GC',
     'variants/GN', 'variants/HWEAF_P', 'variants/HWE_SLP_I', 'variants/HWE_SLP_P', 'variants/IBC_I', 'variants/IBC_P',
     'variants/ID', 'variants/IOR', 'variants/MAX_IF', 'variants/MIN_IF', 'variants/NM0', 'variants/NM1',
     'variants/NMZ', 'variants/NS_NREF', 'variants/QUAL', 'variants/STZ', 'variants/SVM']'''

    variantsPerIndividual = dict()
    individuals = list(vcf['samples'])
    n = len(individuals)
    partitionSizes = divide(n, numProcesses)
    start,end = getStartAndEnd(partitionSizes, threadID)
    logger.info('threadID = ' + str(threadID) + ' processing from ' + str(start) + ' to ' + str(end))
    logger.debug('looping through ' + str(len(individuals)) + ' in samples in VCF')
    for i in range(start, end):
        variantsPerIndividual[individuals[i]] = dict()
        variantsPerIndividual[individuals[i]]['benign'] = list()
        variantsPerIndividual[individuals[i]]['pathogenic'] = list()
        variantsPerIndividual[individuals[i]]['vus'] = list()

        for variant in range(len(vcf['calldata/GT'])):
            if int(vcf['variants/CHROM'][variant].replace('chr', '')) != int(chromosome):
                logger.warning('wrong chromosome?')
                continue
            if 1 in vcf['calldata/GT'][variant][i]:
                c = int(vcf['variants/CHROM'][variant].replace('chr', ''))
                p = int(vcf['variants/POS'][variant])
                r = str(vcf['variants/REF'][variant])
                a = str(vcf['variants/ALT'][variant][0])
                if getGenesForVariant([c,p,r,a], ensemblRelease, gene) is None:
                    logger.warning('no gene for variant? ' + str([c,p,r,a]))
                    continue

                genotype = str(int(str(vcf['calldata/GT'][variant][i][0]) + str(vcf['calldata/GT'][variant][i][1]), 2))
                seqCenter = "NA"
                study = "NA"
                if not annoDF is None:
                    try:
                        if 'CENTER' in individuals[i]:
                            seqCenter = annoDF[annoDF['sample.id'] == individuals[i]]['CENTER'].iloc[0]
                        elif 'seq_center' in individuals[i]:
                            seqCenter = annoDF[annoDF['sample.id'] == individuals[i]]['seq_center'].iloc[0]
                    except Exception as e:
                        seqCenter = "NA"
                    try:
                        study = annoDF[annoDF['sample.id'] == individuals[i]]['study'].iloc[0]
                    except Exception as e:
                        study = "NA"
                if (c, p, r, a) in benignVariants:
                    variantsPerIndividual[individuals[i]]['benign'].append(((c, p, r, a), genotype, seqCenter, study))
                elif (c, p, r, a) in pathogenicVariants:
                    variantsPerIndividual[individuals[i]]['pathogenic'].append(((c, p, r, a), genotype, seqCenter, study))
                # if not a known VUS, it is a VUS now
                else:
                    variantsPerIndividual[individuals[i]]['vus'].append(((c, p, r, a),  genotype, seqCenter, study))

    q.put(variantsPerIndividual)

def getGenesForVariant(variant, ensemblRelease, geneOfInterest):
    ensembl = pyensembl.EnsemblRelease(release=ensemblRelease)
    chrom = variant[0]
    if type(chrom) is str:
        chrom = chrom.split('chr')[1]
    pos = variant[1]
    try:
        #exons = ensembl.exons_at_locus(contig=int(chrom), position=int(pos))
        genes = ensembl.gene_names_at_locus(contig=int(chrom), position=int(pos))
        g_of_i = set()
        g_of_i.add(geneOfInterest)
        g = set(genes)
        intersectingGenes = g_of_i.intersection(g)
        if len(intersectingGenes) != 1:
            return None
        else:
            # this is pythonic way of returning member of singleton set
            return (intersectingGenes,)
    except Exception as e:
        logger.error('exception: ' + str(e))
        return None

def findIndividualsPerCooccurrence(variantsPerIndividual, ensemblRelease, phased, gene):

    individualsPerPathogenicCooccurrence = defaultdict(list)
    n = defaultdict(int)
    k = defaultdict(int)

    for individual in variantsPerIndividual:
        vusVarList = list(variantsPerIndividual[individual]['vus'])
        pathVarList = list(variantsPerIndividual[individual]['pathogenic'])

        for v in vusVarList:
            n[tuple(v[0])] += 1

        vusCrossPath = list(itertools.product(vusVarList, pathVarList))
        for cross in vusCrossPath:
            if sameGeneSameParent(cross[0], cross[1], phased, ensemblRelease, gene):
                k[tuple(cross[0][0])] += 1
                individualsPerPathogenicCooccurrence[(tuple(cross[0][0]), tuple(cross[1][0]))].append(individual)

    return individualsPerPathogenicCooccurrence, n, k

def sameGeneSameParent(vus, path, phased, ensemblRelease, gene):

    if not phased:
        return getGenesForVariant(vus[0], ensemblRelease, gene) == \
               getGenesForVariant(path[0], ensemblRelease, gene)
    else:
        # looking for vus in cis with path
        # if vus is 1|1, then it's either in cis or both in cis and in trans with path
        # else if vus and path are on opposite chromosomes
        return (getGenesForVariant(vus[0], ensemblRelease, gene) == getGenesForVariant(path[0], ensemblRelease, gene)) \
               and \
                ((vus[1] == '3') or ((vus[1] == '2' and path[1] == '1') or (vus[1] == '1' and path[1] == '2')))

                #((vus[1] == '1|1') or ((vus[1] == '1|0' and path[1] == '0|1') or (vus[1] == '0|1' and path[1] == '1|0')))

if __name__ == "__main__":
    main()


