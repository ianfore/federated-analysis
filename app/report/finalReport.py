import json
import pandas as pd
import logging
import hail as hl
import argparse
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from utils.coord2hgvs import coordinateMapper

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.WARN)

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--ipv', help='individuals per variant json file')
	parser.add_argument('-v', '--vpi', help='variants per individual json file')
	parser.add_argument('-n', '--no', help='variants not in gnomad file')
	parser.add_argument('-y', '--yes', help='variants in gnomad file')
	parser.add_argument('-s', '--sites', help='variants sites file')
	parser.add_argument('-o', '--output', help='output file')
	parser.add_argument('-m', '--map', help='boolean map hgvs coords')
	return parser.parse_args()

def main():
	ipvFileName = parse_args().ipv
	logger.info('reading data from ' + ipvFileName)
	with open(ipvFileName, 'r') as f:
		ipvDict = json.load(f)
	f.close()

	inFileName = parse_args().yes
	logger.info('reading data from ' + inFileName)
	f = open(inFileName, 'r')
	inList = f.readlines()
	inList = [x.strip() for x in inList]
	f.close()

	outFileName = parse_args().no
	logger.info('reading data from ' + outFileName)
	f = open(outFileName, 'r')
	outList = f.readlines()
	outList = [x.strip() for x in outList]
	f.close()

	sitesFileName = parse_args().sites
	logger.info('reading data from ' + sitesFileName)
	f = open(sitesFileName, 'r')
	sitesDF = pd.read_csv(sitesFileName, header=0, sep='\t')
	f.close()

	vpiFileName = parse_args().vpi
	logger.info('reading data from ' + vpiFileName)
	with open(vpiFileName, 'r') as f:
		vpiDict = json.load(f)
	f.close()

	outputFileName = parse_args().output

	mapHgvs = bool(int(parse_args().map))

	# get batch effect info
	studyPerVariant, centersPerHomoVus = getStudyAndCenter(vpiDict)

	variantsDict = getVariantStats(ipvDict, studyPerVariant, centersPerHomoVus, inList, outList)

	variantsDF = pd.DataFrame.from_dict(variantsDict)
	variantsDF = variantsDF.transpose()
	variantsDF['variant'] = variantsDF.index

	coordMapper = coordinateMapper('GRCh38')
	variantsWithInfoDF = addInfo(variantsDF, sitesDF, coordMapper, mapHgvs)

	logger.info('writing output to ' + outputFileName)
	variantsWithInfoDF.to_csv(outputFileName, sep='\t', index=False)

def getVariantStats(ipvDict, studyPerVariant, centersPerHomoVus, inList, outList):
	allVariants = ipvDict.keys()
	variantsDict = dict()
	for v in allVariants:
		vClass = ipvDict[v]['class']
		vPopFreq = '%.4f'%(ipvDict[v]['maxFreq'])
		vCohortFreq = '%.4f'%(ipvDict[v]['cohortFreq'])
		aa = str(ipvDict[v]['aa'])
		Aa = str(ipvDict[v]['Aa'])
		AA = str(ipvDict[v]['AA'])
		F = str(ipvDict[v]['F'])
		Z = str(ipvDict[v]['Z'])
		p = (2 * int(AA) + int(Aa)) / (2 * (int(AA) + int(Aa) + int(aa)))
		q = 1 - p
		exonic = str(ipvDict[v]['exonic'])
		chisquare = str(ipvDict[v]['chisquare'])
		if len(ipvDict[v]['homozygous individuals']) == 0:
			homoSample = "None"
		else:
			homoSample = ipvDict[v]['homozygous individuals'][0]
		if len(ipvDict[v]['heterozygous individuals']) == 0:
			heteroSample = "None"
		else:
			heteroSample = ipvDict[v]['heterozygous individuals'][0]
		v = v.replace(' ', '')
		v = v.replace("'", "")
		study = studyPerVariant[v]
		if v in inList:
			vIn = 'True'
		elif v in outList:
			vIn = 'False'
		else:
			vIn = 'NA'
		variantsDict[v] = dict()
		variantsDict[v]['class'] = vClass
		variantsDict[v]['popFreq'] = vPopFreq
		variantsDict[v]['cohortFreq'] = vCohortFreq

		variantsDict[v]['homozygousSample'] = homoSample
		variantsDict[v]['heterozygousSample'] = heteroSample
		variantsDict[v]['inGnomad'] = vIn
		variantsDict[v]['homo_alt'] = aa
		variantsDict[v]['hetero'] = Aa
		variantsDict[v]['homo_ref'] = AA
		variantsDict[v]['hail_hweafp'] = hl.eval(hl.hardy_weinberg_test(int(AA),int(Aa),int(aa))).p_value
		variantsDict[v]['F'] = F
		variantsDict[v]['Z'] = Z
		variantsDict[v]['p'] = p
		variantsDict[v]['q'] = q
		variantsDict[v]['chisquare'] = chisquare
		variantsDict[v]['sequenceCenter'] = str(centersPerHomoVus[v]).replace(" ", "")
		variantsDict[v]['exonic'] = exonic
		variantsDict[v]['study'] = study

	return variantsDict

def getStudyAndCenter(vpiDict):
	# get batch effect info
	centersPerHomoVus = dict()
	studyPerVariant = dict()

	for individual in vpiDict:
		for vus in vpiDict[individual]['benign']:
			variant = vus[0]
			seqCenter = vus[2]
			v = str(tuple(variant)).replace("'", "").replace(" ", "")
			if not v in centersPerHomoVus:
				centersPerHomoVus[v] = set()
			centersPerHomoVus[v].add(seqCenter)
			studyPerVariant[v] = vus[3]
		for vus in vpiDict[individual]['pathogenic']:
			variant = vus[0]
			seqCenter = vus[2]
			v = str(tuple(variant)).replace("'", "").replace(" ", "")
			if not v in centersPerHomoVus:
				centersPerHomoVus[v] = set()
			centersPerHomoVus[v].add(seqCenter)
			studyPerVariant[v] = vus[3]
		for vus in vpiDict[individual]['vus']:
			variant = vus[0]
			seqCenter = vus[2]
			v = str(tuple(variant)).replace("'", "").replace(" ", "")
			if not v in centersPerHomoVus:
				centersPerHomoVus[v] = set()
			centersPerHomoVus[v].add(seqCenter)
			studyPerVariant[v] = vus[3]

	return studyPerVariant, centersPerHomoVus

def addInfo(variantsDF, sitesDF, coordMapper, hgvs):
	variants = list(variantsDF['variant'])
	brca_dict = dict()
	pass_dict = dict()
	for index, row in sitesDF.iterrows():
		var = str((str(row['#CHROM']).split('chr')[1], row['POS'], row['REF'], row['ALT']))
		var = var.replace("'", "").replace(" ", "")
		if var in variants:
			brca_dict[var] = row['INFO']
			pass_dict[var] = row['FILTER']

	info_df = pd.DataFrame(brca_dict.items())
	info_df.columns = ['variant', 'INFO']
	interDF = pd.merge(variantsDF, info_df, on='variant', how='left')

	pass_df = pd.DataFrame(pass_dict.items())
	pass_df.columns = ['variant', 'FILTER']
	finalDF = pd.merge(interDF, pass_df, on='variant', how='left')

	# now iterate through the INFO column and pull out each var=val pair
	# we'll make new cols based on these pairs

	infoDict = dict()
	for i in range(len(finalDF.index)):
		if not type(finalDF.iloc[i]['INFO']) is str:
			logger.info('info column N/A: ' + str(finalDF.iloc[i]))
			continue
		try:
			infoPairs = finalDF.iloc[i]['INFO'].split('|')[0].split(';')
		except Exception as e:
			continue
		infoDict[i] = dict()
		for pair in infoPairs:
			logger.debug('pair = ' + pair)
			vv = pair.split('=')
			infoDict[i][vv[0]] = vv[1]
		if hgvs:
			infoDict[i]['hgvs'] = translate_to_hgvs(finalDF.iloc[i]['variant'], coordMapper)

	infoDF = pd.DataFrame.from_dict(infoDict).transpose()

	for k in infoDF.keys():
		finalDF[k] = infoDF[k]

	finalDF = finalDF.drop(columns=['INFO'])

	return finalDF
		

def translate_to_hgvs(vartokens, coordMapper):
	varArray = vartokens.split(',')
	chrom = int(varArray[0].split('(')[1])
	pos = int(varArray[1])
	ref = varArray[2]
	alt = varArray[3].split(')')[0]
	coords = (chrom, pos, ref, alt)
	return coordMapper.translate_to_hgvs(coords)

if __name__ == "__main__":
    main()


