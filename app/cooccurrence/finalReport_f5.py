import json
import sys
import pandas as pd
import logging
import hail as hl

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.INFO)

def main():
	if len(sys.argv) != 7:
		print('ipv-f.json in.txt not.txt sites.tsv vpi.json output.tsv')
		sys.exit(1)


	ipvFileName = sys.argv[1]
	logger.info('reading data from ' + ipvFileName)
	with open(ipvFileName, 'r') as f:
		ipvDict = json.load(f)
	f.close()

	inFileName = sys.argv[2]
	logger.info('reading data from ' + inFileName)
	f = open(inFileName, 'r')
	inList = f.readlines() 
	f.close()
	inList = [x.strip() for x in inList]

	outFileName = sys.argv[3]
	logger.info('reading data from ' + outFileName)
	f = open(outFileName, 'r')
	outList = f.readlines() 
	f.close()
	outList = [x.strip() for x in outList]

	sitesFileName = sys.argv[4]
	logger.info('reading data from ' + sitesFileName)
	f = open(sitesFileName, 'r')
	sitesDF = pd.read_csv(sitesFileName, header=0, sep='\t')
	f.close()

	vpiFileName = sys.argv[5]
	logger.info('reading data from ' + vpiFileName)
	with open(vpiFileName, 'r') as f:
		vpiDict = json.load(f)
	f.close()

	outputFileName = sys.argv[6]

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
			#studyPerVariant[v] = vus[3]
		for vus in vpiDict[individual]['pathogenic']:
			variant = vus[0]
			seqCenter = vus[2]
			v = str(tuple(variant)).replace("'", "").replace(" ", "")
			if not v in centersPerHomoVus:
				centersPerHomoVus[v] = set()
			centersPerHomoVus[v].add(seqCenter)
			#studyPerVariant[v] = vus[3]
		for vus in vpiDict[individual]['vus']:
			variant = vus[0]
			seqCenter = vus[2]
			v = str(tuple(variant)).replace("'", "").replace(" ", "")
			if not v in centersPerHomoVus:
				centersPerHomoVus[v] = set()
			centersPerHomoVus[v].add(seqCenter)
			#studyPerVariant[v] = vus[3]

	allVariants = ipvDict.keys()
	variantsDict = dict()
	#print('variant\tclass\tpopFreq\tcohortFreq\taa\tAa\tAA\thomozygousSample\tinGnomad')
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
			homoSamples = "None"
		else:
			homoSamples = ipvDict[v]['homozygous individuals']
		if len(ipvDict[v]['heterozygous individuals']) == 0:
			heteroSamples = "None"
		else:
			heteroSamples = ipvDict[v]['heterozygous individuals']
		v = v.replace(' ', '')	
		v = v.replace("'", "")
		#study = studyPerVariant[v]
		if v in inList:
			vIn = 'True'
		elif v in outList:
			vIn = 'False'
		else:
			print('neither in in nor out?')
			vIn = 'NA'
		#print(v + '\t' + vClass + '\t' + vPopFreq + '\t' + vCohortFreq + \
		# '\t' + aa + '\t' + Aa + '\t' + AA + '\t' + homoSample + '\t' + vIn)
		variantsDict[v] = dict()
		variantsDict[v]['class'] = vClass
		variantsDict[v]['popFreq'] = vPopFreq
		variantsDict[v]['cohortFreq'] = vCohortFreq
		variantsDict[v]['homozygousSamples'] = homoSamples
		variantsDict[v]['heterozygousSamples'] = heteroSamples
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
		#variantsDict[v]['study'] = study


	variantsDF = pd.DataFrame.from_dict(variantsDict)
	variantsDF = variantsDF.transpose()
	variantsDF['variant'] = variantsDF.index
	#variantsDF.to_csv('/tmp/brcaDF.tsv', sep='\t', index=True)

	variantsWithInfoDF = addInfo(variantsDF, sitesDF)

	logger.info('writing output to ' + outputFileName)
	variantsWithInfoDF.to_csv(outputFileName, sep='\t', index=False)


def addInfo(variantsDF, sitesDF):
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
		for pair in infoPairs:
			print('pair = ' + pair)
			vv = pair.split('=')
			infoDict[i][vv[0]] = vv[1]
	infoDF = pd.DataFrame.from_dict(infoDict).transpose()

	for k in infoDF.keys():
		finalDF[k] = infoDF[k]

	finalDF = finalDF.drop(columns=['INFO'])

	return finalDF
		
	
if __name__ == "__main__":
    main()


