import json
import sys
import pandas as pd
import logging
import hail as hl
logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.INFO)

def main():
	if len(sys.argv) != 5:
		print('ipv-f.json sites.tsv vpi.json output.tsv')
		sys.exit(1)


	ipvFileName = sys.argv[1]
	logger.info('reading data from ' + ipvFileName)
	with open(ipvFileName, 'r') as f:
		ipvDict = json.load(f)
	f.close()

	sitesFileName = sys.argv[2]
	logger.info('reading data from ' + sitesFileName)
	f = open(sitesFileName, 'r')
	sitesDF = pd.read_csv(sitesFileName, header=0, sep='\t')
	f.close()

	vpiFileName = sys.argv[3]
	logger.info('reading data from ' + vpiFileName)
	with open(vpiFileName, 'r') as f:
		vpiDict = json.load(f)
	f.close()

	outputFileName = sys.argv[4]

	# get batch effect info
	centersPerHomoVus = dict()
	for individual in vpiDict:
		for vus in vpiDict[individual]['benign']:
			variant = vus[0]
			seqCenter = vus[2]
			v = str(tuple(variant)).replace("'", "").replace(" ", "")
			if not v in centersPerHomoVus:
				centersPerHomoVus[v] = set()
			centersPerHomoVus[v].add(seqCenter)
		for vus in vpiDict[individual]['pathogenic']:
			variant = vus[0]
			seqCenter = vus[2]
			v = str(tuple(variant)).replace("'", "").replace(" ", "")
			if not v in centersPerHomoVus:
				centersPerHomoVus[v] = set()
			centersPerHomoVus[v].add(seqCenter)
		for vus in vpiDict[individual]['vus']:
			variant = vus[0]
			seqCenter = vus[2]
			v = str(tuple(variant)).replace("'", "").replace(" ", "")
			if not v in centersPerHomoVus:
				centersPerHomoVus[v] = set()
			centersPerHomoVus[v].add(seqCenter)

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
		# calculate p = (2 x Obs(AA) + Obs(Aa)) / (2 x (Obs(AA) + Obs(Aa) + Obs(aa))
		p = (2 * int(AA) + int(Aa)) / (2 * (int(AA) + int(Aa) + int(aa)))
		q = 1 - p
		F = str(ipvDict[v]['F'])
		Z = str(ipvDict[v]['Z'])
		exonic = str(ipvDict[v]['exonic'])
		#fisher = str(ipvDict[v]['fisher'])
		chisquare = str(ipvDict[v]['chisquare'])

		if len(ipvDict[v]['homozygous individuals']) == 0:
			homoSample = "None"
		else:
			homoSample = ipvDict[v]['homozygous individuals'][0]
		v = v.replace(' ', '')	
		v = v.replace("'", "")

		#print(v + '\t' + vClass + '\t' + vPopFreq + '\t' + vCohortFreq + \
		# '\t' + aa + '\t' + Aa + '\t' + AA + '\t' + homoSample + '\t' + vIn)
		variantsDict[v] = dict()
		variantsDict[v]['class'] = vClass
		variantsDict[v]['popFreq'] = vPopFreq
		variantsDict[v]['cohortFreq'] = vCohortFreq
		variantsDict[v]['homozygousSample'] = homoSample
		variantsDict[v]['homo_alt'] = aa
		variantsDict[v]['hetero'] = Aa
		variantsDict[v]['homo_ref'] = AA
		variantsDict[v]['hail_hweafp'] = hl.eval(hl.hardy_weinberg_test(int(AA),int(Aa),int(aa))).p_value
		variantsDict[v]['p'] = p
		variantsDict[v]['q'] = q
		variantsDict[v]['F'] = F
		variantsDict[v]['Z'] = Z
		variantsDict[v]['chisquare'] = chisquare
		variantsDict[v]['sequenceCenter'] = str(centersPerHomoVus[v]).replace(" ", "")
		variantsDict[v]['exonic'] = exonic
		#variantsDict[v]['fisher'] = fisher


	variantsDF = pd.DataFrame.from_dict(variantsDict)
	variantsDF = variantsDF.transpose()
	variantsDF['variant'] = variantsDF.index
	variantsDF.to_csv('/tmp/brcaDF.tsv', sep='\t', index=True)

	variantsWithInfoDF = addInfo(variantsDF, sitesDF)

	logger.info('writing output to ' + outputFileName)
	variantsWithInfoDF.to_csv(outputFileName, sep='\t', index=False)

def addInfo(variantsDF, sitesDF):


	variants = list(variantsDF['variant'])
	brca_dict = dict()
	for index, row in sitesDF.iterrows():
		var = str((str(row['#CHROM']).split('chr')[1], row['POS'], row['REF'], row['ALT']))
		var = var.replace("'", "").replace(" ", "")
		if var in variants:
			brca_dict[var] = row['INFO']

	info_df = pd.DataFrame(brca_dict.items())
	info_df.columns = ['variant', 'INFO']


	finalDF = pd.merge(variantsDF, info_df, on='variant', how='left')


	# now iterate through the INFO column and pull out each var=val pair
	# we'll make new cols based on these pairs

	infoDict = dict()
	for i in range(len(finalDF.index)):
		infoDict[i] = dict()

		infoPairs = finalDF.iloc[i]['INFO'].split('|')[0].split(';')
		for pair in infoPairs:
			vv = pair.split('=')
			infoDict[i][vv[0]] = vv[1]


	infoDF = pd.DataFrame.from_dict(infoDict).transpose()

	finalDF['FIBC_I'] = infoDF['FIBC_I']
	finalDF['HWEAF_P'] = infoDF['HWEAF_P']
	finalDF['HWE_SLP_I'] = infoDF['HWE_SLP_I']
	finalDF['HWE_SLP_P'] = infoDF['HWE_SLP_P']
	finalDF['AC'] = infoDF['AC']
	finalDF['AF'] = infoDF['AF']
	finalDF['AN'] = infoDF['AN']

	finalDF = finalDF.drop(columns=['INFO'])

	return finalDF
		
	
if __name__ == "__main__":
    main()


