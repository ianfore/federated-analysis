import json
import sys
import pandas as pd
import logging

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

def main():
	if len(sys.argv) != 4:
		print('ipv.json in.txt not.txt')
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
		if len(ipvDict[v]['homozygous individuals']) == 0:
			homoSample = "None"
		else:
			homoSample = ipvDict[v]['homozygous individuals'][0]
		v = v.replace(' ', '')	
		v = v.replace("'", "")
		if v in inList:
			vIn = 'True'
		elif v in outList:
			vIn = 'False'
		else:
			print('neither in in nor out?')
			vIn = 'False'
		#print(v + '\t' + vClass + '\t' + vPopFreq + '\t' + vCohortFreq + \
		# '\t' + aa + '\t' + Aa + '\t' + AA + '\t' + homoSample + '\t' + vIn)
		variantsDict[v] = dict()
		variantsDict[v]['class'] = vClass
		variantsDict[v]['popFreq'] = vPopFreq
		variantsDict[v]['cohortFreq'] = vCohortFreq
		variantsDict[v]['homo_alt'] = aa
		variantsDict[v]['hetero'] = Aa
		variantsDict[v]['homo_ref'] = AA
		variantsDict[v]['homozygousSample'] = homoSample
		variantsDict[v]['inGnomad'] = vIn

	print(variantsDict)


def addInfo():
	brca1DF = pd.read_csv('/Users/jcasaletto/Desktop/TOPMED/f5_copd_brca1_hmb_phased_report.tsv', header=0, sep='\t')
	brca2DF = pd.read_csv('/Users/jcasaletto/Desktop/TOPMED/f5_copd_brca2_hmb_phased_report.tsv', header=0, sep='\t')

	brca1DF_sites = pd.read_csv('/Users/jcasaletto/Desktop/TOPMED/f5_chr17_brca1_sites.tsv', header=0, sep='\t')
	brca2DF_sites = pd.read_csv('/Users/jcasaletto/Desktop/TOPMED/f5_chr13_brca2_sites.tsv', header=0, sep='\t')

	brca1Variants = list(brca1DF['variant'])
	brca1_dict = dict()
	for index, row in brca1DF_sites.iterrows():
		var = str((str(row['#CHROM']).split('chr')[1], row['POS'], row['REF'], row['ALT']))
		var = var.replace("'", "").replace(" ", "")
		if var in brca1Variants:
			brca1_dict[var] = row['INFO']

	brca2Variants = list(brca2DF['variant'])
	brca2_dict = dict()
	for index, row in brca2DF_sites.iterrows():
		var = str((str(row['#CHROM']).split('chr')[1], row['POS'], row['REF'], row['ALT']))
		var = var.replace("'", "").replace(" ", "")
		if var in brca2Variants:
			brca2_dict[var] = row['INFO']

	brca1_info_df = pd.DataFrame(brca1_dict.items())
	brca1_info_df.columns = ['variant', 'INFO']
	brca2_info_df = pd.DataFrame(brca2_dict.items())
	brca2_info_df.columns = ['variant', 'INFO']

	brca1DF_withInfo = pd.merge(brca1DF, brca1_info_df, on='variant', how='left')
	brca2DF_withInfo = pd.merge(brca2DF, brca2_info_df, on='variant', how='left')

	# when we run SQL commands on a DF, SQL will barf if one col named 'aa' and another 'Aa'
	# b/c column names in SQL are case-INsensitive
	# rename the cols ['aa', 'Aa', 'AA'] to ['homo_alt', 'het', 'homo_ref']
	intermediateCols = ['variant', 'class', 'popFreq', 'cohortFreq', 'homo_alt', 'het', 'homo_ref', 'homozygousSample',
						'inGnomad', 'INFO']

	# now iterate through the INFO column and pull out each var=val pair
	# we'll make new cols based on these pairs
	varValPairs = brca1DF_withInfo.iloc[0]['INFO'].split(';')
	infoCols = list()
	for vv in varValPairs:
		infoCols.append(vv.split('=')[0])

	finalCols = intermediateCols + infoCols

	infoDict_1 = dict()
	for i in range(len(brca1DF_withInfo.index)):
		infoDict_1[i] = dict()
		infoPairs = brca1DF_withInfo.iloc[i]['INFO'].split('|')[0].split(';')
		for pair in infoPairs:
			vv = pair.split('=')
			infoDict_1[i][vv[0]] = vv[1]

	# brca1DF_withInfo[vv[0]] = pd.Series(infoDict[i], index=brca1DF_withInfo.index)

	infoDF_1 = pd.DataFrame.from_dict(infoDict_1).transpose()
	print(infoDF_1.columns)
	print(infoDF_1['ABE'])

	brca1DF_withInfo['FIBC_I'] = infoDF_1['FIBC_I']
	brca1DF_withInfo['HWEAF_p'] = infoDF_1['HWEAF_P']
	brca1DF_withInfo['HWE_SLP_I'] = infoDF_1['HWE_SLP_I']
	brca1DF_withInfo['HWE_SLP_P'] = infoDF_1['HWE_SLP_P']
	brca1DF_withInfo['AC'] = infoDF_1['AC']
	brca1DF_withInfo['AF'] = infoDF_1['AF']
	brca1DF_withInfo['AN'] = infoDF_1['AN']

	brca1DF_noInfo = brca1DF_withInfo.drop(columns=['INFO'])

	brca1DF_noInfo.to_csv('/tmp/brca1DF.tsv', sep='\t', index=False)
		
	
if __name__ == "__main__":
    main()


