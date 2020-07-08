import json
import sys
import pandas
import logging

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

def main():
	if len(sys.argv) != 6:
		print('ipv.json homo.txt cooc.txt in.txt not.txt')
		sys.exit(1)


	ipvFileName = sys.argv[1]
	logger.info('reading data from ' + ipvFileName)
	with open(ipvFileName, 'r') as f:
		ipvDict = json.load(f)
	f.close()

	homoFileName = sys.argv[2]
	logger.info('reading data from ' + homoFileName)
	f = open(homoFileName, 'r')
	homoList = f.readlines() 
	f.close()
	homoList = [x.strip() for x in homoList]

	coocFileName = sys.argv[3]
	logger.info('reading data from ' + coocFileName)
	f = open(coocFileName, 'r')
	coocList = f.readlines() 
	f.close()
	coocList = [x.strip() for x in coocList]

	inFileName = sys.argv[4]
	logger.info('reading data from ' + inFileName)
	f = open(inFileName, 'r')
	inList = f.readlines() 
	f.close()
	inList = [x.strip() for x in inList]

	outFileName = sys.argv[5]
	logger.info('reading data from ' + outFileName)
	f = open(outFileName, 'r')
	outList = f.readlines() 
	f.close()
	outList = [x.strip() for x in outList]

	allVariants = ipvDict.keys()
	print('variant\tclass\tpopFreq\tcohortFreq\taa\tAa\tAA\thomozygousSample\tinGnomad')
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
		print(v + '\t' + vClass + '\t' + vPopFreq + '\t' + vCohortFreq + \
			 '\t' + aa + '\t' + Aa + '\t' + AA + '\t' + homoSample + '\t' + vIn) 
		
		
	
if __name__ == "__main__":
    main()


