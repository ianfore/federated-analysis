import sys
import json

cohortSize = 137977
CUTOFF=0.01

if len(sys.argv) != 2:
   print('need chrom number')
   sys.exit(1)

chrom=sys.argv[1]

pathFileName = chrom + '/' + chrom + '-path.json'
with open(pathFileName, 'r') as f:
   vDict = json.load(f)[0]


print('total co-occurring variants = ' + str(len(vDict)))
for k in vDict.keys():
   n = vDict[k][2]
   freq = float(n) / float(cohortSize)
   if freq < CUTOFF:
      print('variant = ' + str(k) + '(freq = ' + str(freq) + ')')

f.close

homoFileName = chrom + '/' + chrom + '-homo.json'
with open(homoFileName, 'r') as f:
   vDict = json.load(f)[0]

print('total homozygous variants = ' + str(len(vDict)))
for k in vDict.keys():
   freq = vDict[k][1][1]
   if freq < CUTOFF:
      print('variant = ' + str(k) + '(freq = ' + str(freq) + ')')

f.close()