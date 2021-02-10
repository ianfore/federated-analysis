import json
import sys

allFileName=sys.argv[1]
intersectionFileName=sys.argv[2]

with open(allFileName, 'r') as f:
   myAll = json.load(f)
f.close()

with open(intersectionFileName, 'r') as f:
   myIntersection = json.load(f)
f.close()

homoList = list()
for h in myIntersection['homozygous']:
   pos=eval(h)[1]
   homoList.append(pos)

homoCounts=dict()
for v in range(len(myAll['vus'])):
   pos=myAll['vus'][v][0][1]
   zyg=myAll['vus'][v][1]
   if pos in homoList and zyg=='3':
      if not pos in homoCounts:
         homoCounts[pos] = 0
      homoCounts[pos] += 1

print('homoCounts ' + str(homoCounts))

pathCounts=list()
for p in range(len(myAll['pathogenic'])):
   pos=myAll['pathogenic'][p][0][1]
   if not pos in pathCounts:
      pathCounts.append(pos)

print('pathCounts = ' + str(len(pathCounts)))