import sys
import json



def main():
   if len(sys.argv) != 2:
      print('provide path to json file as arg')
      sys.exit(1)

   fileName=sys.argv[1]

   with open(fileName, 'r') as f:
      vpiDict = json.load(f)

   genotypeCounts = {'benign': {'homo':0, 'hetero': 0},
                     'pathogenic': {'homo': 0, 'hetero': 0},
                     'vus': {'homo': 0, 'hetero': 0}}
   for individual in vpiDict:
      for b in vpiDict[individual]['benign']:
         if b:
            if b[1] == '3':
               genotypeCounts['benign']['homo'] += 1
            else:
               genotypeCounts['benign']['hetero'] += 1
      for p in vpiDict[individual]['pathogenic']:
         if p:
            if p[1] == '3':
               genotypeCounts['pathogenic']['homo'] += 1
            else:
               genotypeCounts['pathogenic']['hetero'] += 1
      for v in vpiDict[individual]['vus']:
         if v:
            if v[1] == '3':
               genotypeCounts['vus']['homo'] += 1
            else:
               genotypeCounts['vus']['hetero'] += 1

   print(genotypeCounts)

if __name__ == "__main__":
    main()