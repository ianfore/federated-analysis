import json
import sys
import logging

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.INFO)

def main():
    if len(sys.argv) != 6:
        print('vpi.json ipv.json rareFreq commonFreq output.json')
        sys.exit(1)

    vpiFileName = sys.argv[1]
    ipvFileName = sys.argv[2]
    rareFreq = float(sys.argv[3])
    commonFreq = float(sys.argv[4])
    outputJson = sys.argv[5]

    logger.info('reading data from ' + vpiFileName)
    with open(vpiFileName, 'r') as f:
        vpiDict = json.load(f)
    f.close()

    logger.info('reading data from ' + ipvFileName)
    with open(ipvFileName, 'r') as f:
        ipvDict = json.load(f)
    f.close()

    contradictoryPerIndividual= dict()
    for individual in vpiDict:
        contradictoryPerIndividual[individual] = dict()
        contradictoryPerIndividual[individual]['common benign het'] = list()
        contradictoryPerIndividual[individual]['rare homo vus'] = list()
        for vus in vpiDict[individual]['benign']:
            variant = tuple(vus[0])
            gt = int(vus[1])
            freq = ipvDict[str(variant)]['maxFreq']
            if freq > commonFreq and (gt == 1 or gt == 2):
                contradictoryPerIndividual[individual]['common benign het'].append(variant)

        for vus in vpiDict[individual]['vus']:
            for vus in vpiDict[individual]['benign']:
                variant = tuple(vus[0])
                gt = int(vus[1])
                freq = ipvDict[str(variant)]['maxFreq']
                if freq < rareFreq and (gt == 3):
                    contradictoryPerIndividual[individual]['rare homo vus'].append(variant)

    # now find intersection of non-empty lists in these 2 dictionaries
    for individual in contradictoryPerIndividual:
        if len(contradictoryPerIndividual[individual]['common benign het']) !=0 and len(contradictoryPerIndividual[individual]['rare homo vus']) !=0:
            print(individual)

if __name__ == "__main__":
    main()