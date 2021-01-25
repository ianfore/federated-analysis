import pandas as pd
import argparse
import json

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outputDir', help='output directory')
    parser.add_argument('-i', '--ipvFile', help='input ipv json file')
    options = parser.parse_args()
    return options


def main():
    outputDir = parse_args().outputDir
    ipvFile = parse_args().ipvFile

    with open(ipvFile, 'r') as f:
        ipvData = json.load(f)
    f.close()

    ipvDF = pd.DataFrame.from_dict(ipvData)
    onlyHomoList = list()
    onlyHeteroList = list()
    bothList = list()

    for variant in ipvDF:
        numHom = len(ipvDF[variant]['homozygous individuals'])
        numHet = len(ipvDF[variant]['heterozygous individuals'])
        if numHom == 0 and numHet != 0:
            onlyHeteroList.append(variant)
        elif numHom !=0 and numHet == 0:
            onlyHomoList.append(variant)
        elif numHom !=0 and numHet != 0:
            bothList.append(variant)
        else:
            print('wtf? ' + str(variant))

    with open(outputDir + '/only-homo.txt', 'w') as f:
        for item in onlyHomoList:
            f.write("%s\n" % item)
    f.close()

    with open(outputDir + '/only-het.txt', 'w') as f:
        for item in onlyHeteroList:
            f.write("%s\n" % item)
    f.close()

    with open(outputDir + '/both.txt', 'w') as f:
        for item in bothList:
            f.write("%s\n" % item)
    f.close()


if __name__ == "__main__":
    main()