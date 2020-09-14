import json
import pandas as pd
import sys

def main():
    if len(sys.argv) != 4:
        print(
            'chr-out.json chr-ipv.json shuffle.tsv ')
        sys.exit(1)
    variantsFile = sys.argv[1]
    ipvFile = sys.argv[2]
    pathologyFile = sys.argv[3]

    with open(variantsFile, 'r') as f:
        variantsDF = json.load(f)
    f.close()

    with open(ipvFile, 'r') as f:
        ipvDF = json.load(f)
    f.close()

    pathologyDF = pd.read_csv(pathologyFile, sep='\t', header=0)

    pathologyPerCoocIndividual = dict()
    for variant in variantsDF['cooccurring vus']:
        for pathogenicVariant in variantsDF['cooccurring vus'][variant]['pathogenic variants']:
            pv = str(tuple(pathogenicVariant))
            heterozygousIndividuals = ipvDF[pv]['heterozygous individuals']
            pathologyPerCoocIndividual[pv] = list()
            pathologies = dict()
            for hi in heterozygousIndividuals:
                hiInt = int(hi)
                row = pathologyDF.loc[pathologyDF['ID'] == hiInt ]
                aao = row['Age at onset'].tolist()
                if aao:
                    pathologies['Age at onset'] = aao[0]
                pathologies['Ovarian cancer history'] = row['Ovarian cancer history'].tolist()
                pathologies['Bilateral breast cancer'] = row['Bilateral breast cancer'].tolist()
                pathologies['Tissue type (3 groups)'] = row['Tissue type (3 groups)'].tolist()
                pathologyPerCoocIndividual[pv].append(pathologies)

    print(json.dumps(pathologyPerCoocIndividual, indent=4, sort_keys=True))

if __name__ == "__main__":
    main()