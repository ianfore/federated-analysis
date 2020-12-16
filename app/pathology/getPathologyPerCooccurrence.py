import json
import pandas as pd
import sys

def main():
    if len(sys.argv) != 5:
        print(
            'chr-out.json chr-ipv.json shuffle.tsv output.json')
        sys.exit(1)
    variantsFile = sys.argv[1]
    ipvFile = sys.argv[2]
    pathologyFile = sys.argv[3]
    outputFile = sys.argv[4]

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
                else:
                    pathologies['Age at onset'] = 0.0
                pathologies['Ovarian cancer history'] = row['Ovarian cancer history'].tolist()
                pathologies['Bilateral breast cancer'] = row['Bilateral breast cancer'].tolist()
                pathologies['Tissue type (3 groups)'] = row['Tissue type (3 groups)'].tolist()
                pathologies['TMN classification / T'] = row['TMN classification / T'].tolist()
                pathologies['TNM classification / N'] = row['TNM classification / N'].tolist()
                pathologies['TNM classification / M'] = row['TNM classification / M'].tolist()
                pathologies['ER'] = row['ER'].tolist()
                pathologies['PgR'] = row['PgR'].tolist()
                pathologies['HER2'] = row['HER2'].tolist()
                
                pathologyPerCoocIndividual[pv].append(pathologies)

    pathologyPerHomoIndividual = dict()
    for variant in variantsDF['homozygous vus']:
        for pathogenicVariant in variantsDF['cooccurring vus'][variant]['pathogenic variants']:
            pv = str(tuple(pathogenicVariant))
            homozygousIndividuals = ipvDF[pv]['homozygous individuals']
            pathologyPerHomoIndividual[pv] = list()
            pathologies = dict()
            for hi in homozygousIndividuals:
                hiInt = int(hi)
                row = pathologyDF.loc[pathologyDF['ID'] == hiInt]
                aao = row['Age at onset'].tolist()
                if aao:
                    pathologies['Age at onset'] = aao[0]
                else:
                    pathologies['Age at onset'] = 0.0
                pathologies['Ovarian cancer history'] = row['Ovarian cancer history'].tolist()
                pathologies['Bilateral breast cancer'] = row['Bilateral breast cancer'].tolist()
                pathologies['Tissue type (3 groups)'] = row['Tissue type (3 groups)'].tolist()
                pathologies['TMN classification / T'] = row['TMN classification / T'].tolist()
                pathologies['TNM classification / N'] = row['TNM classification / N'].tolist()
                pathologies['TNM classification / M'] = row['TNM classification / M'].tolist()
                pathologies['ER'] = row['ER'].tolist()
                pathologies['PgR'] = row['PgR'].tolist()
                pathologies['HER2'] = row['HER2'].tolist()

                pathologyPerHomoIndividual[pv].append(pathologies)

    pathologyPerAllIndividuals = dict()
    #pathologyPerAllIndividuals.update(pathologyPerHomoIndividual)
    #pathologyPerAllIndividuals.update(pathologyPerCoocIndividual)
    pathologyPerAllIndividuals['homozygous'] = pathologyPerHomoIndividual
    pathologyPerAllIndividuals['cooccurring'] = pathologyPerCoocIndividual

    '''json_dump = json.dumps(pathologyPerAllIndividuals, sort_keys=True)
    with open(outputFile, 'w') as f:
        f.write(json_dump)
    f.close()'''

if __name__ == "__main__":
    main()