#!/usr/bin/env python

import requests
import pprint
from collections import defaultdict
import argparse
import sys

prettyprint = pprint.PrettyPrinter(indent=2).pprint

def fetch(jsondata, url="https://gnomad.broadinstitute.org/api"):
    # The server gives a generic error message if the content type isn't explicitly set
    headers = {"Content-Type": "application/json"}
    response = requests.post(url, json=jsondata, headers=headers)
    json = response.json()
    if "errors" in json:
        raise Exception(str(json["errors"]))
    return json

def get_variant_list(gene_symbol, dataset, reference):
    # Note that this is GraphQL, not JSON.
    fmt_graphql = """
    {
        gene(gene_symbol: "%s", reference_genome: %s) {
          variants(dataset: %s) {
            variant_id: variantId
            exome {
		ac
		an
		ac_hom
	    }
            genome {
		ac
		an
		ac_hom
	    }

          }
        }
      }
    """
    # This part will be JSON encoded, but with the GraphQL part left as a glob of text.
    req_variantlist = {
        "query": fmt_graphql % (gene_symbol, reference, dataset),
        "variables": {}
        }
    response = fetch(req_variantlist)
    return response["data"]["gene"]["variants"]

def generate_variant_dict(gene_symbol, with_dataset, without_dataset, reference):
    theDict = defaultdict(dict)
    with_topmed_list = get_variant_list(gene_symbol, with_dataset, reference)
    without_topmed_list = get_variant_list(gene_symbol, without_dataset, reference)

    all_list = with_topmed_list + without_topmed_list
    for v in all_list:
        theDict[v['variant_id']]['genome_ac_hom_with'] = '-'
        theDict[v['variant_id']]['genome_ac_hom_without'] = '-'
        theDict[v['variant_id']]['exome_ac_hom_with'] = '-'
        theDict[v['variant_id']]['exome_ac_hom_without'] = '-'
        theDict[v['variant_id']]['genome_ac_with'] = '-'
        theDict[v['variant_id']]['genome_ac_without'] = '-'
        theDict[v['variant_id']]['exome_ac_with'] = '-'
        theDict[v['variant_id']]['exome_ac_without'] = '-'


    # { 'exome': None, 'genome': {'ac': 2, 'ac_hom': 0, 'an': 80488}, 'variant_id': '17-43125258-G-A'}

    # construct dict of variants that ARE from topmed
    for v in with_topmed_list:
        if not v['genome'] is None:
            theDict[v['variant_id']]['genome_ac_hom_with'] = v['genome']['ac_hom']
            theDict[v['variant_id']]['genome_ac_with'] = v['genome']['ac']
        if not v['exome'] is None:
            theDict[v['variant_id']]['exome_ac_hom_with'] = v['exome']['ac_hom']
            theDict[v['variant_id']]['exome_ac_with'] = v['exome']['ac']

    # construct dict of variants that are NOT from topmed
    for v in without_topmed_list:
        if not v['genome'] is None:
            theDict[v['variant_id']]['genome_ac_hom_without'] = v['genome']['ac_hom']
            theDict[v['variant_id']]['genome_ac_without'] = v['genome']['ac']
        if not v['exome'] is None:
            theDict[v['variant_id']]['exome_ac_hom_without'] = v['exome']['ac_hom']
            theDict[v['variant_id']]['exome_ac_without'] = v['exome']['ac']

    # construct delta of counts per variant
    for v in theDict:
        theDict[v]['genome_ac_hom_delta'] = get_delta(theDict[v], 'genome_ac_hom_with', 'genome_ac_hom_without')
        theDict[v]['exome_ac_hom_delta'] = get_delta(theDict[v], 'exome_ac_hom_with', 'exome_ac_hom_without')
        theDict[v]['genome_ac_delta'] = get_delta(theDict[v], 'genome_ac_with', 'genome_ac_without')
        theDict[v]['exome_ac_delta'] = get_delta(theDict[v], 'exome_ac_with', 'exome_ac_without')

    return theDict

def get_delta(d, k1, k2):
    try:
        return d[k1] - d[k2]
    except:
        return '-'

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', help='output tsv file result')
    parser.add_argument('-g', '--gene', help='gene (brca1 or brca2)')
    parser.add_argument('-v', '--version', help='gnomad version (r2_1 or r3)')
    parser.add_argument('-r', '--reference', help='hg reference (GRCh37 or GRCh38)')
    options = parser.parse_args()
    return options


def main():
    geneName = parse_args().gene
    version = parse_args().version
    reference = str(parse_args().reference)
    outputFile = parse_args().output + '_' + geneName + '_' + version + '_' + reference + '.tsv'
    dataset = "gnomad_" + version
    theDict = generate_variant_dict(geneName, dataset, dataset + "_non_topmed", reference)
    print('writing output')
    with open(outputFile, 'w') as f:
        f.write('variant' + '\t' +
                'genome_ac_hom_delta' + '\t' +
                'exome_ac_hom_delta' + '\t'
                'genome_ac_delta' + '\t' +
                'exome_ac_delta')
        f.write('\n')
        for v in theDict:
            f.write(str(v) + '\t' +
                    str(theDict[v]['genome_ac_hom_delta']) + '\t' +
                    str(theDict[v]['exome_ac_hom_delta']) + '\t' +
                    str(theDict[v]['genome_ac_delta']) + '\t' +
                    str(theDict[v]['exome_ac_delta']))
            f.write('\n')
    f.close()

if __name__ == "__main__":
    main()