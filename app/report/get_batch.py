#!/usr/bin/env python

import requests
import pprint
from collections import defaultdict

prettyprint = pprint.PrettyPrinter(indent=2).pprint

def fetch(jsondata, url="https://gnomad.broadinstitute.org/api"):
    # The server gives a generic error message if the content type isn't
    # explicitly set
    headers = {"Content-Type": "application/json"}
    response = requests.post(url, json=jsondata, headers=headers)
    json = response.json()
    if "errors" in json:
        raise Exception(str(json["errors"]))
    return json

def get_variant_list(gene_id, dataset):
    # Note that this is GraphQL, not JSON.
    fmt_graphql = """
    {
        gene(gene_id: "%s", reference_genome: GRCh38) {
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
    # This part will be JSON encoded, but with the GraphQL part left as a
    # glob of text.
    req_variantlist = {
        "query": fmt_graphql % (gene_id, dataset),
        "variables": {}
        }
    response = fetch(req_variantlist)
    return response["data"]["gene"]["variants"]

# brca1
b1_dict = defaultdict() 

b1_with_list = get_variant_list("ENSG00000012048", "gnomad_r3")
# { 'exome': None, 'genome': {'ac': 2, 'ac_hom': 0, 'an': 80488}, 'variant_id': '17-43125258-G-A'} 
for v in b1_with_list:
    b1_dict[v['variant_id']] = dict()
    b1_dict[v['variant_id']]['ac_hom_with'] = v['genome']['ac_hom']  

 
b1_without_list = get_variant_list("ENSG00000012048", "gnomad_r3_non_topmed")
for v in b1_without_list:
    if v['variant_id'] not in b1_dict:
        b1_dict[v['variant_id']] = dict()
    b1_dict[v['variant_id']]['ac_hom_without'] = v['genome']['ac_hom']  

print(b1_dict)

# brca2

b2_dict = defaultdict() 

b2_with_list = get_variant_list("ENSG00000139618", "gnomad_r3")
# { 'exome': None, 'genome': {'ac': 2, 'ac_hom': 0, 'an': 80488}, 'variant_id': '17-43125258-G-A'} 
for v in b2_with_list:
    b2_dict[v['variant_id']] = dict()
    b2_dict[v['variant_id']]['ac_hom_with'] = v['genome']['ac_hom']  

 
b2_without_list = get_variant_list("ENSG00000139618", "gnomad_r3_non_topmed")
for v in b2_without_list:
    if v['variant_id'] not in b1_dict:
        b2_dict[v['variant_id']] = dict()
    b2_dict[v['variant_id']]['ac_hom_without'] = v['genome']['ac_hom']  

print(b2_dict)
