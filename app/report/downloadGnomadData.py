#!/usr/bin/env python
import requests
import json
import numpy as np
import pandas as pd
import argparse




def fetch(jsondata, url="https://gnomad.broadinstitute.org/api"):
    # The server gives a generic error message if the content type isn't
    # explicitly set
    headers = {"Content-Type": "application/json"}
    response = requests.post(url, json=jsondata, headers=headers)
    json = response.json()
    if "errors" in json:
        raise Exception(str(json["errors"]))
    return json


def transcript_to_variants(transcript_id, dataset, reference_genome):
    """
    Given a transcript, return the list of variants that map to the exons
    of the transcript, and were observed in samples from the indicated
    dataset.
    """
    fmt_graphql = """
    {
        transcript(transcript_id: "%s", reference_genome: %s) {
          variants(dataset: %s) {
            variantId: variantId
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
    req_variantlist = {
        "query": fmt_graphql % (transcript_id, reference_genome, dataset),
        "variables": {}
    }
    response = fetch(req_variantlist)
    return response["data"]["transcript"]["variants"]

def gene_to_coords(gene_id, reference_genome):
    """                                                                         
    Given a gene symbol, return the coordinates.                                
    """
    graphql_query = """                                                         
    {                                                                           
        gene(gene_symbol: "%s", reference_genome: %s) {                         
            chrom                                                               
            start                                                               
            stop                                                                
        }                                                                       
    }"""
    graphql_request = {
        "query": graphql_query % (gene_id, reference_genome),
        "variables": {}
    }
    response = fetch(graphql_request)
    return response["data"]["gene"]

def coords_to_variants(chrom, start, stop, dataset_id, reference_genome):
    region_query = """
                {   region(chrom: "%s", start: %s, stop: %s, reference_genome: %s) {
                    variants(dataset: %s) {
                        variantId
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
    response = requests.post(
        'http://gnomad.broadinstitute.org/api',
        json={ "query": region_query % (chrom, start, stop, reference_genome, dataset_id),
                 "variables": {}},
            headers={"content-type": "application/json"})

    try:
        parse = json.loads(response.text)
        variants = parse['data']['region']['variants']
    except Exception as e:
        print(e)
        return None

    return variants

def gene_to_region_variants(gene_symbol, dataset_id, reference_genome):
    """
    Given a gene name, return the list of variants via a region
    query.  These will mostly be the intronic variants.
    """
    coords = gene_to_coords(gene_symbol, reference_genome)
    variantList = coords_to_variants(coords["chrom"], coords["start"],
                                     coords["stop"], dataset_id, reference_genome)
    #return(set(variantList))
    return variantList


def add_delta_one_assay(fvd, svd, ome):
    myDict = dict()
    #if svd is not None:
    try:
        myDict[ome + "_ac_delta"] = float(fvd["ac"] - svd["ac"])
        myDict[ome + "_an_delta"] = float(fvd["an"] - svd["an"])
        myDict[ome + "_ac_hom_delta"] = float(fvd["ac_hom"] - svd["ac_hom"])
    #else:
    except Exception as e:
        myDict[ome + "_ac_delta"] = '-'
        myDict[ome + "_an_delta"] = '-'
        myDict[ome + "_ac_hom_delta"] = '-'
    return myDict

def add_deltas(full_variant_data, subset_variant_data, this_variant):
    deltas = list()
    deltas.append(add_delta_one_assay(full_variant_data[this_variant]["exome"], subset_variant_data[this_variant]["exome"], 'exome'))
    deltas.append(add_delta_one_assay(full_variant_data[this_variant]["genome"], subset_variant_data[this_variant]["genome"], 'genome'))
    return deltas


def getDeltas(full_variant_data, subset_variant_data):
    variant_details = dict()
    for this_variant in full_variant_data:
        if this_variant in subset_variant_data:
            variant_details[this_variant] = dict()
            v = this_variant.split('-')
            variant_details[this_variant]['chrom'] = str(v[0])
            variant_details[this_variant]['pos'] = str(v[1])
            variant_details[this_variant]['ref'] = str(v[2])
            variant_details[this_variant]['alt'] = str(v[3])
            deltas = add_deltas(full_variant_data, subset_variant_data, this_variant)
            for d in deltas:
                for k in d.keys():
                    variant_details[this_variant][k] = d[k]
    return variant_details


def getVariants(transcript, gene, dataset, reference):
    exonic_set = set()
    exonic_variants_non_topmed = transcript_to_variants(transcript, dataset, reference)
    for l in exonic_variants_non_topmed:
        exonic_set.add(json.dumps(l, sort_keys=True))
    intronic_set = set()
    intronic_variants_non_topmed = gene_to_region_variants(gene, dataset, reference)
    combined_set = exonic_set | intronic_set
    return combined_set

def convertSetToDict(mySet):
    myDict = dict()
    for elt in mySet:
        myElt = eval(elt.replace('null', 'None'))
        myDict[myElt['variantId']] = {'genome': myElt['genome'], 'exome': myElt['exome']}
    return myDict

def writeToOutputFile(myDict, outputFile):
    header= 'alt' + '\t' + 'chrom' + '\t' + 'pos' + '\t' + 'ref' + \
            '\t' + 'exome_ac_hom_delta' + '\t' + 'genome_ac_hom_delta' + \
            '\t' + 'exome_ac_delta' + '\t' + 'genome_ac_delta'
    with open(outputFile, 'w') as f:
        f.write(header + '\n')
        for k in myDict:
            line = myDict[k]
            f.write(line['alt'] + '\t' + line['chrom'] + '\t' + line['pos'] + '\t' + line['ref'] + \
                '\t' + str(line['exome_ac_hom_delta']) + '\t' + str(line['genome_ac_hom_delta'])+ \
                 '\t' + str(line['exome_ac_delta']) + '\t' + str(line['genome_ac_delta']))
            f.write('\n')
    f.close()

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', help='output file base name')
    parser.add_argument('-r', '--release', help='gnomad release version (r2_1 or r3)')
    parser.add_argument('-v', '--version', help='human genome version (GRCh37 or GRCh38)')

    options = parser.parse_args()
    return options

def main():
    release = parse_args().release
    reference_genome = parse_args().version
    outputFile = parse_args().output + "_" + release + "_" + reference_genome + ".tsv"
    non_topmed_dataset = "gnomad_" + release + "_non_topmed"
    full_dataset = "gnomad_" + release
    brca1_transcript ="ENST00000357654"
    #brca1_gene = "ENSG00000012048"
    brca1_gene = "BRCA1"
    brca2_transcript = "ENST00000544455"
    #brca2_gene = "ENSG00000139618"
    brca2_gene = "BRCA2"

    # organize brca1 request
    brca1_variants_non_topmed = getVariants(brca1_transcript, brca1_gene, non_topmed_dataset, reference_genome)
    b1_nontopmed_dict = convertSetToDict(brca1_variants_non_topmed)
    brca1_variants_topmed = getVariants(brca1_transcript, brca1_gene, full_dataset, reference_genome)
    b1_topmed_dict = convertSetToDict(brca1_variants_topmed)
    b1_deltas = getDeltas(b1_topmed_dict, b1_nontopmed_dict)

    # organize brca2 request
    brca2_variants_non_topmed = getVariants(brca2_transcript, brca2_gene, non_topmed_dataset, reference_genome)
    b2_nontopmed_dict = convertSetToDict(brca2_variants_non_topmed)
    brca2_variants_topmed = getVariants(brca2_transcript, brca2_gene, full_dataset, reference_genome)
    brca2_topmed_dict = convertSetToDict(brca2_variants_topmed)
    b2_deltas = getDeltas(brca2_topmed_dict, b2_nontopmed_dict)

    # put brca1 and brca2 results into single dictionary
    deltas = dict()
    deltas.update(b1_deltas)
    deltas.update(b2_deltas)

    # write dictionary to output
    writeToOutputFile(deltas, outputFile)

if __name__ == "__main__":
    main()
