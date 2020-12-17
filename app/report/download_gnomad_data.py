#!/usr/bin/env python
import requests
import json
import time
import numpy as np
import pandas as pd
import argparse
from math import floor, log10, isnan


def fetch(jsondata, url="https://gnomad.broadinstitute.org/api"):
    # The server gives a generic error message if the content type isn't
    # explicitly set
    headers = {"Content-Type": "application/json"}
    response = requests.post(url, json=jsondata, headers=headers)
    json = response.json()
    if "errors" in json:
        raise Exception(str(json["errors"]))
    return json

def unique_variant_set(variant_list_of_dicts):
    """ 
    Given a list of dicts, for which the items are variant IDs (as returned by the gnomAD API), 
    generate a set of unique set of variant IDs
    """
    variant_set = set()
    for item in variant_list_of_dicts:
        variant_set.add(list(item.values())[0])
    return variant_set

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
            variant_id: variantId
          }                                                                     
        }                                                                       
      }                                                                         
    """
    req_variantlist = {
        "query": fmt_graphql % (transcript_id, reference_genome, dataset),
        "variables": {}
    }
    response = fetch(req_variantlist)
    return unique_variant_set(response["data"]["transcript"]["variants"])

def gene_to_coords(gene_id, reference_genome):
    """                                                                         
    Given a gene symbol, return the coordinates.                                
    """
    graphql_query = """                                                         
    {                                                                           
        gene(gene_id: "%s", reference_genome: %s) {                         
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
        query GnomadRegion($chrom: String!, $start: Int!, $stop: Int!, 
                           $dataset_id: DatasetId!) {
            region(chrom: $chrom, start: $start, stop: $stop, reference_genome: GRCh38) {
                variants(dataset: $dataset_id) {
                   variantId
                   
                }
            }
         }
         """
    region_variables = {
        "chrom": chrom,
        "start": start,
        "stop": stop,
        "dataset_id": dataset_id,
        "reference_genome": reference_genome
        }
    headers = { "content-type": "application/json" }
    response = requests.post(
        'http://gnomad.broadinstitute.org/api',
        json={ "query": region_query, "variables": region_variables },
        headers=headers)
    try:
        parse = json.loads(response.text)
        variants = parse['data']['region']['variants']
    except json.decoder.JSONDecodeError:
        print('json decode error')
        return None

    return(unique_variant_set(variants))

def gene_to_region_variants(gene_id, dataset_id, reference_genome):
    """
    Given a gene name, return the list of variants via a region
    query.  These will mostly be the intronic variants.
    """
    coords = gene_to_coords(gene_id, reference_genome)
    variantList = coords_to_variants(coords["chrom"], coords["start"],
                                     coords["stop"], dataset_id, reference_genome)
    return(set(variantList))

def fetch_data_for_one_variant(variant_id, dataset, reference_genome, max_retries=5):
    variant_detail_query = """
        query GnomadVariant($variantId: String!, $datasetId: DatasetId!) {
            variant(variantId: $variantId, dataset: $datasetId) {
                alt
                chrom
                pos
                ref
                variantId                                            
                ... on VariantDetails {
                        flags
                        sortedTranscriptConsequences {
                            transcript_id
                            hgvsc
                        }
                    }
                    exome {
                        ac
                        an
                        ac_hom
                        faf95 {
                            popmax
                            popmax_population
                       }
                       filters
                            populations {
                            id
                            ac
                            an
                            ac_hom
                       }                 
                    }
                    genome {
                        ac
                        an
                        ac_hom
                        faf95 {
                            popmax
                            popmax_population
                       }
                       filters
                            populations {
                            id
                            ac
                            an
                            ac_hom
                       }                 
                    }
                }
             }"""
    headers = { "content-type": "application/json" }
    print("Fetching", variant_id)
    variant_detail_variables = {
        "variantId": variant_id,
        "datasetId": dataset,
        "referenceGenome": reference_genome,
    }
    retries = 0
    while retries < max_retries:
        parse = None
        try:
            # https://stackoverflow.com/questions/49064398/requests-exceptions-chunkedencodingerror-connection-broken-incompleteread0
            with requests.post(
                'http://gnomad.broadinstitute.org/api',
                json={
                    "query": variant_detail_query,
                    "variables": variant_detail_variables
                },
                headers=headers) as response:
                parse = json.loads(response.text)
                if response.status_code == 400:
                    print(response.reason)
        except json.decoder.JSONDecodeError:
            print('json decode error')
            retries += 1
            time.sleep(0.1)
            continue
        except requests.exceptions.ConnectionError:
            print('connection error')
            retries += 1
            time.sleep(0.1)
            continue
        except requests.exceptions.StreamConsumedError:
            print('stream consumed error')
            retries += 1
            time.sleep(0.1)
            continue
        except requests.exceptions.RequestException:
            print('request exception')
            retries += 1
            time.sleep(0.1)
            continue
        except Exception as e:
            print('generic exception: ' + str(e))
            return None

        print(parse)
        if not parse is None:
            return(parse['data']['variant'])
        else:
            return None


def add_delta_one_assay(fvd, svd):
    if fvd is not None:
        if svd is not None:
            fvd["ac_delta"] = fvd["ac"] - svd["ac"]
            fvd["an_delta"] = fvd["an"] - svd["an"]
            fvd["ac_hom_delta"] = fvd["ac_hom"] - svd["ac_hom"]
            if "populations" in fvd:
                for counter  in range(0, len(fvd["populations"])):
                    full_pop = fvd["populations"][counter]
                    subset_pop = svd["populations"][counter]
                    full_pop = add_delta_one_assay(full_pop, subset_pop)
        else:
            fvd["ac_delta"] = 0
            fvd["an_delta"] = 0
            fvd["ac_hom_delta"] = 0
            if "populations" in fvd:
                for counter  in range(0, len(fvd["populations"])):
                    full_pop = fvd["populations"][counter]
                    full_pop = add_delta_one_assay(full_pop, None)
    return fvd

def add_deltas(full_variant_data, subset_variant_data):
    if "exome" in full_variant_data:
        full_variant_data["exome"] = add_delta_one_assay(full_variant_data["exome"], subset_variant_data["exome"])
    if "genome" in full_variant_data:
        full_variant_data["genome"] = add_delta_one_assay(full_variant_data["genome"], subset_variant_data["genome"])
    return full_variant_data


def variant_set_to_variant_data(variants, full_dataset, subset_dataset, reference_genome):
    variant_details = []
    for this_variant in variants:
        time.sleep(0.5)
        full_variant_data = fetch_data_for_one_variant(this_variant, full_dataset, reference_genome)
        if full_variant_data is not None:
            subset_variant_data = fetch_data_for_one_variant(this_variant, subset_dataset, reference_genome)
            if subset_variant_data is not None:
                full_variant_data = add_deltas(full_variant_data, subset_variant_data)
                variant_details.append(full_variant_data)
                #time.sleep(0.01)
        print(full_variant_data)
    return(variant_details)


def compile_allele_values(df):
    populations = ['AFR', 'AFR_FEMALE', 'AFR_MALE', 'AMR', 'AMR_FEMALE', 'AMR_MALE', 'ASJ', 'ASJ_FEMALE', 'ASJ_MALE', 'EAS',
                   'EAS_JPN', 'EAS_KOR', 'EAS_OEA', 'EAS_FEMALE', 'EAS_MALE', 'FIN', 'FIN_FEMALE', 'FIN_MALE', 'NFE', 'NFE_BGR',
                   'NFE_EST', 'NFE_NWE', 'NFE_ONF', 'NFE_SEU', 'NFE_SWE', 'NFE_FEMALE', 'NFE_MALE', 'OTH', 'OTH_FEMALE', 'OTH_MALE',
                   'SAS', 'SAS_FEMALE', 'SAS_MALE', 'FEMALE', 'MALE']
    for population in populations:
        df['genome_' + population + '_af'] = calculate_frequency(df['genome_' + population + '_ac'], df['genome_' + population + '_an'])
        df['exome_' + population + '_af'] = calculate_frequency(df['exome_' + population + '_ac'], df['exome_' + population + '_an'])
    df['exome_af'] = calculate_frequency(df['exome_ac'], df['exome_an'])
    df['genome_af'] = calculate_frequency(df['genome_ac'], df['genome_an'])
    return df


def round_popmax(df):
    df['exome_popmax'] = pd.to_numeric(df['exome_popmax'], errors='coerce').apply(round_four_sigfigs)
    df['genome_popmax'] = pd.to_numeric(df['genome_popmax'], errors='coerce').apply(round_four_sigfigs)
    return df


def calculate_frequency(ac, an):
    freq = pd.to_numeric(ac, errors='coerce').divide(pd.to_numeric(an, errors='coerce'))
    return freq.apply(round_four_sigfigs)


def round_four_sigfigs(num):
    if isnan(num):
        return num
    elif num == 0 or num == 0.0:
        return 0
    else:
        return round(num, -int(floor(log10(abs(num))) - (3)))


def flatten(variant, field, genome_or_exome):
    for f in field:
        if f not in ['populations', 'filters']:
            if f == 'faf95':
                for p in field[f]:
                    variant[genome_or_exome + '_' + p] = field[f][p]
            else:
                variant[genome_or_exome + '_' + f] = field[f]
    populations = field['populations']
    for population in populations:
        name = population['id']
        keys = population.keys()
        for key in keys:
            if name != key:
                variant[genome_or_exome + '_' + name + '_' + key] = population[key]
    return variant


def flatten_populations(variants):
    for variant in variants:
        genome = variant['genome']
        exome = variant['exome']
        if genome:
            variant = flatten(variant, genome, 'genome')
        if exome:
            variant = flatten(variant, exome, 'exome')
        del variant['genome']
        del variant['exome']
    return variants


def find_correct_hgvs(variants, transcripts):
    """
    Given the set of transcript IDs that we queried for (one per gene), and
    given the data for one particular variant, return the cDNA HGVS string
    for that variant, corresponding to one of the transcripts we're intereted
    in.  There should be one and only one.  If there is no such transcript
    and HGVS string for this variant, then make a note and throw out the variant.
    """
    variants_in_expected_transcripts = []
    for variant in variants:
        sortedTranscriptConsequences = variant["sortedTranscriptConsequences"]
        for transcript_hgvs in sortedTranscriptConsequences:
            if transcript_hgvs["transcript_id"] in transcripts:
                variant["hgvs"] = transcript_hgvs["hgvsc"].split(':')[1]
                variant["transcript"] = transcript_hgvs["transcript_id"]
                del variant["sortedTranscriptConsequences"]
                variants_in_expected_transcripts.append(variant)
                break
        if "hgvs" not in variant:
            print("Warning: variant %s falls outside the expected transcripts"
                  % variant["variantId"])
    return variants_in_expected_transcripts


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', help='Ouput tsv file result.')
    options = parser.parse_args()
    return options


def main():
    f_out = parse_args().output
    #non_topmed_dataset = "gnomad_r2_1_non_topmed"
    non_topmed_dataset = "gnomad_r3_non_topmed"
    #full_dataset = "gnomad_r2_1"
    full_dataset = "gnomad_r3"
    brca1_transcript ="ENST00000357654"
    brca1_gene = "ENSG00000012048"
    brca2_transcript = "ENST00000544455"
    brca2_gene = "ENSG00000139618"
    reference_genome = "GRCh38"
    transcripts = (brca1_transcript, brca2_transcript)

    # organize brca1 request
    brca1_exonic_variants = transcript_to_variants(brca1_transcript, non_topmed_dataset, reference_genome)
    brca1_intronic_variants = gene_to_region_variants(brca1_gene, non_topmed_dataset, reference_genome)
    brca1_variants = brca1_intronic_variants | brca1_exonic_variants

    # organize brca2 request
    brca2_exonic_variants = transcript_to_variants(brca2_transcript, non_topmed_dataset, reference_genome)
    brca2_intronic_variants = gene_to_region_variants(brca2_gene, non_topmed_dataset, reference_genome)
    brca2_variants = brca2_intronic_variants | brca2_exonic_variants

    # combine requests and get brca1 and brca2 data from gnomAD
    brca12_variants = brca1_variants | brca2_variants
    brca12_variant_data = variant_set_to_variant_data(brca12_variants, full_dataset, non_topmed_dataset, reference_genome)

    # find hgvs, flatten, convert to dataframe, compute allele frequencies, and normalize
    variants_with_hgvs = find_correct_hgvs(brca12_variant_data, transcripts)
    variants_with_flattened_populations = flatten_populations(variants_with_hgvs)
    variants_df = pd.json_normalize(variants_with_flattened_populations)
    variants_df['flags'] = variants_df['flags'].apply(', '.join)
    df_with_allele_values = compile_allele_values(variants_df)
    df_with_rounded_popmax = round_popmax(df_with_allele_values)
    stringified_df_with_allele_values = df_with_rounded_popmax.replace(np.nan, '-', regex=True).replace('', '-', regex=True)

    # output to .tsv
    stringified_df_with_allele_values.to_csv(f_out, sep='\t', index=False)


if __name__ == "__main__":
    main()
