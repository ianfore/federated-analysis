#!/usr/bin/env python
"""
Given container output, add the cDNA HGVS nomenclature strings for each variant
"""
import argparse
import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.edit, hgvs.posedit
import hgvs.sequencevariant
import json
import re

def variant_string_to_list(variant_string):
    """
    Given a variant in string representation, translate it to a list
    """
    varlist = list()
    vartokens = re.split(", ", re.sub("[()]", "", variant_string))
    chrom = int(vartokens[0])
    varlist.append(chrom)
    pos = int(vartokens[1])
    varlist.append(pos)
    ref = re.sub("[']", "", vartokens[2])
    varlist.append(ref)
    alt = re.sub("[']", "", vartokens[3])
    varlist.append(alt)
    return(varlist)
    
    
def translate_to_hgvs(vartokens, varmapper):
    """ 
    Given the variant expressed as a list translated to a string, 
    return the cDNA HGVS string for that variant
    """
    chrom = vartokens[0]
    pos = vartokens[1]
    ref = vartokens[2]
    alt = vartokens[3]
    if chrom == 17:
        accessioned_chrom = "NC_000017.10"
        ref_cdna = "NM_007294.3"
    else:
        accessioned_chrom = "NC_000013.10"
        ref_cdna = "NM_000059.3"
    start = hgvs.location.BaseOffsetPosition(base=pos)
    end = hgvs.location.BaseOffsetPosition(base=pos + len(ref) - 1)
    iv = hgvs.location.Interval(start=start,end=end)
    edit = hgvs.edit.NARefAlt(ref=ref,alt=alt)
    posedit = hgvs.posedit.PosEdit(pos=iv,edit=edit)
    var_g = hgvs.sequencevariant.SequenceVariant(ac=accessioned_chrom,
                                                 type='g', posedit = posedit)
    var_c = varmapper.g_to_c(var_g,ref_cdna)
    return(str(var_c))

def add_hgvs_to_variant_tuple(variant, varmapper):
    varlist = variant_string_to_list(variant)
    hgvs_str = translate_to_hgvs(varlist, varmapper)
    varlist.append(hgvs_str)
    updated_variant = varlist
    return(updated_variant)
    

def cooccurring_variant(variant, variant_data, varmapper, debug=False):
    if debug:
        print("working on", variant)
    updated_variant = add_hgvs_to_variant_tuple(variant, varmapper)
    new_variant_data = dict()
    new_variant_data["pathogenic_variants"] = list()
    for this_pathogenic in variant_data["pathogenic variants"]:
        if debug:
            print("Working on pathogenic variant", this_pathogenic)
        pathogenic_variant_hgvs = translate_to_hgvs(this_pathogenic,
                                                    varmapper)
        this_pathogenic.append(pathogenic_variant_hgvs)
        new_variant_data["pathogenic_variants"].append(this_pathogenic)
    return(updated_variant, new_variant_data)

def cooccurring_variant_set(variant_set, varmapper, debug=False):
    new_vardata = dict()
    for variant in variant_set.keys():
        (updated_variant,
         updated_variant_data) = cooccurring_variant(variant,
                                                     variant_set[variant],
                                                     varmapper, debug)
        new_vardata[str(tuple(updated_variant))] = updated_variant_data
    return(new_vardata)


def homozygous_variant_set(variant_set, varmapper, debug=False):
    new_vardata = list()
    for variant in variant_set:
        if debug:
            print("working on", variant)
        updated_variant = add_hgvs_to_variant_tuple(variant, varmapper)
        new_vardata.append(updated_variant)
    return(new_vardata)

        
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
                        help="input JSON file")
    parser.add_argument("-o", "--output",
                        help="Output JSON file")
    parser.add_argument("-f", "--field",
                        help="Field to remove")
    parser.add_argument("-d", "--debug", default=False,
                        help="Print debugging info")
    args = parser.parse_args()

    hdp = hgvs.dataproviders.uta.connect()
    varmapper = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh37',
                                                   alt_aln_method='splign')
    with open(args.input, 'r') as fin:
        vardata = json.load(fin)
        new_cooccur_set = cooccurring_variant_set(vardata["cooccurring vus"],
                                                  varmapper, args.debug)
        vardata["cooccurring vus"] = new_cooccur_set
        new_homozygous_set = homozygous_variant_set(vardata["homozygous vus"],
                                                    varmapper, args.debug)
        vardata["homozygous vus"] = new_homozygous_set

    with open(args.output, 'w') as fout:
        json.dump(vardata, fout, indent=4)
        

if __name__ == "__main__":
    main()
