import json
import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.edit, hgvs.posedit
import hgvs.sequencevariant
import logging
import argparse
import pandas as pd

class coordinateMapper():
    def __init__(self, genomeAssembly):
        self.genomeAssembly = genomeAssembly
        self.hdp = hgvs.dataproviders.uta.connect()
        self.varmapper = hgvs.assemblymapper.AssemblyMapper(self.hdp, assembly_name=self.genomeAssembly,
                                                   alt_aln_method='splign')

        logging.basicConfig()
        self.logger = logging.getLogger()
        self.logger.setLevel(logging.WARN)

    def translate_to_hgvs(self, vartokens):
        # expected input is a tuple of the form (chrom, pos, ref, alt)
        # eg (13, 139929939, 'A', 'C')
        if not isinstance(vartokens, tuple) or len(vartokens) != 4:
            return None
        elif not isinstance(vartokens[0], int) or not isinstance(vartokens[1], int) or \
            not isinstance(vartokens[2], str) or not isinstance(vartokens[3], str):
            return None
        chrom = int(vartokens[0])
        pos = int(vartokens[1])
        ref = vartokens[2]
        alt = vartokens[3]
        if chrom == 17:
            accessioned_chrom = "NC_000017.10"
            ref_cdna = "NM_007294.4"
        else:
            accessioned_chrom = "NC_000013.10"
            ref_cdna = "NM_000059.3"
        var_c = None
        try:
            start = hgvs.location.BaseOffsetPosition(base=pos)
            end = hgvs.location.BaseOffsetPosition(base=pos + len(ref) - 1)
            iv = hgvs.location.Interval(start=start,end=end)
            edit = hgvs.edit.NARefAlt(ref=ref,alt=alt)
            posedit = hgvs.posedit.PosEdit(pos=iv,edit=edit)
            var_g = hgvs.sequencevariant.SequenceVariant(ac=accessioned_chrom, type='g', posedit = posedit)
            var_c = self.varmapper.g_to_c(var_g,ref_cdna)
        except Exception as e:
            self.logger.warning('exception processing ' + str((chrom, start, end)) + ' : got ' + str(var_c))
            pass
        return(str(var_c))

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--intersection', help='intersection json file')
    parser.add_argument('-o', '--out', help='out json file')
    parser.add_argument('-b', '--bayesdel', help='bayesdel vcf file')
    parser.add_argument('-r', '--report', help='report tsv file')
    options = parser.parse_args()
    return options


def main():
    intersectionFileName = parse_args().intersection
    outFileName = parse_args().out
    bayesdelFileName = parse_args().bayesdel
    reportFileName = parse_args().report

    with open(intersectionFileName, 'r') as f:
        intersectionDict = json.load(f)
    f.close()

    with open(outFileName, 'r') as f:
        outDict = json.load(f)
    f.close()

    bayesdelDF = pd.read_csv(bayesdelFileName, header=0, sep='\t')

    mapper = coordinateMapper('GRCh37')

    for coocVariant in outDict['cooccurring vus']:
        if int(outDict['cooccurring vus'][coocVariant]['likelihood data']['k']) < 2:
            continue
        pos = eval(coocVariant)[1]
        hgvsCoord = mapper.translate_to_hgvs(eval(coocVariant))
        caseFreq = intersectionDict['cooccurring'][coocVariant]['caseFreq']
        controlFreq = intersectionDict['cooccurring'][coocVariant]['controlFreq']
        af = outDict['cooccurring vus'][coocVariant]['allele frequencies']
        likelihoodData = outDict['cooccurring vus'][coocVariant]['likelihood data']
        bayesInfo = bayesdelDF[bayesdelDF['POS'] == pos]['INFO']

    for homoVariant in intersectionDict['homozygous']:
        hgvsCoord = mapper.translate_to_hgvs(eval(homoVariant))

if __name__ == "__main__":
    main()