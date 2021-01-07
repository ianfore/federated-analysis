import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.edit, hgvs.posedit
import hgvs.sequencevariant

class coordinateMapper():
    def __init__(self, genomeAssembly):
        self.genomeAssembly = genomeAssembly
        self.hdp = hgvs.dataproviders.uta.connect()
        self.varmapper = hgvs.assemblymapper.AssemblyMapper(self.hdp, assembly_name=self.genomeAssembly,
                                                   alt_aln_method='splign')

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
            accessioned_chrom = "NC_000017.11"
            ref_cdna = "NM_007294.4"
        else:
            accessioned_chrom = "NC_000013.11"
            ref_cdna = "NM_000059.3"
        start = hgvs.location.BaseOffsetPosition(base=pos)
        end = hgvs.location.BaseOffsetPosition(base=pos + len(ref) - 1)
        iv = hgvs.location.Interval(start=start,end=end)
        edit = hgvs.edit.NARefAlt(ref=ref,alt=alt)
        posedit = hgvs.posedit.PosEdit(pos=iv,edit=edit)
        var_g = hgvs.sequencevariant.SequenceVariant(ac=accessioned_chrom, type='g', posedit = posedit)
        var_c = self.varmapper.g_to_c(var_g,ref_cdna)
        return(str(var_c))

if __name__ == "__main__":
    main()