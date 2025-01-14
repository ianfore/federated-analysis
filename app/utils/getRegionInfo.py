import json
import logging
import argparse
from pyliftover import LiftOver


class regionInfo():
    # genomeAssembly = 'GRCh38' or 'GRCh37'
    def __init__(self, regionsFileName, hg):
        with open(regionsFileName, 'r') as f:
            self.regionsDict = json.load(f)
        f.close()
        self.lo = None
        if hg != 'hg38':
            self.lo = LiftOver(hg, 'hg38')


    def get(self, chr, gene, strand, pos):

        # get splice donors
        '''"brca1RefSpliceAcceptorBounds" : {
            "exon2": {
                "start": 43124135,
                "end": 43124113
            },'''
        if strand == 'negative':
            myStrand = '-'
        else:
            myStrand = '+'
        if not self.lo is None:
            pos = self.lo.convert_coordinate(chr, pos ,myStrand)[0][1]
        spliceDonors = list()
        spliceDonorBounds = gene + 'RefSpliceDonorBounds'
        for exon in self.regionsDict[spliceDonorBounds]:
            if strand == 'positive':
                start = int(self.regionsDict[spliceDonorBounds][exon]['start'])
                end = int(self.regionsDict[spliceDonorBounds][exon]['end'])
            elif strand == 'negative':
                start = int(self.regionsDict[spliceDonorBounds][exon]['end'])
                end = int(self.regionsDict[spliceDonorBounds][exon]['start'])
            else:
                print('unknown strand: ' + strand)
                return None
            if pos >= start and pos <= end:
                spliceDonors.append(exon)

        # get splice acceptors
        '''"brca2RefSpliceAcceptorBounds" : {
            "exon2": {
                "start": 32316402,
                "end": 32316424
            },'''
        spliceAcceptorBounds = gene + 'RefSpliceAcceptorBounds'
        spliceAcceptors = list()
        for exon in self.regionsDict[spliceAcceptorBounds]:
            if strand == 'positive':
                start = int(self.regionsDict[spliceAcceptorBounds][exon]['start'])
                end = int(self.regionsDict[spliceAcceptorBounds][exon]['end'])
            elif strand == 'negative':
                start = int(self.regionsDict[spliceAcceptorBounds][exon]['end'])
                end = int(self.regionsDict[spliceAcceptorBounds][exon]['start'])
            else:
                print('unknown strand: ' + strand)
                return None
            if pos >= start and pos <= end:
                spliceAcceptors.append(exon)

        # get exons
        '''"brca2Exons" : {
                "exon1": {
                    "end": 32315667,
                    "start": 32315479
                },'''
        exonBounds = gene + 'Exons'
        exons = list()
        for exon in self.regionsDict[exonBounds]:
            if strand == 'positive':
                start = int(self.regionsDict[exonBounds][exon]['start'])
                end = int(self.regionsDict[exonBounds][exon]['end'])
            elif strand == 'negative':
                start = int(self.regionsDict[exonBounds][exon]['end'])
                end = int(self.regionsDict[exonBounds][exon]['start'])
            else:
                print('unknown strand: ' + strand)
                return None
            if pos >= start and pos <= end:
                exons.append(exon)


        # get CI domains
        geneCIDomains = gene + 'CIDomains'
        inCIdomain = dict()
        for org in self.regionsDict[geneCIDomains]:
            inCIdomain[org] = list()
            for domain in self.regionsDict[geneCIDomains][org]['domains']:
                name = domain['name']
                if strand == 'positive':
                    start = int(domain['start'])
                    end = int(domain['end'])
                elif strand == 'negative':
                    start = int(domain['end'])
                    end = int(domain['start'])
                else:
                    print('unknown strand: ' + strand)
                    return None
                if pos >= start and pos <= end:
                    inCIdomain[org].append(name)

        return {"ciDomains": inCIdomain,
                "spliceDonors": spliceDonors,
                "spliceAcceptors": spliceAcceptors,
                "exons": exons}

def main():
    parser = argparse.ArgumentParser(usage="getRegionInfo")
    parser.add_argument("--r", dest="regionsFile", help="json regions file", default=None)
    parser.add_argument("--g", dest="gene", help="gene name", default=None)
    parser.add_argument("--p", dest="position", help="variant position", default=None)
    parser.add_argument("--s", dest="strand", help="strand (positive or negative", default=None)
    parser.add_argument("--v", dest="version", help="genome assembly version", default=None)
    parser.add_argument("--c", dest="chrom", help="chromosome", default=None)

    options = parser.parse_args()
    regionsFile = options.regionsFile
    chr = options.chrom
    gene = options.gene
    position = int(options.position)
    strand = options.strand
    version = options.version

    ri = regionInfo(regionsFile, version)
    infoDict = ri.get(chr, gene, strand, position)
    print(infoDict)


if __name__ == "__main__":
    main()