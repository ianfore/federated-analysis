import subprocess
import os
import argparse
import pandas as pd
import subprocess
import json

account = 'jcasalet@ucsc.edu'
project = 'bdcat-cohort-2-fellow'



def getReadDepths(account, project, readRange, bucket, sample, variant):
    login = subprocess.Popen(['gcloud', 'config', 'set', 'account', account], stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)

    config = subprocess.Popen(['gcloud', 'config', 'set', 'project', project], stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)

    token = subprocess.check_output(['gcloud', 'auth', 'application-default', 'print-access-token']).decode(
        "utf-8").rstrip()

    os.environ['GCS_OAUTH_TOKEN'] = token

    readDepth = subprocess.check_output(['samtools', 'depth', '-r', readRange, bucket]).decode("utf-8").rstrip()

    return readDepth

def main():
    parser = argparse.ArgumentParser(usage="getReads account project bucket chromosome position")
    parser.add_argument("--a", dest="account", help="google account name", default=None)
    parser.add_argument("--p", dest="project", help="google project name", default=None)
    parser.add_argument("--b", dest="bucket", help="google bucket name", default=None)
    parser.add_argument("--v", dest="variant", help="variant", default=None)
    options = parser.parse_args()
    account = options.account
    project = options.project
    bucket = options.bucket
    variant = options.variant


    # NWD825720	chr13:32329240:G:A	chr13:32329140-32329340	gs://fc-secure-98adf1aa-0499-46bf-8bb5-d99bc1f4359e/test-crai-cram/NWD825720.b38.irc.v1.cram	gs://fc-secure-98adf1aa-0499-46bf-8bb5-d99bc1f4359e/test-crai-cram/NWD825720.b38.irc.v1.cram.crai
    igvTable = open('igv_table.tsv', 'r')
    lines = igvTable.readlines()

    readsDict = dict()
    for row in lines:
        sample = row.split('\t')[0]
        variant = row.split('\t')[1]
        bucket = row.split('\t')[3]
        if not variant in readsDict:
            readsDict[variant] = dict()
        readsDict[variant][sample] = getReads(account, project, bucket, variant)

    igvTable = open('igv_table.tsv', 'r')
    lines = igvTable.readlines()
    account = 'jcasalet@ucsc.edu'
    project = 'bdcat-cohort-2-fellow'

    # NWD825720	chr13:32329240:G:A	chr13:32329140-32329340	gs://fc-secure-98adf1aa-0499-46bf-8bb5-d99bc1f4359e/test-crai-cram/NWD825720.b38.irc.v1.cram	gs://fc-secure-98adf1aa-0499-46bf-8bb5-d99bc1f4359e/test-crai-cram/NWD825720.b38.irc.v1.cram.crai

    readDepthDict = dict()
    for row in lines:
        variant = row.split('\t')[1]
        if not variant in readDepthDict:
            readDepthDict[variant] = list()
        sample = row.split('\t')[0]
        readRange = row.split('\t')[2]
        bucket = row.split('\t')[3]
        readDepthDict[variant].append([sample, getReadDepths(account, project, readRange, bucket, sample, variant)])

    # "NWD377807": ["chr17:43077166:G:A", "chr17\t43077066\t32\nchr17\t43077067\t32\nchr17\t43077068\t32\n"],

    readDepthFreqDict = dict()
    minDepth = float("inf")
    maxDepth = -1
    for variant in readDepthDict:
        readDepthFreqDict[variant] = dict()
        for i in range(len(readDepthDict[variant])):
            sample = readDepthDict[variant][i][0]
            variantPos = int(variant.split(':')[1])
            depths = readDepthDict[variant][i][1]
            depthArray = depths.split('\n')
            depthSum = 0
            variantDepth = -1
            for depth in depthArray:
                d = int(depth.split('\t')[2])
                if d < minDepth:
                    minDepth = d
                if d > maxDepth:
                    maxDepth = d
                depthSum += d
                pos = int(depth.split('\t')[1])
                if pos == variantPos:
                    variantDepth = d

            avgDepth = float(depthSum) / float(len(depthArray))
            readDepthFreqDict[variant][sample] = {"depth": variantDepth, "minDepth": minDepth,
                                                  "maxDepth": maxDepth, "meanDepth": avgDepth}

    print(readDepthFreqDict)
    with open('readDepthFreqs.json', 'w') as f:
        json.dump(readDepthFreqDict, f)
    f.close()

    # join read dict with depth dict
    for variant in readDepthFreqDict:
        for sample in readDepthFreqDict[variant]:
            readDepthFreqDict[variant][sample]['reads'] = readsDict[variant][sample]

    print(readDepthFreqDict)

def getReads(account, project, bucket, variant):


    # account=jcasalet@ucsc.edu
    # project = bdcat-cohort-2-fellow

    login = subprocess.Popen(['gcloud', 'config', 'set', 'account', account], stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)

    config = subprocess.Popen(['gcloud', 'config', 'set', 'project', project], stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)

    token = subprocess.check_output(['gcloud', 'auth', 'application-default', 'print-access-token']).decode(
        "utf-8").rstrip()

    os.environ['GCS_OAUTH_TOKEN'] = token

    # determine readRange from variant
    chrom = variant.split(':')[0]
    pos = int(variant.split(':')[1])
    ref = variant.split(':')[2]
    alt = variant.split(':')[3]

    lref = len(ref)
    lalt = len(alt)

    if lref == lalt and lref == 1:
        # snp
        readRange = chrom + ":" + str(pos) + "-" + str(pos)

    elif lref > lalt and lalt == 1:
        # deletion
        end = pos + lref - 1
        readRange = chrom + ":" + str(pos) + "-" + str(end)
    elif lref < lalt and lref == 1:
        # insertion
        end = pos + lalt - 1
        readRange = chrom + ":" + str(pos) + "-" + str(end)
    elif lref == lalt and lref != 1:
        # other
        end = pos + lalt - 1
        readRange = chrom + ":" + str(pos) + "-" + str(end)
    else:
        print('huh?')

    viewString = subprocess.check_output(['samtools', 'view', '-F 1024', bucket, readRange]).decode('utf-8')

    viewArray = viewString.split('\n')

    startPos = int(readRange.split(':')[1].split('-')[0])
    endPos = int(readRange.split(':')[1].split('-')[1])
    length = abs(endPos - startPos) + 1

    reads = dict()
    for i in range(len(viewArray) - 1):
        readArray = viewArray[i].split('\t')
        readPos = int(readArray[3])
        offset = abs(readPos - startPos)
        sequence = readArray[9]
        if offset < len(sequence):
            nuc = sequence[offset:offset + length]
            if not nuc in reads:
                reads[nuc] = 0
            reads[nuc] += 1

    return reads

if __name__ == "__main__":
    main()