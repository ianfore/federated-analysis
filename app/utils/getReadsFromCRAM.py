import subprocess
import os
import argparse


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

    viewString = subprocess.check_output(['samtools', 'view', bucket, readRange]).decode('utf-8')

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

    print(variant + ':' + str(reads))


if __name__ == "__main__":
    main()