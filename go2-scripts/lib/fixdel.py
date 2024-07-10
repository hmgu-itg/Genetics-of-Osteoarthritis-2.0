#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import os
import pysam
import gzip

def check_file_writable(fnm):
    if os.path.exists(fnm):
        # path exists
        if os.path.isfile(fnm): # is it a file or a dir?
            # also works when file is a link and the target is writable
            return os.access(fnm, os.W_OK)
        else:
            return False # path is a dir, so cannot write as a file
    # target does not exist, check perms on parent dir
    pdir = os.path.dirname(fnm)
    if not pdir: pdir = '.'
    # target is creatable if parent dir is writable
    return os.access(pdir, os.W_OK)

def fix_record(record: pysam.libcbcf.VariantRecord) -> pysam.libcbcf.VariantRecord:
    """
    If the ALT column of the records contains a `-`, the variant position is shift one position to the left.
    For reaching this the reference base for the position is determined and prepended to the current REF bases.
    The ALT value `-` is replaced with the determined ref base on position before. In case there are multiple ALT
    values the determined ref base is prepended to them.

    :param record: vcf file record
    :return: fixed vcf record
    """
    global iter
    iter+=1
    if iter % 1000 == 0:
        sys.stderr.write("processed "+str(iter)+" lines.     \r")
        sys.stderr.flush()
    if "<DEL>" in record.alts:
        ref = fasta.fetch(region=f"{record.chrom}:{record.pos - 1}-{record.pos - 1}").upper()
        record.pos = record.pos - 1
        record.ref = f"{ref}{record.ref}"
        record.alts = tuple(map(lambda alt: ref if alt == "-" else f"{ref}{alt}", record.alts))
    if record.ref == "<DEL>":
        ref = fasta.fetch(region=f"{record.chrom}:{record.pos - 1}-{record.pos - 1}").upper()
        record.pos = record.pos - 1
        record.ref = f"{ref}"
        record.alts = tuple(map(lambda alt: f"{ref}{alt}", record.alts))

    return record


vcf=sys.argv[1]
if (not os.path.exists(vcf)) or os.path.getsize(vcf) == 0:
    print("Error: Input file: "+vcf+" does not exist or is not empty")
    exit(1)

try:
    exec_path=os.path.dirname(os.path.realpath(__file__))
    config_file=open(exec_path+"/../config_liftnorm.txt", "r");
    ref=dict(v.split(";") for v in config_file.read().split("\n"))['reference']
    config_file.close();
except Exception as e:
    print("Unable to extract the reference field from the config file. Make sure it is located in the parent directory of this script.")
    exit(1)

if (not os.path.exists(ref)) or os.path.getsize(ref) == 0:
    print("Error: Reference file : "+ref+" found in config but it does not exist.")
    exit(1)


# out=sys.argv[2]
# if (not check_file_writable(out)) :
#     print("Error: Output file: "+out+" is not writeable.")
#     exit(1)

sys.stdin = open(vcf, "r")
# sys.stdout = gzip.open(out, "wt")
vcf = pysam.VariantFile(sys.stdin)
fasta = pysam.FastaFile(ref)
print(vcf.header, end="")
iter=0

print(*map(fix_record, vcf), sep="", end="")
sys.stderr.write("\nDone.\n")
