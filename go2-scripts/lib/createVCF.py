#!/usr/bin/env python3

## Dependencies
try:
    import re
    import os
    import sys
    import gzip
    import argparse
    import pandas as pd
    import csv
    from functools import partial
except ImportError:
    print("Check and install modules: os, sys, re, argparse, sh, gzip,\
            functools, csv, sh")
    sys.exit('ImportError')

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="Generate VCF from summary statisctics file from EasyQC")
parser.add_argument('statsfile',help='stats summary file to convert to vcf file (REQUIRED)')
parser.add_argument('-o', '--outdir', required=True, help="output directory (must exist and be writable)")
parser.add_argument('-w', '--overwrite', action="store_true", help='overwrite existing file', default=False)
parser.add_argument('-n', '--do-not-sort', action="store_true",
                   help='Assume that the input file is already sorted by chromosome, then position.')
parser.add_argument('-u', '--useheader', action="store_true", help='do not use list2 of expected headers two')
args = parser.parse_args()


def template_vcf(outFilename, typeChr='autosome'):
    """
    Creates a dummy vcf file with header and default/fixed columns
    VCF structure based on VCFV4.2 specifications
    CrossMAp requires version 4.2 for lifting over

    outFilename: filename to write output
                 text sting from args
    """

    ## Write header and contigs
    with gzip.open(outFilename, 'wt') as fout:

        # ## Write VCF header
        write = partial(print, file=fout)
        write('##fileformat=VCFv4.2')
        write('##FILTER=<ID=PASS,Description="All filters passed">')
        write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">')
        write('##INFO=<ID=FLIPPED,Number=1,Type=Integer,Description="REF allele flip check">')
        write('##INFO=<ID=FLIPPED_lo,Number=1,Type=Integer,Description="REF allele flip check for liftedover variants">')
        write('##INFO=<ID=STRAND_lo,Number=1,Type=Integer,Description="Allele fixed to opposing strand (+ to - by CrossMap)">')
        write('##INFO=<ID=OLD,Number=1,Type=String,Description="REF/ALT pre setref with bcftools norm step1">')
        write('##contig=<ID=chr1>')
        write('##contig=<ID=chr2>')
        write('##contig=<ID=chr3>')
        write('##contig=<ID=chr4>')
        write('##contig=<ID=chr5>')
        write('##contig=<ID=chr6>')
        write('##contig=<ID=chr7>')
        write('##contig=<ID=chr8>')
        write('##contig=<ID=chr9>')
        write('##contig=<ID=chr10>')
        write('##contig=<ID=chr11>')
        write('##contig=<ID=chr12>')
        write('##contig=<ID=chr13>')
        write('##contig=<ID=chr14>')
        write('##contig=<ID=chr15>')
        write('##contig=<ID=chr16>')
        write('##contig=<ID=chr17>')
        write('##contig=<ID=chr18>')
        write('##contig=<ID=chr19>')
        write('##contig=<ID=chr20>')
        write('##contig=<ID=chr21>')
        write('##contig=<ID=chr22>')
        write('##contig=<ID=chrX>')

def check_header(inFilename, headerTwo=False):
    if not headerTwo:
        expected_header = ['CPTID', 'CHR', 'POS', 'EA', 'NEA', 'EAF']
    else:
        expected_header = ['CHROM', 'POS', 'CPTID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    with gzip.open(inFilename, 'rt') as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        header_stats = reader.fieldnames
        header_stats=[x.upper() for x in header_stats]
        if not set(expected_header).issubset(set(header_stats)):
            print('Warning: Missing expected columns in header')
            print(f'Input file:{inFilename} \n')
            print(f'Columns Missing: ', list(set(expected_header) - set(header_stats)))
            print(f'Headers provided: ', header_stats)
            sys.exit(1)



def append_vcf(outFilename, statsData, noSortFile=False):
    """
    Append variants in stats summary data to vcf file template (outFilename)

    outFilename: filename to write output
                 text string from args (output from template_vcf)

    statsData: .gz compressed file
             : tab delimited
    """

    # Sort by CHROM and POS
    statsData.columns = [x.upper() for x in statsData.columns]
    statsData = statsData.loc[:, ['CHR', 'POS', 'CPTID', 'EA', 'NEA', 'EAF']]
    statsData['EAF'] = statsData['EAF'].round(3)

    if not noSortFile:
        print("Info: sorting by chr and pos...")
        statsData.sort_values(['CHR', 'POS'], inplace=True)

    statsData['AF'] = "AF="+statsData['EAF'].astype(str)

    # checking
    chroms_unique=[str(i) for i in list(statsData.CHR.unique())]
    if not any([re.match("^chr", c) for c in chroms_unique]):
        # none of the chr have the chr prefix
        statsData['CHR']= 'chr' + statsData['CHR'].astype(str)
    chroms_unique=[str(i) for i in list(statsData.CHR.unique())]
    if not all([re.match("^chr", c) for c in chroms_unique]):
        print("Error: some chromosomes have the \"chr\" prefix but not others:", [i for idx,i in enumerate(a) if [not bool(re.match("^chr", b)) for b in a][idx]])
        exit(1)

    # rename columns
    statsData.rename(columns={"CHR": "#CHROM", "CPTID": "ID", "NEA": "REF",
                       "EA": "ALT", "AF": "INFO"}, inplace=True)
    # Add QUAL and PASS
    statsData['QUAL'] = 100
    statsData['FILTER'] = "PASS"
    statsData['POS'] = statsData['POS'].astype(int)
    # re-order columns
    statsData = statsData[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']]

    print('Info: writing to VCF file')
    statsData.to_csv(f'{outFilename}', sep='\t', index=False, mode='a', header=True,
             compression='gzip')
    print('Info: VCF conversion completed')
    print(f'VCF file name:{outFilename}')
    print('Done\n')


#----------------------------------------------------------------------
def main(args):
    """
    Call and create vcf template, populate
    template with variants data
    """
    ### Create output filename
    if args.statsfile:
        if re.search(r'.gz$', args.statsfile):
            outFilename=re.split(r'.gz$', args.statsfile)[0] + '.vcf.gz'
        else:
            print('Warning: Summary statistics file name should end with ".gz"')
            sys.exit(NameError)
    else:
        sys.exit("ERROR: No summary statistics file provided.")

    ### Check if directory is writable
    if not os.access(args.outdir, os.W_OK):
        print("ERROR: directory "+args.outdir+" is not writable.")
        exit()

    ### Resolve path, extract filename only and redirect to outdir
    outFilename=os.path.basename(outFilename)
    outFilename=args.outdir+"/"+outFilename
    print("Info: writing to "+outFilename)
    if args.overwrite:
        print("Info: will overwrite "+outFilename+" if it exists.")
    if args.do_not_sort:
        print("Info: will not re-sort the summary statistics before creating VCF.")
    ### create VCF template
    if (not os.path.isfile(outFilename)) or args.overwrite:
        #print("Info: checking headers.")
        check_header(args.statsfile, headerTwo=args.useheader)
        #print("Info: headers OK.")
        #print("Info: reading in summary statistics.")
        datain=pd.read_csv(args.statsfile,sep='\t', low_memory=False)
        #print(f'\n****  Creating VCF template  ****')
        print(f'inputfile (summary stats): {args.statsfile}')

        duplicates = datain[datain['CPTID'].duplicated(keep=False)]
        if duplicates.shape[0] > 0:
            print(f'Warning: {duplicates.shape[0]} duplicate CPTID variants detected. Removing..')
            duplicates_file = re.split(r'.gz$', args.statsfile)[0] + '.duplicates.vcf.gz'
            duplicates_file = os.path.basename(duplicates_file)
            duplicates_file = args.outdir+"/"+duplicates_file
            print(f'Warning: Writing duplicate CPTIDs to {duplicates_file}')
            template_vcf(duplicates_file)
            append_vcf(duplicates_file, duplicates, noSortFile=args.do_not_sort)
            
            datain = datain[~datain['CPTID'].duplicated(keep=False)]

        template_vcf(outFilename)

        ## populate template
        append_vcf(outFilename, statsData=datain, noSortFile=args.do_not_sort)
    else:
        print(f'Warning: A file named {outFilename} already exists.')
        print('Program halted')
        sys.exit(NameError)

if __name__ == '__main__':
    main(args)
