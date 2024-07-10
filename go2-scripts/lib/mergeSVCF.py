#!/usr/bin/env python3

##########################################################################################################
# Dependencies
import re
import os
import sys
import gzip
import argparse
import pandas as pd
import numpy as np
##########################################################################################################
# Arguments
parser = argparse.ArgumentParser(
    description='Merge original summary stats file from EasyQC to normalizde VCF output',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-s','--statsfile',help='stats summary file to convert to vcf file (REQUIRED)')
parser.add_argument('-v', '--vcfnormfile',help='fixed columns from normalized VCF output (REQUIRED)')
parser.add_argument('-u', '--unmapfile',help='unmapped variants from liftOver output', nargs="?",
                    default=None)
parser.add_argument('-o', '--outdir', required=True, help="output directory (must exist and be writable)")
args = parser.parse_args()
##########################################################################################################
## Functions
def mergeSVCF(statsfile, vcfnormfile, outdir, unmapfile=None, fileName=None, 
              liftOver=False, unmap=False):
    """
    Merge original summary stats file from EasyQC to normalizde VCF output
    
    input1:  summary stats file
             DataFrame, tab delimitted
    input2:  normalized VCF fixed columns output
             DataFrame, tab delimitted
    input3:  unmapped variants in the case where liftOver was performed
             Dataframe
    """
    # set outfilenames
    if fileName is None:
        outFilename1= outdir + "/" + 'merge_sumstats_vcfnormalized_out'
        outFilename2= outdir + "/" + 'summary_' + outFilename1
    else:
        ### Resolve path, extract filename only and redirect to outdir
        outFilename1=os.path.basename(fileName)
        outFilename2= 'summary_' + outFilename1
        
        outFilename1=outdir + "/" + outFilename1
        outFilename2=outdir + "/" + outFilename2
    
    # Data
    data1=statsfile
    data2=vcfnormfile
    if (unmap==1):
        data3=unmapfile

    ## Get summary of Input records and Output records 
    # ** Add Filename to output
    nSumstats= data1.shape[0]
    nNormalized= data2.shape[0]
    diff_StatsNorm= nSumstats-nNormalized
    
    # check liftOver
    if (liftOver==1) & (unmap==1):
        nUnmap= data3.shape[0]
    else:
        nUnmap=np.nan
    out_summary=pd.DataFrame({'fileName':[outFilename1], 'nSumstats':nSumstats, 'nNormalized':nNormalized,
                              'diff_StatsNorm':diff_StatsNorm, 'nUnmap':nUnmap, 'liftOver':liftOver})

    # Merge data
    # Rename some headers and check duplicates
    data2.rename(columns={'POS':'POS1', 'INFO':'INFO1'}, inplace=True)
    
    print('Info: merging data')
    out1=data1.merge(data2, left_on='CPTID', right_on='ID')
    keepcols= ['CPTID', 'CHR', 'POS', 'EA', 'NEA', 'EAF', 'BETA', 'SE', 'P', 'NCASES',
       'NCONTROLS', 'N', 'Neff', 'EAC', 'INFO', 'REF', 'ALT']
    out2=out1.loc[:, keepcols]
    
    
    # Write output
    print(f'Info: writing merged data to outFilename1: {outFilename1}.gz')
    out2.to_csv(f'{outFilename1}.gz', sep='\t', index=False, compression='gzip')
    
    print(f'Info: writing summary to outFilename2: {outFilename2}')
    out_summary.to_csv(f'{outFilename2}', sep='\t', index=False)
    print('Info: merging completed')
    print('Done \n')
     
#=================================================================================================
def main(args):
    """
    Check passed arguments and merge data
    
    """
    
    ### Create output filename
    if args.statsfile:
        if re.search(r'.gz$', args.statsfile):
            filename=re.split(r'.gz$', args.statsfile)[0] + '.normalized'
        else:
            print('Warning: statsfile name should end with ".gz"')
            sys.exit(NameError)
    else:
        sys.exit("Please provide stats file")
        
    ### Check & read input files
    if args.statsfile:
        print(f'Info: reading statsfile: {args.statsfile}')
        input1=pd.read_csv(args.statsfile, sep='\t', low_memory=False)
    else:
        sys.exit("Please provide stats file")
    
    if args.vcfnormfile:
        print(f'Info: reading vcfnormfile: {args.vcfnormfile}')
        input2=pd.read_csv(args.vcfnormfile, sep='\t', low_memory=False)
    else:
        sys.exit("Please provide vcfnorm file")
        
    ### Check expected headers
    print("Info: checking headers")
    # data1
    header_dat1= input1.columns.to_list()
    header_dat1=[x.upper() for x in input1]
    
    expected_header1= ['CPTID', 'CHR', 'POS', 'EA', 'NEA', 'EAF', 'BETA', 'SE', 
                       'P', 'NCASES','NCONTROLS', 'N', 'NEFF', 'EAC', 'INFO']
    msg = f'Warning: wrong header, missing expected headers in statsfile '
    assert all(x in header_dat1 for x in expected_header1 ), msg

    # data2
    header_dat2= input2.columns.to_list()
    header_dat2=[x.upper() for x in input2]
    
    expected_header2=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    msg = f'Warning: wrong header, missing expected headers in vcfnomfile'
    assert all(x in header_dat2 for x in expected_header2 ), msg 
    print("Info: headers OK")

    ### Check for unmap file and set liftOver
    if args.unmapfile:
        liftOver=True
        unmap=True
        input3=pd.read_csv(args.unmapfile, sep='\t', comment='#',low_memory=False,
                          header=False)
        len_input3=input3.shape[0]
    else:
        liftOver=False
        len_input3=0
        

    ### merge Data
    print(f'\n****  Merge sumstats file and normalized vcf file  ****')
    if (len_input3 > 0) & (liftOver==1):
        print('TRUE: unmap file found')
        mergeSVCF(statsfile=input1, vcfnormfile=input2, outdir=args.outdir,
                  unmapfile=input3,fileName=filename, liftOver=liftOver, unmap=unmap)
    elif len_input3==0:
        mergeSVCF(statsfile=input1, vcfnormfile=input2, outdir=args.outdir,
                  fileName=filename)
        

if __name__ == '__main__':
    main(args)
