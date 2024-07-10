import os
import operator
import itertools
from pathlib import Path

import numpy as np
import pandas as pd

container: 'worker_3.2.3'

metal = config['metal']
BFILE = config['bfile']

INPUT_PHENOTYPES_DF = pd.read_csv(config['input_phenotypes'], header=None)
PHENOTYPES = INPUT_PHENOTYPES_DF[0].to_list()
OUTDIR = config['outdir']

def bfile(ethnicity, chrom):
    return multiext(BFILE.format(chrom=chrom), '.bed', '.bim', '.fam')

def bfile_params(ethnicity, chrom):
    return BFILE.format(chrom=chrom)

rule make_cojofile:
    input:
        metal=metal
    output:
        '{output}/{ethnicity}/{phenotype}/cojofile.ma'
    shell: '''
        workflow/scripts/make_cojofile.sh {input.metal} {output}
    '''

def input_get_to_extract(w):
    d = pd.read_csv('{output}/{ethnicity}/across-pheno/{outdir}/to_run.csv'.format(
                    output=w.output, ethnicity=w.ethnicity, outdir=OUTDIR))
    left_pheno = d.loc[d['pairID']==w.pairID, 'leftPheno'].iloc[0]
    sumstats = metal.format(ethnicity=w.ethnicity, phenotype=left_pheno)
    return sumstats

def info_to_extract(w):
    d = pd.read_csv('{output}/{ethnicity}/across-pheno/{outdir}/to_run.csv'.format(
                    output=w.output, ethnicity=w.ethnicity, outdir=OUTDIR))
    row = d.loc[d['pairID']==w.pairID].iloc[0]
    return row['leftVar'], row['leftPheno'], row['rightVar'], row['rightPheno'], row['group_chrom'], row['group_start'], row['group_end']

rule get_to_extract:
    input:
        sumstats=input_get_to_extract,
        bfile=lambda w: bfile(w.ethnicity, info_to_extract(w)[4])
    params:
        this_var=lambda w: info_to_extract(w)[0],
        this_pheno=lambda w: info_to_extract(w)[1],
        other_var=lambda w: info_to_extract(w)[2],
        other_pheno=lambda w: info_to_extract(w)[3],
        chrom=lambda w: info_to_extract(w)[4],
        start=lambda w: info_to_extract(w)[5],
        end=lambda w: info_to_extract(w)[6],
        bfile_in=lambda w: bfile_params(w.ethnicity, info_to_extract(w)[4]),
        bfile_out='{output}/{ethnicity}/across-pheno/{outdir}/conditional/{pairID}/t',
        bfile_out2='{output}/{ethnicity}/across-pheno/{outdir}/conditional/{pairID}/t2'
    output:
        to_extract='{output}/{ethnicity}/across-pheno/{outdir}/conditional/{pairID}/to_extract',
        other_var_status='{output}/{ethnicity}/across-pheno/{outdir}/conditional/{pairID}/other-var-status.txt'
    log:
        '{output}/{ethnicity}/across-pheno/{outdir}/conditional/{pairID}/to_extract.log'
    shell: '''
        touch {log}
        echo {params.this_var} > {output.to_extract}
        # we wrap the zgrep in brackets with "|| true" to prevent rule failure 
        # for cases when the other var doesn't exist in the sumstats
        echo $(zgrep -F -m1 -w {params.other_var} {input.sumstats} || true) | cut -d' ' -f1 >> {output.to_extract}
        echo -e "{wildcards.pairID}\t{params.other_var}\t" > {output.other_var_status}
        
        non_empty_lines=$(cat {output.to_extract} | sed '/^\s*$/d' | wc -l)
        if [[ "$non_empty_lines" -eq 1 ]] ; then
            # If the other variant does not exist in the phenotype sumstats,
            # extract a proxy variant in LD with the other variant
            echo RightVar {params.other_var}, lead for rightPheno {params.other_pheno}, not found in leftPheno {params.this_pheno} sumstats. Linking via LD. > {log}
            plink \
                --bfile {params.bfile_in} \
                --chr {params.chrom} \
                --from-bp {params.start} \
                --to-bp {params.end} \
                --make-bed \
                --threads {threads} \
                --memory {resources.mem_mb} \
                --out {params.bfile_out} >&2
            plink \
                --bfile {params.bfile_out} \
                --r2 \
                --ld-snp {params.other_var} \
                --ld-window 1000000 \
                --ld-window-kb 10000000 \
                --ld-window-r2 0.8 \
                --threads {threads} \
                --memory {resources.mem_mb} \
                --out {params.bfile_out} >&2
            echo Found $(( $(cat {params.bfile_out}.ld | wc -l) - 2 )) variants in LD with {params.other_var} >> {log}
            # index_var=$(fgrep -w -f <(tail -n+2 {params.bfile_out}.ld | awk '{{print $6}}' | fgrep -v {params.other_var}) {input.sumstats} | sort -k7,7h | head -1 | cut -f1)
            fgrep -w -f <(tail -n+2 {params.bfile_out}.ld | awk '{{print $6}}' | fgrep -v {params.other_var}) {input.sumstats} > {params.bfile_out}.tmp || true
            index_var=$(sort -k7,7h {params.bfile_out}.tmp | head -1 | cut -f1)
            echo $index_var selected as the proxy variant for rightVar {params.other_var} >> {log}
            echo $index_var >> {output.to_extract}
            echo -e "{wildcards.pairID}\t$index_var\tother_var proxy" > {output.other_var_status}
        fi

        non_empty_lines=$(cat {output.to_extract} | sed '/^\s*$/d' | wc -l)
        if [[ "$non_empty_lines" -eq 1 ]] ; then
            # If the proxy of the other variant also doesn't exist,
            # then extract the proxy variant in LD with the main variant
            echo Proxy of rightVar {params.other_var} not found. Finding proxy for leftVar {params.this_var} >> {log}
            plink \
                --bfile {params.bfile_out} \
                --r2 \
                --ld-snp {params.this_var} \
                --ld-window 1000000 \
                --ld-window-kb 10000000 \
                --ld-window-r2 0.8 \
                --threads {threads} \
                --memory {resources.mem_mb} \
                --out {params.bfile_out2} >&2
            echo Found $(( $(cat {params.bfile_out2}.ld | wc -l) - 2 )) variants in LD with {params.this_var} >> {log}
            fgrep -w -f <(tail -n+2 {params.bfile_out2}.ld | awk '{{print $6}}' | fgrep -v {params.this_var}) {input.sumstats} > {params.bfile_out2}.tmp || true
            index_var=$(sort -k7,7h {params.bfile_out2}.tmp | head -1 | cut -f1)
            echo $index_var selected as the proxy variant for {params.this_var} >> {log}
            echo $index_var >> {output.to_extract}
            echo -e "{wildcards.pairID}\t$index_var\tthis_var proxy" > {output.other_var_status}
        fi

        non_empty_lines=$(cat {output.to_extract} | sed '/^\s*$/d' | wc -l)
        if [[ "$non_empty_lines" -eq 1 ]] ; then
            echo Could not find proxy for {params.this_var}! >> {log}
            echo -e "{wildcards.pairID}\t\tno proxy" > {output.other_var_status}
        fi
    '''




def input_run_all_conditional(w):
    d = pd.read_csv(f'{w.output}/{w.ethnicity}/across-pheno/{OUTDIR}/to_run.csv')
    return expand('{output}/{ethnicity}/across-pheno/{outdir}/conditional/{pairID}/other-var-status.txt',
                        output=w.output, ethnicity=w.ethnicity, outdir=OUTDIR, pairID=set(d['pairID']))

rule to_run_cond:
    input:
        to_run=f'{{output}}/{{ethnicity}}/across-pheno/{OUTDIR}/to_run.csv',
        var_status=input_run_all_conditional
    output:
        to_run_cond=f'{{output}}/{{ethnicity}}/across-pheno/{OUTDIR}/to_run_cond.csv'
    run:
        d = pd.read_csv(input.to_run)
        var_status = pd.concat([pd.read_csv(f, sep='\t', header=None) for f in input.var_status])
        var_status.columns = ['pairID', 'to-condition', 'cond-flag']

        d2 = d.merge(var_status, how='left')
        d2.to_csv(output.to_run_cond, index=False)

