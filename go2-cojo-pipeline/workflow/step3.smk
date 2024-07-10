import os
import operator
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


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



def input_subset_bfile_for_cojo(w):
    d = pd.read_csv(f'{w.output}/{w.ethnicity}/across-pheno/{OUTDIR}/to_run_cond.csv')
    chrom = d.loc[d['pairID']==w.pairID, 'group_chrom'].iloc[0]
    return bfile_params(w.ethnicity, chrom), bfile(w.ethnicity, chrom)

rule subset_bfile_for_cojo:
    input:
        to_extract=f'{{output}}/{{ethnicity}}/across-pheno/{OUTDIR}/conditional/{{pairID}}/to_extract', # rules.get_to_extract.output[0],
        bfile=lambda w: input_subset_bfile_for_cojo(w)[1]
    params:
        bfile=lambda w: input_subset_bfile_for_cojo(w)[0],
        minbfile=f'{{output}}/{{ethnicity}}/across-pheno/{OUTDIR}/conditional/{{pairID}}/minbfile'
    resources:
        mem_mb=10000
    output:
        multiext(f'{{output}}/{{ethnicity}}/across-pheno/{OUTDIR}/conditional/{{pairID}}/minbfile', '.bed', '.bim', '.fam')
    shell: '''
        plink \
            --bfile {params.bfile} \
            --extract {input.to_extract} \
            --out {params.minbfile} \
            --make-bed \
            --threads {threads} \
            --memory {resources.mem_mb} \
            --silent
    '''


def get_left_pheno(w):
    d = pd.read_csv(f'{w.output}/{w.ethnicity}/across-pheno/{OUTDIR}/to_run_cond.csv')
    left_pheno = d.loc[d['pairID']==w.pairID, 'leftPheno'].iloc[0]
    return left_pheno

def info_to_extract(w):
    d = pd.read_csv(f'{w.output}/{w.ethnicity}/across-pheno/{OUTDIR}/to_run.csv')
    row = d.loc[d['pairID']==w.pairID].iloc[0]
    return row['leftVar'], row['leftPheno'], row['rightVar'], row['rightPheno'], row['group_chrom'], row['group_start'], row['group_end']

rule run_conditional:
    input:
        cojofile=lambda w: f'{w.output}/{w.ethnicity}/{get_left_pheno(w)}/cojofile.ma', # rules.make_cojofile.output[0].format(phenotype=get_left_pheno(w)),
        minbfile=rules.subset_bfile_for_cojo.output,
        to_extract=f'{{output}}/{{ethnicity}}/across-pheno/{OUTDIR}/conditional/{{pairID}}/to_extract' # rules.get_to_extract.output[0]
    params:
        other_var=lambda w: info_to_extract(w)[2],
        minbfile=rules.subset_bfile_for_cojo.params.minbfile,
        out=f'{{output}}/{{ethnicity}}/across-pheno/{OUTDIR}/conditional/{{pairID}}/cond'
    output:
        minss=f'{{output}}/{{ethnicity}}/across-pheno/{OUTDIR}/conditional/{{pairID}}/minss',
        to_cond=f'{{output}}/{{ethnicity}}/across-pheno/{OUTDIR}/conditional/{{pairID}}/to_condition',
        cma=f'{{output}}/{{ethnicity}}/across-pheno/{OUTDIR}/conditional/{{pairID}}/cond.cma.cojo',
        cma2=f'{{output}}/{{ethnicity}}/across-pheno/{OUTDIR}/conditional/{{pairID}}/cond.cma.cojo2'
    shell: '''
        fgrep -w -f <(cat <(echo SNP) {input.to_extract}) {input.cojofile} > {output.minss}
        echo {params.other_var} > {output.to_cond}

        this_var=$(grep -v -w {params.other_var} {input.to_extract})
        gcta64 \
			--bfile {params.minbfile} \
			--cojo-file {output.minss} \
			--cojo-collinear 0.99 \
            --threads {threads} \
			--cojo-cond {output.to_cond} \
			--out {params.out} 2>&1 > /dev/null || true

        if [ ! -f "{params.out}.freq.badsnps" ] ; then
            # If the run was successful
            paste <(echo -e "pairID\n{wildcards.pairID}") {output.cma} <(echo -e "cond-flag\n") > {output.cma2}
        else
            # If the run failed because of big freq difference of other_var between cojofile and bfile
            echo -e "Chr SNP bp refA freq b se p n freq_geno bC bC_se pC" | tr ' ' '\t' > {output.cma}
            echo -e "pairID Chr SNP bp refA freq b se p n freq_geno bC bC_se pC cond-flag" | tr ' ' '\t' > {output.cma2}
            grep -w $this_var {output.minss} | awk -F' ' '{{print "{wildcards.pairID}","",$1,"",$2,$4,$5,$6,$7,$8,"","","","","cojo-cond-failed"}}' | tr ' ' '\t' >> {output.cma2}
        fi
    '''

def input_run_all_conditional(w):
    d = pd.read_csv(f'{w.output}/{w.ethnicity}/across-pheno/{OUTDIR}/to_run_cond.csv')
    return expand(f'{w.output}/{w.ethnicity}/across-pheno/{OUTDIR}/conditional/{{pairID}}/cond.cma.cojo2', pairID=set(d.loc[d['to-condition'].notna(), 'pairID']))
    

rule run_all_conditional:
    input:
        input_run_all_conditional
    output:
        f'{{output}}/{{ethnicity}}/across-pheno/{OUTDIR}/conditional/all.cma.cojo'
    run:
        all_cma = pd.concat([pd.read_csv(f, sep = '\t') for f in input])
        all_cma.to_csv(output[0], sep = '\t', index = False)
