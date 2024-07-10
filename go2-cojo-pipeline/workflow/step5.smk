'''
To run conditional against known-signals.
Example
-------
$ snakemake \
  --snakefile workflow/step5 \
  --cores 100 \
  --config \
    signals='GO2.20230201.txt' \
    positive_controls='GO2_full_positiveControls_1Feb2023.nospaces.txt' \
    window_span=1000000 \
  --use-singularity \
  --singularity-args="-B /" \
  --keep-going \
  output/ALL/pairwise/all.cma.cojo -n
'''
container: 'worker_3.2.3'

import pandas as pd

metal = ''
BFILE_ALL = ''
BFILE_ETHNICITY = ''


COJOFILE = '/file/to/output/{ethnicity}/{phenotype}/cojofile.ma'


SIGNALS = config['signals']
POSITIVE_CONTROLS = config['positive_controls']
WINDOW_SPAN = config['window_span']


checkpoint create_runlist:
    input:
        signals=SIGNALS,
        positive_controls=POSITIVE_CONTROLS
    output:
        '{output}/{ethnicity}/pairwise_runlist.csv'
    run:
        signals = pd.read_csv(input.signals, sep = '\t')
        subset = signals[['groups', 'subgroup', 'phenotype', 'ID', 'chrom', 'pos']]

        KNOWN = pd.read_csv(input.positive_controls, header=None)
        KNOWN.columns = ['known']
        split_string = KNOWN['known'].str.split('(:|_)')
        KNOWN['chrom'] = split_string.str[0].astype(int)
        KNOWN['pos'] = split_string.str[2].astype(int)

        data = list()
        for _, row in subset.iterrows():
            chrom = row['chrom']
            pos = row['pos']
            start = pos - WINDOW_SPAN
            end = pos + WINDOW_SPAN
            
            known_variants = KNOWN.loc[(KNOWN['chrom']==chrom) & (start < KNOWN['pos']) & (KNOWN['pos'] < end), 'known']
            
            
            for known_var in known_variants:
                data.append([*row[:4], known_var])

        out_df = pd.DataFrame(data, columns = [*signals.columns[:4], 'known_var'])
        out_df.to_csv(output[0], index=False)


rule check_known_exists:
    input:
        rules.create_runlist.output[0],
        COJOFILE
    output:
        '{output}/{ethnicity}/checklist/{phenotype}_checklist.csv'
    run:
        runlist = pd.read_csv(input[0])
        cojofile = pd.read_csv(input[1], sep = ' ')
        subset_runlist = runlist[(runlist['phenotype']==wildcards.phenotype) & (runlist['known_var'].isin(cojofile['SNP']))]

        subset_runlist.to_csv(output[0], index=False)



def check_all_known_variants_input(w):
    input_f = checkpoints.create_runlist.get(output=w.output, ethnicity=w.ethnicity).output[0]
    runlist = pd.read_csv(input_f)
    runlist_subset = runlist[runlist['ID']!=runlist['known_var']]
    return expand(rules.check_known_exists.output[0], output=w.output, ethnicity=w.ethnicity, phenotype=set(runlist_subset['phenotype']))


checkpoint check_all_known_variants:
    input:
        check_all_known_variants_input
    output:
        '{output}/{ethnicity}/pairwise_runlist.checklist.csv'
    run:
        out_df = pd.concat([pd.read_csv(f) for f in input])\
                   .sort_values(['groups', 'subgroup'])
        out_df.to_csv(output[0], index=False)



rule create_extract:
    input:
        rules.check_all_known_variants.output[0]
    output:
        '{output}/{ethnicity}/pairwise/{group}/{subgroup}/{phenotype}/{ID}/{known_var}/extract.txt'
    shell: '''
        echo {wildcards.ID} > {output}
        echo {wildcards.known_var} >> {output}
    '''


def bfile(ethnicity, chrom):
    if ethnicity=='ALL':
        return multiext(BFILE_ALL.format(chrom=chrom),
                        '.bed', '.bim', '.fam')
    else:
        return multiext(BFILE_ETHNICITY.format(ethnicity=ethnicity, chrom=chrom),
                        '.bed', '.bim', '.fam')

def bfile_params(ethnicity, chrom):
    if ethnicity=='ALL':
        return BFILE_ALL.format(chrom=chrom)
    else:
        return BFILE_ETHNICITY.format(ethnicity=ethnicity, chrom=chrom)

def input_subset_bfile_for_cojo(w):
    chrom = w.known_var.split(':')[0]
    return bfile_params(w.ethnicity, chrom), bfile(w.ethnicity, chrom)

rule subset_bfile_for_cojo:
    input:
        to_extract=rules.create_extract.output[0],
        bfile=lambda w: input_subset_bfile_for_cojo(w)[1]
    params:
        bfile=lambda w: input_subset_bfile_for_cojo(w)[0],
        minbfile='{output}/{ethnicity}/pairwise/{group}/{subgroup}/{phenotype}/{ID}/{known_var}/minbfile'
    resources:
        mem_mb=10000
    threads:
        50
    output:
        multiext('{output}/{ethnicity}/pairwise/{group}/{subgroup}/{phenotype}/{ID}/{known_var}/minbfile', '.bed', '.bim', '.fam')
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


def select_cojofile(ethnicity, phenotype):
    return COJOFILE.format(ethnicity=ethnicity, phenotype=phenotype)

rule run_conditional:
    input:
        cojofile=lambda w: select_cojofile(w.ethnicity, w.phenotype),
        minbfile=rules.subset_bfile_for_cojo.output,
        to_extract=rules.create_extract.output[0]
    threads:
        workflow.cores * 0.1
    params:
        minbfile=rules.subset_bfile_for_cojo.params.minbfile,
        out='{output}/{ethnicity}/pairwise/{group}/{subgroup}/{phenotype}/{ID}/{known_var}/cond'
    output:
        minss='{output}/{ethnicity}/pairwise/{group}/{subgroup}/{phenotype}/{ID}/{known_var}/minss',
        to_cond='{output}/{ethnicity}/pairwise/{group}/{subgroup}/{phenotype}/{ID}/{known_var}/to_condition',
        cma='{output}/{ethnicity}/pairwise/{group}/{subgroup}/{phenotype}/{ID}/{known_var}/cond.cma.cojo',
        cma2='{output}/{ethnicity}/pairwise/{group}/{subgroup}/{phenotype}/{ID}/{known_var}/cond.cma.cojo2'
    shell: '''
        fgrep -w -f <(cat <(echo SNP) {input.to_extract}) {input.cojofile} > {output.minss}
        echo {wildcards.known_var} > {output.to_cond}

        this_var=$(grep -v -w {wildcards.known_var} {input.to_extract})
        known_var_not_exists=$(grep '{wildcards.known_var}' {params.minbfile}.bim || true)
        if [ -z "$known_var_not_exists" ] ; then
            # If the run failed because the known_variant doesn't exist in the reference panel genotype data
            echo -e "Chr SNP bp refA freq b se p n freq_geno bC bC_se pC" | tr ' ' '\t' > {output.cma}
            echo -e "groups subgroup phenotype ID known_var Chr SNP bp refA freq b se p n freq_geno bC bC_se pC cond-flag" | tr ' ' '\t' > {output.cma2}
            grep -w $this_var {output.minss} | awk -F' ' '{{print "{wildcards.group}","{wildcards.subgroup}","{wildcards.phenotype}",$1,"{wildcards.known_var}","",$1,"",$2,$4,$5,$6,$7,$8,"","","","","known_var-absent"}}' | tr ' ' '\t' >> {output.cma2}
        else
            gcta64 \
                --bfile {params.minbfile} \
                --cojo-file {output.minss} \
                --cojo-collinear 0.99 \
                --threads {threads} \
                --cojo-cond {output.to_cond} \
                --out {params.out} 2>&1 > /dev/null || true

            if [ ! -f "{params.out}.freq.badsnps" ] ; then
                # If the run was successful
                paste <(echo -e "groups\n{wildcards.group}") \
                      <(echo -e "subgroup\n{wildcards.subgroup}") \
                      <(echo -e "phenotype\n{wildcards.phenotype}") \
                      <(echo -e "ID\n$this_var") \
                      <(echo -e "known_var\n{wildcards.known_var}") \
                      {output.cma} <(echo -e "cond-flag\n") > {output.cma2}
            else
                # If the run failed because of big freq difference of other_var between cojofile and bfile
                echo -e "Chr SNP bp refA freq b se p n freq_geno bC bC_se pC" | tr ' ' '\t' > {output.cma}
                echo -e "groups subgroup phenotype ID known_var Chr SNP bp refA freq b se p n freq_geno bC bC_se pC cond-flag" | tr ' ' '\t' > {output.cma2}
                grep -w $this_var {output.minss} | awk -F' ' '{{print "{wildcards.group}","{wildcards.subgroup}","{wildcards.phenotype}",$1,"{wildcards.known_var}","",$1,"",$2,$4,$5,$6,$7,$8,"","","","","cojo-cond-failed"}}' | tr ' ' '\t' >> {output.cma2}
            fi
        fi
    '''    



def input_run_all_conditional(w):
    runlist_f = checkpoints.check_all_known_variants.get(output=w.output, ethnicity=w.ethnicity).output[0]
    runlist = pd.read_csv(runlist_f)
    subset_runlist = runlist[runlist['ID']!=runlist['known_var']]
    return [rules.run_conditional.output.cma2.format(output=w.output,
                                                     ethnicity=w.ethnicity,
                                                     group=row['groups'],
                                                     subgroup=row['subgroup'],
                                                     phenotype=row['phenotype'],
                                                     ID=row['ID'],
                                                     known_var=row['known_var'])
            for _, row in subset_runlist.iterrows()]
    

rule run_all_conditional:
    input:
        input_run_all_conditional
    output:
        '{output}/{ethnicity}/pairwise/all.cma.cojo'
    run:
        all_cma = pd.concat([pd.read_csv(f, sep = '\t') for f in input])
        all_cma.to_csv(output[0], sep = '\t', index = False)

