"""
2022-08-23 by YC: Notes for later when resuming the development

This snakefile is for setting up the plink binary file for the 1kGP LD reference data for clumping in ethnic groups other than Europeans.
We will be using the UKBB genotype data Kostas prepared for the European cohorts.
For others, it's not ideal, but we'll be using 1kGP (not ideal because 1kGP has small sample size and cojo needs n>4000).


Last thing I was doing was to resolve the multiallelic variants preventing the merge across all chromosome to create a single plink binary fileset. 
Plink 2 provides --rm-dup with different modes, and I was wondering whether it would be best to run with `force-first` or `exclude-all`.
See https://www.cog-genomics.org/plink/2.0/filter#rm_dup
I think force-first mode would keep most SNPs, as a quick check looked like plink2 tends put SNP variant before INDELs, so `force-first` would retain mostly SNPs.
Although I need to confirm how Kostas prepared his UKBB plink binary fileset, whether he kept one or excluded all.
"""


container: '/lustre/groups/itg/shared/containers/worker_3.2.3'

one_kgp_file = '/lustre/groups/itg/shared/referenceData/1kG/20190312_SNV_INDEL/ALL.chr{chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz'

rule setup:
    input:
        expand('output/1kGP/chr{chrom}.bed', chrom = list(range(1, 23)) + ['X'])
    params:
        out='output/1kGP/all'
    output:
        mergelist='output/1kGP/mergelist',
        binary=multiext('output/1kGP/all', '.bed', '.bim', '.fam')
    threads:
        workflow.cores
    shell: '''
        for i in {{1..22}} X ; do echo output/1kGP/chr$i ; done > {output.mergelist}

        plink \
            --merge-list {output.mergelist} \
            --make-bed \
            --out {params.out} \
            --memory {resources.mem_mb} \
            --threads {threads} \
            --silent
    '''

rule convert_vcf2bed:
    input:
        one_kgp_file
    params:
        out='output/1kGP/chr{chrom}'
    threads:
        workflow.cores
    output:
        bfile=multiext('output/1kGP/chr{chrom}', '.bed', '.bim', '.fam')
    shell: '''
        plink2 \
            --vcf {input[0]} \
            --vcf-half-call missing \
            --vcf-require-gt \
            --make-bed \
            --new-id-max-allele-len 100 missing \
            --set-all-var-ids chr@:# \
            --maj-ref force \
            --memory {resources.mem_mb} \
            --threads {threads} \
            --out {params.out} \
            --silent
        '''


rule deduplicate:
    input:
        rules.convert_vcf2bed.output.bfile
    params:
        input=rules.convert_vcf2bed.params.out,
        out='output/1kGP/chr{chrom}-dedup'
    output:
        bfile=multiext('output/1kGP/chr{chrom}-dedup', '.bed', '.bim', '.fam')
    shell: '''
        plink2 \
            --bfile {params.input} \
            --make-bed \
            --rm-dup exclude-all \
            --memory {resources.mem_mb} \
            --threads {threads} \
            --out {params.out} \
            --silent
    '''