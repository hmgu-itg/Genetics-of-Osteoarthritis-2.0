import pandas as pd

sumstats = pd.read_csv(snakemake.input[0], sep = '\t')

for chrom in set(sumstats['CHR']):
    subset = sumstats[sumstats['CHR']==chrom]
    subset.to_csv(snakemake.params.output.format(chrom), sep = '\t', index = False)


