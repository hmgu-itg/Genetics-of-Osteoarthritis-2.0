import pandas as pd

input_f = snakemake.input[0]
phenotype = snakemake.wildcards['phenotype']
output_f = snakemake.output[0]

sigs = pd.read_csv(input_f, sep = ' ')\
         .sort_values(['CHR', 'BP'])\
         .reset_index(drop=True)

sigs['overlap'] = ((sigs['CHR'].rolling(window=2, min_periods=1).max()==sigs['CHR'].rolling(window=2, min_periods=1).min())) \
                  & ((sigs['BP'].rolling(window=2, min_periods=1).max()-sigs['BP'].rolling(window=2, min_periods=1).min())<=1_000_000)

group = 1
out = list()
for overlap in sigs['overlap']:
    if overlap is True:
        out.append(group)
    else:
        group+=1
        out.append(group)
sigs['groups'] = out



regions = list()
for group in set(sigs['groups']):
    subset = sigs[sigs['groups']==group]
    chrom = subset['CHR'].iloc[0]
    start, end = subset['BP'].min() - 1_000_000, subset['BP'].max() + 1_000_000
    start = 0 if start < 0 else start

    indeps = ','.join(subset['SNP'])

    regions.append([chrom, start, end, indeps])

regions = pd.DataFrame(regions, columns = ['chrom', 'start', 'end', 'indeps'])
regions.insert(0, 'phenotype', phenotype)

regions.to_csv(output_f, index=False)