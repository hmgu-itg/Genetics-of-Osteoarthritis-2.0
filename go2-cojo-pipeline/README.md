```bash
# Step 1
snakemake \
  --use-singularity \
  --singularity-args="-B /data" \
  --cores 70 \
  'output/across-pheno/AllOA_FingerOA_HipOA_KneeOA/to_run.csv' -n




# Step 2
snakemake \
  --snakefile workflow/step2.smk \
  --use-singularity \
  --singularity-args="-B /data" \
  --cores 70 \
  output/across-pheno/AllOA_FingerOA_HipOA_KneeOA/graph-plots -n
```



```bash
snakemake \
  --use-singularity \
  --singularity-args="-B /data" \
  --keep-going \
  --cores 50 \
  "output/across-pheno/AllOA_FingerOA_HipOA_KneeOA/to_run.csv" -n

snakemake \
  --snakefile workflow/step2.smk \
  --use-singularity \
  --singularity-args="-B /data" \
  --keep-going \
  --cores 50 \
  output/across-pheno/AllOA_FingerOA_HipOA_KneeOA/graph-plots -n

```