# go2-scripts
scripts for pre-processing data files prior to meta-analysis

## 1. `filename_checker.py`
Checks whether file name follows the GO2 naming convention prior to analysis.

For detailed information, refer to the [filename_checker documentation wiki page](https://github.com/hmgu-itg/go2-scripts/wiki/filename_checker-documentation).

**Example call**
```
./filename_checker.py \
  TestDat/filepaths.list.txt \
  -v 'config.txt' \
  -o 'output_filename'

### Check usage 
./filename_checker.py -h
```
**output**:
* Returns table with all fields of file name checked and outcome
* output stored as .csvs


## 2. `preQC`

Reads in the table produced by `filename_checker.py` and examines the files that it points to. 
For each sumstats within the table, it ensures that required columns exists, processes columns, adds new columns, and many more and creates a new sumstats file.  
The list of processes `preQC` performs can be found in the [preQC documentation wiki page](https://github.com/hmgu-itg/go2-scripts/wiki/preQC-documentation)

It is recommended to examine the table produced by `filename_checker` and fix any filename error. `preQC` will by default ignore lines with errors unless strict mode is enabled.

**Example call**
Typing `./preQC -h` will print help message

```
./preQC \
  --infile /path/to/filename_checker/files.checkout.csv \
  --outdir /path/to/output/directory \
  --header-file /path/to/header.dict.txt \
  --linear /path/to/linear.software.list.txt \
  [--strict]
```

### header file format example
```
# the columns will be converted to uppercase, so the names below are case insensitive.
# commas are not allowed in column names.
# .* indicate a partial match. Any character following the pattern will be matched as well.
# Use parsimoniously: if multiple columns match, an error will be thrown.
# ^ indicate a beginning of string match. E.g. info will match ^info.* but not frequentist_add_info

POS ^POS.*,POSB37
CHR CHR.*
OR ALL_OR
P frequentist_add_pvalue,PVAL.*
BETA frequentist_add_beta_1,beta.*
INFO ^info.*
EA effect_allele,alleleB
NEA other_allele,alleleA
SE ^se.*
OR_L95 all_OR_lower,ALL_OR_lower,ALL_OR_Lower,L95
OR_U95 all_OR_upper,ALL_OR_upper,ALL_OR_Upper,U95
```


## 3. `runliftnorm`

`runliftnorm` accepts a summary statistic file to create a normalised sumstats file. The purpose of this script is to harmonise different sumstat files to simplify downstream analysis such as meta-analysis.

For detailed information of the process, refer to the [runliftnorm documentation wiki page](https://github.com/hmgu-itg/go2-scripts/wiki/runliftnorm-documentation).

**Example run**
```bash
./runliftnorm \
  -f /path/to/EasyQCed/CLEANED.cohort_37_sumstats.gz \
  -C ./config_liftnorm.txt \
  -O /path/to/output-directory \
  --overwrite
```

Typing `./runliftnorm -h` will print help message.

**Input file format**
* Input file must be a summary statistics file processed by EasyQC
* Filename should end with `.gz`
* Filename must contain either `_37_` or `_38_` depending on variant position's build
* File content should be tab-delimitted
