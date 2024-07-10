#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(argparse))

.BASES = list(A='T', T='A', G='C', C='G')
.complement = function(base) .BASES[[base]]
.strReverse = function(x) sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
flip.strand = function(allele) paste(sapply(strsplit(.strReverse(allele), '')[[1]], .complement), collapse='')

main = function(before.vcf, after.vcf, output.vcf, excluded.vcf) {
    cat('Starting Liftover fix...\n')
    cat('Reading input VCF files\n')
    input_vcf = fread(before.vcf, skip="#CHROM\tPOS")
    output_vcf = fread(after.vcf, skip="#CHROM\tPOS")

    input_subset = input_vcf[,.(ID, REF, ALT)]
    setnames(input_subset, c('ID', 'REF_old', 'ALT_old'))

    cat('Merging..\n')
    output_vcf2 = merge(output_vcf, input_subset, by = 'ID', all.x=T)

    cat('Fixing alleles..\n')
    as_is = output_vcf2[REF!=ALT & REF==REF_old & ALT==ALT_old]
    as_is[,`:=`(FLIPPED=0, Strand=0)]

    output_vcf2 = output_vcf2[!(ID %in% as_is$ID)]

    output_vcf2[, `:=`(
        REF_old_flipped=sapply(REF_old, flip.strand),
        ALT_old_flipped=sapply(ALT_old, flip.strand)
    )]

    as_is_minus = output_vcf2[REF!=ALT & REF==REF_old_flipped & ALT==ALT_old_flipped]
    as_is_minus[,`:=`(FLIPPED=0, Strand=1)]

    output_vcf2 = output_vcf2[!(ID %in% as_is_minus$ID)]

    fixed_ref = output_vcf2[REF==ALT_old & ALT==ALT_old]
    fixed_ref[,`:=`(ALT=REF_old, FLIPPED=1, Strand=0)]

    fixed_ref_minus = output_vcf2[REF==ALT_old_flipped & ALT==ALT_old_flipped]
    fixed_ref_minus[,`:=`(ALT=REF_old_flipped, FLIPPED=1, Strand=1)]



    cannot_determine = intersect(fixed_ref$ID, fixed_ref_minus$ID)
    fixed_ref = fixed_ref[!ID %in% cannot_determine]
    fixed_ref_minus = fixed_ref_minus[!ID %in% cannot_determine]


    as_is_minus[,`:=`(REF_old_flipped=NULL, ALT_old_flipped=NULL)]
    fixed_ref[,`:=`(REF_old_flipped=NULL, ALT_old_flipped=NULL)]
    fixed_ref_minus[,`:=`(REF_old_flipped=NULL, ALT_old_flipped=NULL)]


    to_exclude = output_vcf2[!(ID %in% unique(c(fixed_ref$ID, fixed_ref_minus$ID)))]


    fixed_count = nrow(as_is) + nrow(as_is_minus) + nrow(fixed_ref) + nrow(fixed_ref_minus)
    total_count = nrow(output_vcf)

    cat(paste0(fixed_count, ' / ', total_count, ' variants fixed.\n'))
    cat(paste(nrow(as_is), 'variants with same allele info\n'))
    cat(paste(nrow(as_is_minus), 'variants with strand flipped\n'))
    cat(paste(nrow(fixed_ref), 'variants with allele flipped\n'))
    cat(paste(nrow(fixed_ref_minus), 'variants with strand and allele flipped\n'))
    cat(paste(nrow(to_exclude), 'variants to exclude as we are unable to determine how to fix\n'))

    fixed.combined = rbind(as_is, as_is_minus, fixed_ref, fixed_ref_minus)
    fixed.combined[,ORDER:=paste0(`#CHROM`, ':', POS)]
    fixed.combined = fixed.combined[stringr::str_order(ORDER, numeric=TRUE)]
    fixed.combined[,ORDER:=NULL]

    cat('Creating outputs\n')

    conn = gzfile(after.vcf)
    header=grep("^#", readLines(conn, 5000), value=T)
    if(!grepl( "^#CHROM\tPOS", header[length(header)])){
        cat(paste0("ERROR: Header of flipped file too long. Increase hardcoded header length value.\n"))
        stop(call.=F)
    }
    close(conn)

    conn=gzfile(output.vcf, "w")
    writeLines(header, con=conn)
    close(conn)

    fixed.combined[,`:=`(
        INFO=paste0(INFO,';FLIPPED_lo=',FLIPPED, ';STRAND_lo=', Strand),
        REF_old=NULL,
        ALT_old=NULL,
        FLIPPED=NULL,
        Strand=NULL
    )]
    fwrite(fixed.combined[,.(`#CHROM`, POS, ID, REF, ALT, QUAL, FILTER, INFO)],
           output.vcf, sep = '\t', append=TRUE, compress = 'gzip')


    conn=gzfile(excluded.vcf, "w")
    writeLines(header, con=conn)
    close(conn)

    to_exclude[,`:=`(
        REF_old=NULL,
        ALT_old=NULL,
        REF_old_flipped=NULL,
        ALT_old_flipped=NULL
    )]
    fwrite(to_exclude[,.(`#CHROM`, POS, ID, REF, ALT, QUAL, FILTER, INFO)],
           excluded.vcf, sep = '\t', append=TRUE, compress = 'gzip')

    cat('Done liftover fix\n')
}



## create parser object
parser = ArgumentParser(description="Fix allele info of lifted over VCF file")

## set options
parser$add_argument("-b", "--before-vcf", required=T,
    help="(REQUIRED): VCF file used as input to CrossMap")
parser$add_argument("-a", "--after-vcf", required=T,
    help="(REQUIRED): VCF file outputted by CrossMap")
parser$add_argument("-o", "--output-vcf", required=T,
    help="(REQUIRED): Output VCF file which contains the fixed allele information")
parser$add_argument("-e", "--excluded-vcf", required=T,
    help="(REQUIRED): Output VCF file which contains variants which could not be fixed")

# get command line options, otherwise if options not found on command line, set defaults
# if help option encountered print help and exit,
if(length(commandArgs(T))==0){
  parser$print_help()
}
args = parser$parse_args()



main(before.vcf = args$before_vcf,
     after.vcf = args$after_vcf,
     output.vcf = args$output_vcf,
     excluded.vcf = args$excluded_vcf)