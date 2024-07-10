#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("R.utils"))

## create parser object
parser = ArgumentParser(description=" 1. Re-assigns REF and ALT where flipping occured before normalization and  2.Merges sumstats file from EasyQC to the normalized output")

## set options
parser$add_argument("-s", "--sumstatsfile", required=T,
    help="(REQUIRED): A tab seperated summary stats file from EasyQC")
parser$add_argument("-n", "--normfile", required=T,
    help="(REQUIRED): A tab seperated file of VCF fixed columns from running bcftools norm with realignment")
parser$add_argument("-d", "--no-diagnostic-file", action='store_false', default=TRUE,
    help="Do not generate a diagnostic file with old EA/NEA and REF/ALT alleles.")
parser$add_argument("-o", "--outdir", required=F, default=getwd(),
        help="Directory in which to generate the processed output")

parser$add_argument("--unmapped-file",
    help="Path to VCF file outputted by CrossMap which contain unmapped variants")

parser$add_argument("--non-atgc-file",
    help="Path to VCF file which contain non-ATGC allele variants")

parser$add_argument("--invalid-chrom-file",
    help="Path to VCF file which contain variants with invalid chromosome value")
    
parser$add_argument("--duplicate-ids",
    help="Path to VCF file which contain variants with duplicate CPTIDs")

parser$add_argument("--liftover-excluded",
    help="Path to VCF file which contain variants excluded by fix_liftover.R")

# get command line options, otherwise if options not found on command line, set defaults
# if help option encountered print help and exit,
if(length(commandArgs(T))==0){
  parser$print_help()
}
args = parser$parse_args()

##################################################################################################################################
## Functions
##################################################################################################################################

inform=function(...){
    cat(paste0("Info: ",...,"\n"))
    flush.console()
}

error=function(...){
    cat(paste0("ERROR: ",...,"\n"))
    stop(call.=F)
}

fixallele = function(sumStats, normDat){
    #* sets EA to ALT and NEA to REF where no flipping occured
    #* sets EA to REF and NEA to ALT where flipping occured

    #* sumStats -> original summary stats from EasyQC used to create VCF
    #* normDat -> fixed columns from bcftools norm with realignment
    #* getDiag -> generate diagnostic file
    #============================================================================

    ## Prepare normed data
    normDat[,CHRNUM:=as.numeric(sub('chr', '', `#CHROM`))]
    normDat[,newid:=paste0(CHRNUM, ":", POS,"_", REF, "_", ALT)]
    
    if (grepl('FLIPPED_lo', normDat$INFO[1])) {
        normDat[,c("AF", "FLIPPED_lo", "STRAND_lo", "FLIPPED", "OLD"):=tstrsplit(INFO, ";")]
        normDat[,`:=`(
            FLIPPED_lo=as.logical(as.integer(tstrsplit(FLIPPED_lo, "=")[[2]])),
            STRAND_lo=as.logical(as.integer(tstrsplit(STRAND_lo, "=")[[2]])),
            FLIPPED=as.logical(as.integer(tstrsplit(FLIPPED, "=")[[2]]))
        )]
    } else {
        normDat[,c("AF", "FLIPPED", "OLD"):=tstrsplit(INFO, ";")]
        normDat[,`:=`(
            FLIPPED=as.logical(as.integer(tstrsplit(FLIPPED, "=")[[2]])),
            FLIPPED_lo=FALSE,
            STRAND_lo=FALSE)]
    }
    
    normDat[,c("#CHROM", "QUAL", "FILTER", "INFO", "AF", "OLD"):=NULL]
    setnames(normDat, c("CHRNUM", "POS"), c("CHR_new", "POS_new"))

    ## Merge data stes by CPTID
    #inform('Merging sumstats data with normalized data ....')
    normDat=merge(sumStats,normDat, by.x='CPTID', by.y="ID")
    gc()
    ## Using FLIPPED  column, reassign EA/NEA
    #inform('Checking and re-assigning flipped sites .....')
    # allows to reassing in place rather than bothering with col order
    normDat[,c("EA_old", "NEA_old"):=list(EA, NEA)]
    normDat[,c("EA", "NEA"):=list(ALT, REF)]
    normDat[(FLIPPED_lo==T & FLIPPED==F)
            | (FLIPPED_lo==F & FLIPPED==T), c("EA", "NEA"):=list(REF, ALT)]

    # Update ID
    normDat[,c("CPTID_old"):=CPTID]
    normDat[,c("CPTID", "newid"):=list(newid, NULL)]

    normDat[,c("CHR_old", "POS_old"):=list(CHR, POS)]
    normDat[,c("CHR", "POS", "CHR_new", "POS_new"):=list(CHR_new, POS_new, NULL, NULL)]

    return(normDat)

}
#Fxn Test
#test_out=fixallele(sumStats=dat_sumstats, normDat=dat_norm, getDiag=TRUE)
## Write final data
#fwrite(test_out, file = paste(outDir, outFilname, sep="/"), sep="\t", quote=F, compress = 'gzip')
#==================================================================================================================================
###################################################################################################################################



read.invalid_chr = function(dat_sumstats, inputfile){
    invalid_chr = fread(inputfile, skip="#CHROM\tPOS")
    invalid_chr[,INFO:=NULL]
    invalid_chr[,ID:=as.character(ID)]
    setnames(invalid_chr, c("#CHROM", "POS"), c("CHR_new", "POS_new"))

    invalid_chr.m = merge(dat_sumstats, invalid_chr, by.x='CPTID', by.y="ID")
    setnames(invalid_chr.m, c('CPTID', 'CHR', 'POS', 'CHR_new', 'POS_new'), c('CPTID_old', 'CHR_old', 'POS_old', 'CHR', 'POS'))
    invalid_chr.m[,`:=`(
                CPTID=NA,
                EA_old=EA,
                NEA_old=NEA,
                FLIPPED=NA,
                FLIPPED_lo=NA,
                STRAND_lo=NA,
                Excluded='Invalid liftover chromosome'
        )]
    return(invalid_chr.m[,c("CPTID_old", "CPTID", "CHR_old", "POS_old", "CHR", "POS", "EA_old", "NEA_old", "EA", "NEA", "FLIPPED", "FLIPPED_lo", "STRAND_lo", "Excluded")])
}

read.unmapped = function(dat_sumstats, inputfile){
    unmapped = fread(inputfile, skip="#CHROM\tPOS")
    setnames(unmapped, c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'Excluded'))
    unmapped[,INFO:=NULL]
    unmapped[,ID:=as.character(ID)]
    setnames(unmapped, c("#CHROM", "POS"), c("CHR_new", "POS_new"))

    unmapped.m = merge(dat_sumstats, unmapped, by.x='CPTID', by.y="ID")
    setnames(unmapped.m, c('CPTID', 'CHR', 'POS', 'CHR_new', 'POS_new'), c('CPTID_old', 'CHR_old', 'POS_old', 'CHR', 'POS'))
    unmapped.m[,`:=`(
                CPTID=NA,
                EA_old=EA,
                NEA_old=NEA,
                FLIPPED=NA,
                FLIPPED_lo=NA,
                STRAND_lo=NA
        )]

    unmapped.m[,Excluded:=sub('Fail', 'Liftover failed ', Excluded)]
    return(unmapped.m[,c("CPTID_old", "CPTID", "CHR_old", "POS_old", "CHR", "POS", "EA_old", "NEA_old", "EA", "NEA", "FLIPPED", "FLIPPED_lo", "STRAND_lo", "Excluded")])
}

read.non_atgc = function(dat_sumstats, inputfile){
    non_atgc = fread(inputfile, skip="#CHROM\tPOS")
    non_atgc[,INFO:=NULL]
    non_atgc[,ID:=as.character(ID)]
    setnames(non_atgc, c("#CHROM", "POS"), c("CHR_new", "POS_new"))
    non_atgc.m = merge(dat_sumstats, non_atgc, by.x='CPTID', by.y="ID")
    setnames(non_atgc.m, c('CPTID', 'CHR', 'POS', 'CHR_new', 'POS_new'), c('CPTID_old', 'CHR_old', 'POS_old', 'CHR', 'POS'))

    non_atgc.m[,`:=`(
                CPTID=NA,
                EA_old=EA,
                NEA_old=NEA,
                FLIPPED=NA,
                FLIPPED_lo=NA,
                STRAND_lo=NA,
                Excluded='non-ATGC allele'
        )]
    return(non_atgc.m[,c("CPTID_old", "CPTID", "CHR_old", "POS_old", "CHR", "POS", "EA_old", "NEA_old", "EA", "NEA", "FLIPPED", "FLIPPED_lo", "STRAND_lo", "Excluded")])
}

read.duplicates = function(dat_sumstats, inputfile){
    duplicates = fread(inputfile, skip="#CHROM\tPOS")
    duplicates[,INFO:=NULL]
    duplicates[,ID:=as.character(ID)]
    setnames(duplicates, c("#CHROM", "POS"), c("CHR_new", "POS_new"))
    duplicates.m = merge(dat_sumstats, duplicates, by.x='CPTID', by.y="ID")
    setnames(duplicates.m, c('CPTID', 'CHR', 'POS', 'CHR_new', 'POS_new'), c('CPTID_old', 'CHR_old', 'POS_old', 'CHR', 'POS'))

    duplicates.m[,`:=`(
                CPTID=NA,
                EA_old=EA,
                NEA_old=NEA,
                FLIPPED=NA,
                FLIPPED_lo=NA,
                STRAND_lo=NA,
                Excluded='duplicated CPTID'
        )]
    return(duplicates.m[,c("CPTID_old", "CPTID", "CHR_old", "POS_old", "CHR", "POS", "EA_old", "NEA_old", "EA", "NEA", "FLIPPED", "FLIPPED_lo", "STRAND_lo", "Excluded")])
}

read.liftover_excluded = function(dat_sumstats, inputfile){
    liftover_excluded = fread(inputfile, skip="#CHROM\tPOS")
    liftover_excluded[,INFO:=NULL]
    liftover_excluded[,ID:=as.character(ID)]
    setnames(liftover_excluded, c("#CHROM", "POS"), c("CHR_new", "POS_new"))
    liftover_excluded.m = merge(dat_sumstats, liftover_excluded, by.x='CPTID', by.y="ID")
    setnames(liftover_excluded.m, c('CPTID', 'CHR', 'POS', 'CHR_new', 'POS_new'), c('CPTID_old', 'CHR_old', 'POS_old', 'CHR', 'POS'))

    liftover_excluded.m[,`:=`(
                CPTID=NA,
                EA_old=EA,
                NEA_old=NEA,
                FLIPPED=NA,
                FLIPPED_lo=NA,
                STRAND_lo=NA,
                Excluded='Could not fix allele of liftedover variant'
        )]
    return(liftover_excluded.m[,c("CPTID_old", "CPTID", "CHR_old", "POS_old", "CHR", "POS", "EA_old", "NEA_old", "EA", "NEA", "FLIPPED", "FLIPPED_lo", "STRAND_lo", "Excluded")])
}

read.to_exclude = function(dat_sumstats, to_exclude_dt){
    to_exclude_dt[,INFO:=NULL]
    to_exclude_dt[,ID:=as.character(ID)]
    setnames(to_exclude_dt, c("#CHROM", "POS"), c("CHR_new", "POS_new"))
    to_exclude_dt.m = merge(dat_sumstats, to_exclude_dt, by.x='CPTID', by.y="ID")
    setnames(to_exclude_dt.m, c('CPTID', 'CHR', 'POS', 'CHR_new', 'POS_new'), c('CPTID_old', 'CHR_old', 'POS_old', 'CHR', 'POS'))

    to_exclude_dt.m[,`:=`(
                CPTID=NA,
                EA_old=EA,
                NEA_old=NEA,
                FLIPPED=NA,
                FLIPPED_lo=NA,
                STRAND_lo=NA,
                Excluded='Special character allele assigned by leftalignment'
        )]
    return(to_exclude_dt.m[,c("CPTID_old", "CPTID", "CHR_old", "POS_old", "CHR", "POS", "EA_old", "NEA_old", "EA", "NEA", "FLIPPED", "FLIPPED_lo", "STRAND_lo", "Excluded")])
}


###############################
## Main
###############################

## check/create directory
outDir=(args$outdir)
if(file.exists(outDir)){
  if(dir.exists(outDir)) {
    inform("The specified output directory ", args$outdir, " exists, outputs will be overwritten")
  }else error("The specified output directory ", args$outdir, " exists and is a file")
} else{
  inform("Creating output directory ", args$outdir)
  dir.create(args$outdir)
}
outDir=normalizePath(outDir)

## get output file name
if(file.exists(args$sumstatsfile)){
    outFilename=paste0(tools::file_path_sans_ext(basename(args$sumstatsfile)), '.norm.txt.gz')
}else{
    error('File ', args$sumstatsfile, ' does not exist')
}

## test existence of input files
if(!file.exists(args$normfile)){
  error("File ", args$normfile, " does not exist")
}

## diagnostic file
if(args$no_diagnostic_file){ #logic is inverted because store_false
  inform('Diagnostic file will be generated')
}else{
  inform('No diagnostic file will be generated')
}

## Read in data files and set column names
#inform('reading normalized data file .....')
dat_norm= fread(args$normfile, skip="#CHROM\tPOS")
if(ncol(dat_norm) < 8){
      error('Normalized data file has fewer than 8 columns.')
    }

#inform('reading in sumstats data file .....')
dat_sumstats=fread(args$sumstatsfile)

# Convert Chromosome X to 23 to prevent coercion of X to NA
dat_sumstats[CHR=='X', CHR:='23']
dat_norm[`#CHROM`=='chrX', `#CHROM`:='chr23']

# Leftaligned alleles are saved as lowercase
to_exclude = dat_norm[!grepl('^[ATGCNatgcn]+$', REF) | !grepl('^[ATGCNatgcn]+$', ALT)]
nrow_before = nrow(dat_norm)
dat_norm = dat_norm[grepl('^[ATGCNatgcn]+$', REF) & grepl('^[ATGCNatgcn]+$', ALT)]
stopifnot(nrow_before == (nrow(to_exclude) + nrow(dat_norm)))

dat_norm[,`:=`(REF=toupper(REF), ALT=toupper(ALT))]

## Fxn Call
outDat=fixallele(sumStats=dat_sumstats, normDat=dat_norm)
outDat = outDat[str_order(CPTID, numeric=TRUE)]

outDat$Excluded = 'FALSE'
outDat[CHR_old!=CHR, Excluded:='Liftover chromosome mismatch']

# Write final outputs
outPath=paste(outDir, outFilename, sep="/")
oldcols=colnames(dat_sumstats)
fwrite(outDat[Excluded=='FALSE',..oldcols], file=outPath, sep="\t", quote=F)

if(args$no_diagnostic_file){
  #inform('writing processed data to: ', outPath)

  diagnostic_data = outDat[,c("CPTID_old", "CPTID", "CHR_old", "POS_old", "CHR", "POS", "EA_old", "NEA_old", "EA", "NEA", "FLIPPED", "FLIPPED_lo", "STRAND_lo", "Excluded")]

  if (!is.null(args$invalid_chrom_file)){
    invalid_chr = read.invalid_chr(dat_sumstats, args$invalid_chrom_file)
    diagnostic_data = rbind(diagnostic_data, invalid_chr)
  }
  if (!is.null(args$unmapped_file)){
    unmapped = read.unmapped(dat_sumstats, args$unmapped_file)
    diagnostic_data = rbind(diagnostic_data, unmapped)
  }
  if (!is.null(args$non_atgc_file)){
    non_atgc = read.non_atgc(dat_sumstats, args$non_atgc_file)
    diagnostic_data = rbind(diagnostic_data, non_atgc)
  }

  if (!is.null(args$duplicate_ids)){
    duplicates = read.duplicates(dat_sumstats, args$duplicate_ids)
    diagnostic_data = rbind(diagnostic_data, duplicates)
  }
  
  if (!is.null(args$liftover_excluded)){
    liftover_excluded = read.liftover_excluded(dat_sumstats, args$liftover_excluded)
    diagnostic_data = rbind(diagnostic_data, liftover_excluded)
  }
  
  if (nrow(to_exclude) > 0){
    to_exclude2 = read.to_exclude(dat_sumstats, to_exclude)
    diagnostic_data = rbind(diagnostic_data, to_exclude2)
  }


  diagnostic_data = diagnostic_data[str_order(CPTID_old, numeric=TRUE)]

  fwrite(diagnostic_data, file=sub(".norm.txt.gz", ".diagnostic.txt.gz", outPath), sep="\t", quote=F)
}
