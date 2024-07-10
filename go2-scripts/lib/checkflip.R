#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("R.utils"))

## create parser object
parser = ArgumentParser(description="Checks and gets info where REF allele was flipped by bcftools norm")

## set options
parser$add_argument("-p", "--preflipfile", required=T,
    help="(REQUIRED): A VCF before alleles were flipped.")
parser$add_argument("-r", "--setreffile", required=T,
    help="(REQUIRED): A VCF after alleles were flipped.")
parser$add_argument("-o", "--outdir", required=F, default=getwd(),
        help="Directory in which to generate the processed output")

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

checkflip = function(preflipDat, setrefDat){
    #* checks if ALT and REF was flipped before norm step 1
    #* sets EA to ALT and NEA to REF where no flipping occured
    #* sets EA to REF and NEA to ALT where flipping occured

    #* setref -> fixed columns from bcftools norm with no realignment
    #* preflipDat-> fixed columns data pre bcftools norm step1
    #====================================================================

    # Set headers
    if(ncol(preflipDat) < 8 | ncol(setrefDat) < 8){
        print(colnames(preflipDat))
        print(colnames(setrefDat))
        error('Unexpected number of columns in input files.')
    }
    n1=nrow(preflipDat)
    n2=nrow(setrefDat)
    # preflipDat=preflipDat[REF!=ALT]
    # setrefDat=setrefDat[REF!=ALT]
    # n3=nrow(preflipDat)
    # n4=nrow(setrefDat)
    # inform(n1-n3, "records with identical ALT and REF removed from pre-flipped file.")
    # inform(n1-n3, "records with identical ALT and REF removed from post-flipped file.")
    # if(n1!=n2 | n3!=n4){
    #   # the two numbers in the infos above should be strictly identical
    #   # because ref flipping does not change alleles, just switch them

    ## the above lines are removed because this is exactly what no-comp-allele does in crossmap...
    
    if(n1!=n2){
      error("Unequal number of rows pre and post REF flip. Something went wrong.")
    }

    ## Merge tables
    preflipDat=preflipDat[,c("ID", "REF", "ALT")]
    setnames(preflipDat, c("REF", "ALT"), c("OLD_REF", "OLD_ALT"))

    preflipDat=merge(preflipDat, setrefDat, by='ID')
    if(n1!=nrow(preflipDat)){
      error("Some rows did not match between pre and post REF flip. Something went wrong.")
    }
    gc()

    ## Add flipped column
    preflipDat[,FLIPPED:=0]
    preflipDat[REF!=OLD_REF, FLIPPED:=1]

    #update INFO column
    preflipDat[, INFO := paste0(INFO,';FLIPPED=',FLIPPED, ';OLD=',OLD_REF, '_', OLD_ALT)]
    #delete FLIPPED and ...
    preflipDat[,c("FLIPPED", "OLD_REF", "OLD_ALT"):=NULL]
    setcolorder(preflipDat, c('#CHROM','POS', 'ID', 'REF','ALT','QUAL','FILTER','INFO') )
    setorder(preflipDat)
    # inform('NrowPreflipDataIn: ', format(nrow(preflipDat), big.mark=","))
    # inform('NrowSetrefDataIn: ', format(nrow(setrefDat), big.mark=","))
    # inform('NrowDataOut: ', format(nrow(preflipDat), big.mark=","))
    chrorder=unique(preflipDat[['#CHROM']])
    chrorder = str_sort(chrorder, numeric=TRUE)
    if(any(is.na(chrorder))){
      print(chrorder)
      error("Some chromosome names are not numeric")
    }
    preflipDat[,CHRNUM:=match(`#CHROM`, chrorder)]
    setorder(preflipDat, CHRNUM, POS)
    preflipDat[,CHRNUM:=NULL]
    return(preflipDat)
}
#Fxn Test
#test_checkflip=checkflip(preflipDat=dat_preflip,setrefDat=dat_setref)
#==================================================================================================================================
###################################################################################################################################

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
if(file.exists(args$preflipfile)){
    #outFilename1=paste0('tab_flipchecked_',tools::file_path_sans_ext(basename(args$preflipfile)), '.gz')
    outF1=basename(args$preflipfile)
    outFilename=str_replace(outF1, '.vcf.gz', ".flipchecked.vcf.gz")
    inform('Outfile name is: ', outFilename)
}else{
    error('File ', args$preflipfile, ' does not exist')
}

## test existence of input files
if(!file.exists(args$setreffile)){
  error("File ", args$setreffile, " does not exist")
}

## Read in data files and set column names
cat('\n')
#inform('reading in preflipData file .....')
dat_preflip=fread(args$preflipfile, skip="#CHROM\tPOS")

#inform('reading in setrefData file .....')
dat_setref=fread(args$setreffile, skip="#CHROM\tPOS")

###################################################################################################################################
## Fxn Call
outDat=checkflip(preflipDat=dat_preflip,setrefDat=dat_setref)
# Write final outputs
outPath=paste(outDir, outFilename, sep="/")
inform('writing output to: ', outPath)

# pipe through header (max 5000 lines of header)
conn=gzfile(args$setreffile)
header=grep("^#", readLines(conn, 5000), value=T)
if(! grepl( "^#CHROM\tPOS", header[length(header)])){
  error("Header of flipped file too long. Increase hardcoded header length value.")
}
close(conn)
conn=gzfile(outPath, "w")
writeLines(header, con=conn)
close(conn)

#write main data table
fwrite(outDat, file=outPath, sep="\t", quote=F, append=T)
