#!/bin/bash

############################################################################################################################
# Help
############################################################################################################################
Help()
{
   # Display Help
   echo "usage: ./liftover.sh [-h] [-f <inputfile>][-c <chainfile>] [-r <refgenome>] [-o <outputfile>] [-d <outputdir>] [-k]"
   echo
   echo "Converts genome coordinates between different reference assemblies for VCF format files
        (e.g.from human hg38 to hg19) using CrossMap utility."
   printf "Chromosomes Coordinates & Reference Alleles are updated to a new assembly,
         and all the other fields are not changed\n\n"
   echo
   echo "required arguments:"
   echo "-f,  --inputfile   file name (REQUIRED)"
   echo "-c,  --chainfile  chain file describes pairwise alignments between two genomes (REQUIRED):
                (https://genome.ucsc.edu/goldenPath/help/chain.html)"
   echo "-r,  --refgenome  chromosome sequences of target assembly in FASTA (REQUIRED):
                (https://en.wikipedia.org/wiki/FASTA_format)"
   echo "-d,  --outputdir output directory normalized files (default: 'liftnorm')"
   echo "optional arguments:"
   echo "-h,  --help   show  this help message and exit"
   echo "-o,  --outputfile  prefix output file name (default: 'lift_out')"
   echo "-k,  --no-comp-alleles  if set, CrossMap does NOT check if the reference allele is different from the
                     alternate allele. Disabled by default"
   printf "\nMore information about CrossMap can be found at:
         http://crossmap.sourceforge.net/#convert-vcf-format-files"
   echo
}
####################################################################################################################################################
# Main program
####################################################################################################################################################
####################################################################################################################################################
# Process input options
####################################################################################################################################################
# Set switch
compalleles=""
#####################################################################################################################################################
PARSED_ARGUMENTS=$(getopt -a -n liftover -o hkc:r:f:d:o: --long help,no-comp-alleles,chainfile:,refgenome:,inputfile:,outputdir:,outputfile: -- "$@")
VALID_ARGUMENTS=$?
if [ "$VALID_ARGUMENTS" != "0" ]; then
  Help
  exit 2
fi
#echo "PARSED_ARGUMENTS : $PARSED_ARGUMENTS"
eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    -k | --no-comp-alleles) compalleles=1 ; shift ;;
    -c | --chainfile)   chainfile="$2"   ; shift 2 ;;
    -r | --refgenome)   refgenome="$2"   ; shift 2 ;;
    -f | --inputfile)   inputfile="$2" ; shift 2 ;;
    -d | --outputdir)   outputdir="$2"   ; shift 2 ;;
    -o | --outputfile)  outputfile="$2"   ; shift 2 ;;
    -h | --help)
    Help
    exit   ;;
    # -- means the end of the arguments; drop this, and break out of the while loop
    --) shift; break ;;
    # If invalid options were passed, then getopt should have reported an error,
    # which we checked as VALID_ARGUMENTS when getopt was called...
    *) echo "Error: Invalid option: $1"
       Help
       exit ;;
  esac
done

#---------------------------
### Verify Input Options
#==========================
## Input file
printf "\n****** input file ********\n"
if [ -s $inputfile ]; then
    echo "input file is: $inputfile"
    # Get filename without suffix
    filename=$(basename ${inputfile} .gz)
else
    #echo "No input file provided"
    echo "Warning: {inputfile} ${inputfile} does not exist"
    exit 1
fi

## Chain file
printf "\n****** chain file ********\n"
if [ -s ${chainfile} ]; then
    echo "chain file is: $chainfile"
else
    echo "Warning: {chainfile} ${chainfile} does not exist"
    exit 1
fi

## refgenome
printf "\n****** reference genome ********\n"
if [ -s ${refgenome} ]; then
    echo "target reference genome is: ${refgenome}"
else
    echo "Warning: {refgenome} ${refgenome} does not exist"
    exit 1
fi

# output file name
printf "\n****** output filename ********\n"
if [[ -z ${outputfile} ]]; then
    # outputfile name prefix
    outputfile="lift_out"
fi
echo "Output file name tag is: $outputfile"

# output  directory created only if it doesn't exist
printf "\n****** output directory ********\n"
if [[ -z "${outputdir}" ]]; then
    outputdir="$(pwd)/liftnorm"
    echo "output directory not given, defaulting to $outputdir"
fi
mkdir -p ${outputdir}
abs_outdir=$(readlink -f $outputdir)
echo "output directory: $abs_outdir"

# keep/Filter variants
printf "\n*****  Keep/Discard variants  *****\n"
discardvariant_arg=$([ -z $discardvariants ] && echo "--no-comp-allele" || echo "")
if [ -z $discardvariant_arg ]; then
    echo "Discard positions where none of REF/ALT matches the reference"
else
    echo "Keep: No variants discarded (use -k to override) "
fi

#-------------------------------------------------------------------------------
##  Perform liftOver to target  genome assembly
#================================================================================
printf "\n*****      Performing LiftOver with CrossMap       *****\n"

compalleles_arg=$([ -z $compalleles ] &&  echo "" || echo "--no-comp-allele")
echo "Running CrossMap.py vcf $chainfile $inputfile $refgenome ${abs_outdir}/${outputfile}_${filename} $compalleles_arg"
python3 $(which CrossMap.py) vcf $chainfile $inputfile $refgenome ${abs_outdir}/${outputfile}_${filename} $compalleles_arg

# Delete liftover output header
if [[ -s ${abs_outdir}/${outputfile}_${filename} ]] && [[ $(grep '##liftOver' ${abs_outdir}/${outputfile}_${filename} | wc -l) -gt 0 ]]; then
    echo 'Cleaning liftover headers'
    sed -i '/^##liftOver/d;/^##originalFile/d;/^##targetRefGenome/d' ${abs_outdir}/${outputfile}_${filename}
else
  echo "Error: file ${abs_outdir}/${outputfile}_${filename} is empty or does not contain the expected headers."
  exit 1
fi

# Count input and output records/lines
if [[ -s ${abs_outdir}/${outputfile}_${filename} ]]; then
    echo "Counting input and output entries"
    echo "inputfile: ${inputfile}"
    echo "outputfile: ${abs_outdir}/${outputfile}_${filename}.gz"
    num_in=$(bcftools view --no-header ${inputfile} | wc -l  )
    num_out=$(bcftools view --no-header ${abs_outdir}/${outputfile}_${filename} | wc -l)
    num_unmap=$(awk '!/^#/' ${abs_outdir}/${outputfile}_${filename}.unmap | wc -l )
    echo
    echo "Lines numIn/numOut/numUnmap:   ${num_in}/${num_out}/${num_unmap}"
    # Compress output
    gzip -f ${abs_outdir}/${outputfile}_${filename}
else
    echo "liftover output ${abs_outdir}/${outputfile}_${filename} does not exist"
fi
echo "Finished Liftover `date`"
echo
######################################################################################################
## Example Call
# ./liftover.sh -c chainfile -r refgenome -f inputfile -d outdir > logfile 2>&1&
######################################################################################################
