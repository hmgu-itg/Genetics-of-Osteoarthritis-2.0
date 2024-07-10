#!/bin/bash

#####################################################################################################
# Help
######################################################################################################
Help()
{
   # Display Help
   echo "usage: ./normalize.sh [-h] [-f <inputfile>][-r <refgenome>][-o <outputfile>] [-d <outputdir>] "
   echo
   echo "Left-align and normalizes indels, checks if REF alleles match the reference using bcftools utility."
   echo "One or both of -l and -s must be enabled. If both are enabled, -s will be processed before -l."

   echo "required arguments:"
   echo "-f,  --inputfile  file name (REQUIRED)"
   echo "-r,  --refgenome  chromosome sequences of target assembly in FASTA (REQUIRED):
                (https://en.wikipedia.org/wiki/FASTA_format)"
   echo
   echo "optional arguments:"
   echo "-h,  --help   show  this help message and exit"
   echo "-o,  --outputfile  prefix for output file name (default: none, just add suffix)"
   echo "-d,  --outputdir output directory (default: 'liftnorm') "
   echo "-s,  --set-ref-allele set or fix reference allele from reference genome, enabled by default"
   echo "-l,  --left-align, left-align and normalize indels"
   echo "-i,  --keep-intermediate-files keep or delete intermediate files"
   echo
}

################################################################################################################################
# Main program
################################################################################################################################
#################################################################################################################################
# Process input options
##################################################################################################################################
# Set switch
#Fix flipped REF/ALT and
setrefallele=""
leftalign=""
####################################################################################################################################
# $# total number of passed arguments the $1 $2 ...
# -a -> alternative mode (e.g, -refgenome == -- refgenome)
PARSED_ARGUMENTS=$(getopt -a -n normalize -o hf:r:d:o:sl --long help,inputfile:,refgenome:,outputdir:,outputfile:,set-ref-allele,left-align -- "$@")
VALID_ARGUMENTS=$?
if [ "$VALID_ARGUMENTS" != "0" ]; then
  Help
  exit 2
fi
eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    -s | --set-ref-alllele)   setrefallele=1;      shift   ;;
    -l | --left-align)  leftalign=1;        shift ;;
    -f | --inputfile) inputfile="$2" ;      shift 2 ;;
    -r | --refgenome)   refgenome="$2"      ; shift 2 ;;
    -d | --outputdir)   outputdir="$2"      ; shift 2 ;;
    -o | --outputfile)   outputfile="$2"   ; shift 2 ;;
    -h | --help)
    Help
    exit   ;;
    --) shift; break ;;
    *) echo "Error: Invalid option: $1"
       Help
       exit ;;
  esac
done

#---------------------------
### Verify Input Options
#==========================
## Input file
#printf "\n****** input file ********\n"
if [ -s $inputfile ]; then
    echo "input file is: $inputfile"
    filename=$(basename ${inputfile} .gz)
else
    echo -e "ERROR: Provided inputfile \"${inputfile}\" does not exist."
    exit 1
fi

if [[ -z "$setrefallele" ]] && [[ -z "$leftalign" ]]; then
  echo "ERROR: neither -s or -l is specified. Nothing to do."
  exit 1
fi

## refgenome
#printf "\n****** reference genome ********\n"
if [ -s ${refgenome} ]; then
    echo "target reference genome is: ${refgenome}"
else
    echo -e "ERROR: Provided reference genome does not exist: \"${refgenome}\"."
    exit 1
fi

# output file name tag
#printf "\n****** output filename tag ********\n"
if [[ -z ${outputfile} ]]; then
    outputfile=""
    echo "Info: output file prefix not given, keeping filename as is."
else
  outputfile=$outputfile.
fi
#echo "Output file name tag: '$outputfile'"

# output  directory created only if it doesn't exist
#printf "\n****** output directory ********\n"
if [[ -z "${outputdir}" ]]; then
    outputdir="$(pwd)/liftnorm"
    echo "Info: output directory not given, defaulting to $outputdir"
fi
mkdir -p ${outputdir}
abs_outdir=$(readlink -f $outputdir)
echo "output directory: $abs_outdir"

#============================================================================================================
### 1. Fix/set REF allele from reference genome
#============================================================================================================
#setrefallele_arg=$([ -z $setrefallele ] && echo "" || echo -s)
if [[ -n "$setrefallele" ]]; then
    #printf "\n*****      Fix/Set REF from reference genome with BCFtools        *****\n"
    # -c > check-ref; what to do when incorrect or missing REF allele is encountered:
    # exit(e), warn(w), exclude(x), or set/fix (s) bad sites
    #outfile_fixref=$abs_outdir/${filename/.vcf.gz/.setref.vcf.gz}
    #bcftools norm --do-not-normalize -cs --fasta-ref $refgenome -Oz -o $outfile_fixref $inputfile
    outfile_setref=$abs_outdir/${filename/.vcf/.setref.vcf.gz}
    #echo "Step 1: Flipping reference."
    echo "Running: bcftools norm -Oz -o $outfile_setref --do-not-normalize -cs --fasta-ref $refgenome $inputfile"
    echo
    bcftools norm -oZ -o $outfile_setref --do-not-normalize -cs --fasta-ref $refgenome $inputfile

    if [[ ! -s $outfile_setref ]]; then
        echo "Error: File $outfile_setref failed to generate. Unable to fix REF."
        exit 1
    fi
    # if both options are enabled:
    inputfile=$outfile_setref
fi


#===============================================================================================
## 2. Re-align Indels
#===============================================================================================
#leftalign_arg=$( [ -z $leftalign ] && echo "" || echo -l )
if [[ -n "$leftalign" ]]; then
  printf "\n*****      Performing Normalization with BCFtools        *****\n"
  echo 'Left-align and Normalizing Indels'
  outfile_norm=$abs_outdir/$outputfile${filename/.vcf/.leftaligned.vcf.gz}
  echo
  echo "Running: bcftools norm -Oz -o $outfile_norm -cw -f $refgenome $inputfile"
  echo
  bcftools norm -Oz -o $outfile_norm -cw -f $refgenome $inputfile

  if [[ ! -s $outfile_norm ]]; then
      echo "Error: File $outfile_norm failed to generate. Unable to realign indels and normalize."
      echo "`date`"
      exit 1
  else
      echo "Finished Normalization `date`"
  fi
fi

echo "Finished Fixing REF `date`"
#############################################################################################
#Example Call
# ./normalize.sh -f inputfile.vcf.gz -r refgenome.fa -d outdir  > logfile 2>&1&
## output
# inputfilename.formatted.setref.vcf.gz
# tab.setref.inputfilename.formatted.setref.gz
# inputfilename.setref.norm.vcf.gz
# tab.norm.inputfilename.norm.gz
##############################################################################################
