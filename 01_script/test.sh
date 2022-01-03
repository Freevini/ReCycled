#!/bin/bash


############################################################
############################################################
# set defaults and parameters                                        #
############################################################
############################################################

# Set defaults
circleriesPATH=$(dirname $(realpath -s $0))
outputFolderName=./
threads=4
verbose="N"
Version="N"
tmpss="N"
# Set Infos
AUTHOR="Vincent Somerville"
EXE="Circleries"
VERSION="V0.0.1"
starttime=$(date)
start=`date +%s`
outName="N"
force="N"
limitOverlaps=5 #not incorporated into parameters

shortreads_1=""
shortreads_2=""

unset genomeFASTAname
unset longreads
#unset shortreads_1
#unset shortreads_2

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Circleries: checks the completness, circularity and bacterial origins of contigs and start aligns them at dnaA if possible"
   #echo -e "This is "$EXE $VERSION
   echo
   echo "minimal syntax: Circleries -i <genome_input.fasta> -l <raw_long_read.fastq.gz>"
   echo "options:"
   echo
   echo "INPUT"
   echo "   -i     input genome name (in fasta format)" #genomeFASTAname
   echo "   -l     long read file (fq or fq.gz)" #longreads
   echo "   -f     short read forward read (read 1) (fq or fq.gz)" #shortreads_1
   echo "   -r     short read reverse read (read 2) (fq or fq.gz)" #shortreads_2
   echo
   echo "OUTPUT"
   echo "   -d     output directory [.]" #outputFolderName
   echo "   -o     output file name  [results]" #threads #outName
   echo
   echo "RUNNING OPTIONS"
   echo "   -p     circleries script path (with all dependencies) [PATH]" #circleriesPATH
   echo "   -t     number of threads to use [4]" #threads
   echo "   -x     keep all tmp files created [N]" #threads
   echo "   -F     Force everything to run again [N]" #threads
   echo
   echo "INFOS"
   echo "   -h     help option" #shortreads_2
   echo "   -v     verbose [N]" #verbose
   echo "   -V     print Version [N]" #Version

   echo
}

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# welcome message

#echo -e "This is "$EXE $VERSION
#echo -e "Written by "$AUTHOR
#echo -e "Local time is" $starttime
#echo -e "Operating system is "$OSTYPE


############################################################
# Process the input options. Add options as needed.        #
############################################################
# Get the options
while getopts ":h :p: :i: :d: :o: :t: :l: :r: :f: :v :V :x :F" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      p) #circleries scirpt path
         circleriesPATH=$OPTARG;;
      i) #input genome name (in fasta format)
         genomeFASTAname=$OPTARG;;
      d) #output direcotry [.]
         outputFolderName=$OPTARG;;
      o) #output file name  [results]
         outName=$OPTARG;;
      t) #number of threads to use [4]" #threads
         threads=$OPTARG;;
      l) #long read file
         longreads=$OPTARG;;
      r) #short read reverse read (read 2)
         shortreads_2=$OPTARG;;
      f) #short read forward read (read 1)
         shortreads_1=$OPTARG;;
      v) #print verbose
         verbose="Y";;
      V) #print version
         Version="Y";;
      x) #print version
         tmpss="Y";;
      F) #print version
         force="Y";;
     \?) # Invalid option
         echo "Error: Invalid option"
         Help
         exit;;
   esac
done
shift $((OPTIND -1)) #removing previously inputed options
#echo $tmpss

############################################################
# If version wanted print version    #
############################################################
#echo $Version
if [ "$Version" == "Y" ]
then
        echo -e "This is" $EXE $VERSION
        exit 1
fi


# ###########################################################
#check for mandatory parameters and inputs and names
############################################################

: ${longreads:?Missing: "-l" which is the long read file. This information is mandatory}
: ${genomeFASTAname:?Missing: "-i" which is the input genome name (\in fasta format). This information is mandatory}

##--------------------------------------------prepare outName if necessary----------------------------------------
if [ "$outName" == "N" ]
then
    outName=$(echo $genomeFASTAname|sed 's/.*\///g'| awk '{print "Start_aligned_"$0}')
fi


# ###########################################################
#Test
############################################################

#echo "hello World"

# ###########################################################
#prep and dependencies
############################################################



# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# welcome message
echo
echo -e "This is "$EXE $VERSION
echo -e "Written by "$AUTHOR
echo -e "Local time is" $starttime
echo -e "Operating system is "$OSTYPE
echo
echo -e "Options:"
echo -e "   Input Genome:         "${genomeFASTAname}
echo -e "   long read file:       "${longreads}
echo -e "   short read forward:   "${shortreads_1}
echo -e "   short read reverse:   "${shortreads_2}
echo -e "   Output directory:     "${outputFolderName}
echo -e "   Output file name:     "${outName}
echo -e "   Circleries path :     "${circleriesPATH}
echo -e "   Threads :             "${threads}
echo -e "   Force rerun :         "${force}


# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# Dependencies
echo
echo -e "Dependencies:"


fasta_shift_path=$(echo -e ${circleriesPATH}"/02_dependencies/fasta_shift")
seqkit_path=$(echo -e ${circleriesPATH}"/02_dependencies/seqkit")
startalining_genes=$(echo -e ${circleriesPATH}"/05_startAlign_data/starting_genes.fasta ")

#${circleriesPATH}/05_startAlign_data/starting_genes.fasta

#$fasta_shift_path --help
#type baloooihok >/dev/null 2>&1 && echo -e "   baloooihok......OK" || { echo >&2 "   baloooihok is not installed.  Aborting."; exit 1; }

type minimap2 >/dev/null 2>&1 && echo -e "   Minimap2......OK" || { echo >&2 "   Minimap2 is not installed.  Aborting."; exit 1; }
#type $fasta_shift_path >/dev/null 2>&1 && echo -e "   fasta_shift......OK" || { echo >&2 "   fasta_shift is not installed.  Aborting."; exit 1; }
#type revseq >/dev/null 2>&1 && echo -e "   revseq......OK" || { echo >&2 "   revseq is not installed.  Aborting."; exit 1; }
#type bedtools >/dev/null 2>&1 && echo -e "   bedtools......OK" || { echo >&2 "   bedtools is not installed.  Aborting."; exit 1; }
type $seqkit_path >/dev/null 2>&1 && echo -e "   seqkit......OK" || { echo >&2 "   seqkit is not installed.  Aborting."; exit 1; }
[ -f $startaligning_genes ] && echo "   startaligning_genes......OK" || { echo "   startaligning_genes is not installed.  Aborting."; exit 1; }


#type $startalining_genes >/dev/null 2>&1 && echo -e "   startalining_genes......OK" || { echo >&2 "   startalining_genes is not installed.  Aborting."; exit 1; }


echo

###########################################################
#genome preperation
############################################################

echo -e "Preparing genome"


##--------------------------------------------set the direcotry----------------------------------------
if [[ "$force" == "Y" ]]; then
  rm -r ${outputFolderName}/tmp_${outName}
fi

mkdir -p ${outputFolderName}/tmp_${outName}/genome/

##--------------------------------------------put every contig on one line----------------------------------------
sed '/>/s/\t/_/g' ${genomeFASTAname}|sed '/>/s/ /_/g'|sed '/>/s/$/\t/g' | tr -d '\n' | sed 's/\t/\n/g' | sed 's/>/\n>/g' | sed '/^[[:space:]]*$/d'  >  ${outputFolderName}/tmp_${outName}/genome/tmp_wide_all.fasta

##--------------------------------------------split every contig into one file------------------------------------
echo -e "       Seperating contigs..."

for header in $(grep ">" ${outputFolderName}/tmp_${outName}/genome/tmp_wide_all.fasta)
do
header_short=$(echo ${header}|sed 's/>//g' )
echo -e "           "$header_short
grep "${header}" -A 1 ${outputFolderName}/tmp_${outName}/genome/tmp_wide_all.fasta > ${outputFolderName}/tmp_${outName}/genome/${header_short}.fasta
done
echo

###########################################################
#dnaA search
############################################################

echo -e "Searching for OriC"

##--------------------------------------------minimap2 dnaA to sequences------------------------------------

#rm -r  ${outputFolderName}/tmp_${genomeFASTAname}//start_alignment_mapping/
mkdir -p  ${outputFolderName}/tmp_${outName}//start_alignment_mapping/
echo -e "       On contig..."

echo -e "#contigName\tContigOrigin\tNumberOfdnaAmappings\tOrientationOfMapping\tStartOfTargetMapping\tContigLength" >  ${outputFolderName}/tmp_${outName}//start_alignment_mapping//StartAlignment_contigs.minimap


for header in $(grep ">" ${outputFolderName}/tmp_${outName}/genome/tmp_wide_all.fasta |sed 's/>//g')
do
#  echo -e "           "$header
##map dnaA genes to the individual references and remove non-mapping or low quality mapping reads (mapq>50)
minimap2 ${outputFolderName}/tmp_${outName}/genome/${header}.fasta $startalining_genes  > ${outputFolderName}/tmp_${outName}//start_alignment_mapping//${header}_unfiltered.minimap 2> ${outputFolderName}/tmp_${outName}//start_alignment_mapping//${header}.minimap.log
awk -F "\t" '{OFS="\t"}{if($12>50) print $0}' ${outputFolderName}/tmp_${outName}//start_alignment_mapping//${header}_unfiltered.minimap > ${outputFolderName}/tmp_${outName}//start_alignment_mapping//${header}.minimap


##--------------------------------------------minimap2 analysis------------------------------------


#find out if there are mappings and in which orientation
number_of_mappings=$(wc -l ${outputFolderName}/tmp_${outName}//start_alignment_mapping//${header}.minimap|cut -d ' ' -f 1)
orientation=$(sort -k10 -n -r ${outputFolderName}/tmp_${outName}//start_alignment_mapping//${header}.minimap|head -1 |cut -f 5)


# create a info table containing information about the presence of dnaA and if present the orientation

if [ $number_of_mappings == "0" ]
then
  echo -e "           "$header"......no dnaA map found"
      contig_length=$(grep -v ">" ${outputFolderName}/tmp_${outName}/genome/${header}.fasta|wc -c )

       echo -e ${header}"\tnon-dnaA-containing-contig\t"${number_of_mappings}"\tNA\tNA\t"${contig_length} >>  ${outputFolderName}/tmp_${outName}//start_alignment_mapping//StartAlignment_contigs.minimap

elif [ $orientation == "+" ]
then
  echo -e "           "$header"......dnaA map found"
  #the most left alignment this shows the start of the gene
  #sort -k8 -n ${outputFolderName}/start_alignment_mapping/${header}.minimap

  start=$(sort -k8 -n ${outputFolderName}/tmp_${outName}//start_alignment_mapping//${header}.minimap|head -1 |awk -F "\t" '{OFS="\t"}{print $8-$3-5}')
        contig_length=$(sort -k8 -n ${outputFolderName}/tmp_${outName}//start_alignment_mapping//${header}.minimap|head -1 |cut -f 7)

        echo -e ${header}"\tdnaA-containing-contig\t"${number_of_mappings}"\t"${orientation}"\t"${start}"\t"${contig_length} >>  ${outputFolderName}/tmp_${outName}//start_alignment_mapping//StartAlignment_contigs.minimap

else
      echo -e "           "$header"......dnaA map found"

    #the most right alignment this shows the start of the gene
    #sort -k9 -n ${outputFolderName}/start_alignment_mapping/${header}.minimap

    start=$(sort -k9 -n -r ${outputFolderName}/tmp_${outName}//start_alignment_mapping//${header}.minimap|head -1 |awk -F "\t" '{OFS="\t"}{print $9+$2-$4+5}')
        contig_length=$(sort -k9 -n -r  ${outputFolderName}/tmp_${outName}//start_alignment_mapping//${header}.minimap|head -1 |cut -f 7)

        echo -e ${header}"\tdnaA-containing-contig\t"${number_of_mappings}"\t"${orientation}"\t"${start}"\t"${contig_length} >>  ${outputFolderName}/tmp_${outName}//start_alignment_mapping//StartAlignment_contigs.minimap

fi

done

echo


###########################################################
#1. Check circularity by Overlaping contig edges
############################################################

echo -e "1. Checking circularity"
echo -e "       Overlaping contig edges..."


#rm -r ${outputFolderName}/tmp_${genomeFASTAname}//Start_end_overlap/
mkdir -p ${outputFolderName}/tmp_${outName}/Start_end_overlap/{fastas,mapping}/

echo -e "OverlapingContigEdges" > ${outputFolderName}/tmp_${outName}/Start_end_overlap/tmp.analysis

for header in $(grep ">" ${outputFolderName}/tmp_${outName}/genome/tmp_wide_all.fasta |sed 's/>//g')
do
#echo -e ${header}
#contig_length=$(sort -k9 -n -r  ${outputFolderName}/tmp_${outName}//start_alignment_mapping//${header}.minimap|head -1 |cut -f 7)
contig_length=$(grep -v ">" ${outputFolderName}/tmp_${outName}/genome/${header}.fasta|wc -c )

##--------------------------------------------extract start of contig----------------------------------

echo -e ${header}"\t1\t1000" > ${outputFolderName}/tmp_${outName}/Start_end_overlap/fastas/${header}_start.bed
$seqkit_path subseq ${outputFolderName}/tmp_${outName}/genome/${header}.fasta --bed ${outputFolderName}/tmp_${outName}/Start_end_overlap/fastas/${header}_start.bed > ${outputFolderName}/tmp_${outName}/Start_end_overlap/fastas/${header}_start.fasta 2> ${outputFolderName}/tmp_${outName}/Start_end_overlap/fastas/${header}_start.log
#end extract

##--------------------------------------------extract end of contig-----------------------------------

echo -e ${header}| awk -F "\t" -v contigLength="$contig_length" '{OFS="\t"}{print $1,contigLength-1001,contigLength-1}' > ${outputFolderName}/tmp_${outName}/Start_end_overlap/fastas/${header}_end.bed
$seqkit_path subseq ${outputFolderName}/tmp_${outName}/genome/${header}.fasta --bed ${outputFolderName}/tmp_${outName}/Start_end_overlap/fastas/${header}_end.bed > ${outputFolderName}/tmp_${outName}/Start_end_overlap/fastas/${header}_end.fasta 2> ${outputFolderName}/tmp_${outName}/Start_end_overlap/fastas/${header}_end.log

##--------------------------------------------mapping start and end----------------------------------

minimap2 ${outputFolderName}/tmp_${outName}/Start_end_overlap/fastas/${header}_start.fasta ${outputFolderName}/tmp_${outName}/Start_end_overlap/fastas/${header}_end.fasta > ${outputFolderName}/tmp_${outName}/Start_end_overlap/mapping/${header}_mapping_unfiltered.paf 2> ${outputFolderName}/tmp_${outName}/Start_end_overlap/mapping/${header}_mapping.log
awk -F "\t" '{OFS="\t"}{if($12>50) print $0}'  ${outputFolderName}/tmp_${outName}/Start_end_overlap/mapping/${header}_mapping_unfiltered.paf > ${outputFolderName}/tmp_${outName}/Start_end_overlap/mapping/${header}_mapping.paf

##--------------------------------------------analysis----------------------------------
readCounts=$(wc -l ${outputFolderName}/tmp_${outName}/Start_end_overlap/mapping/${header}_mapping.paf)
echo -e ${readCounts} | awk  '{if ($1==0)	print "N" ; else 	print "Y"}' >> ${outputFolderName}/tmp_${outName}/Start_end_overlap/tmp.analysis
echo -e ${readCounts} | awk -v contigName="$header" '{if ($1==0)	print "           "contigName"......Not Overlapping" ; else 	print  "           "contigName"......Overlapping"}'


done

##--------------------------------------------Merge with previous analysis----------------------------------

paste -d "\t" ${outputFolderName}/tmp_${outName}//start_alignment_mapping//StartAlignment_contigs.minimap ${outputFolderName}/tmp_${outName}/Start_end_overlap/tmp.analysis >  ${outputFolderName}/tmp_${outName}/Start_end_overlap/analysis_circularity_extended
#cat  ${outputFolderName}/tmp_${outName}/Start_end_overlap/analysis_circularity_extended
echo

###########################################################
#2. Check circularity by aligning long reads
############################################################

echo -e "2. Checking circularity"
echo -e "       Aligning long reads...This may take some time!"



if [ -f "${outputFolderName}/tmp_${outName}/Start_end_readmapping/long_read/tmp.analysis" ]; then
#if [[ "$longreads" == "NA" ]]; then #if longreads are alternative
    echo  "           skip because output file already exists"

elif [[ "$longreads" == "NA" ]]; then
      echo  "           skip because no long reads supplied"


      echo -e "OverlapingLongReads" > ${outputFolderName}/tmp_${outName}/tmp_longReads.analysis
      grep "^#" -v ${outputFolderName}/tmp_${outName}/Start_end_overlap/analysis_circularity_extended | awk '{print "NA"}' >>  ${outputFolderName}/tmp_${outName}/tmp_longReads.analysis
      ##--------------------------------------------Merge with previous analysis----------------------------------
      paste -d "\t" ${outputFolderName}/tmp_${outName}/Start_end_overlap/analysis_circularity_extended ${outputFolderName}/tmp_${outName}/tmp_longReads.analysis >  ${outputFolderName}/tmp_${outName}/Start_end_readmapping/long_read/analysis_circularity_extended


    else

#rm -r ${outputFolderName}/tmp_${outName}/Start_end_readmapping/long_read/
mkdir -p ${outputFolderName}/tmp_${outName}/Start_end_readmapping/long_read/{fastas,mapping}/
echo -e "OverlapingLongReads" > ${outputFolderName}/tmp_${outName}/Start_end_readmapping/long_read/tmp.analysis

for header in $(grep ">" ${outputFolderName}/tmp_${outName}/genome/tmp_wide_all.fasta |sed 's/>//g')
do
#echo -e ${header}

##--------------------------------------------merge start and end divided by NNNNs----------------------------------

echo -e ">"${header}"_EndAndStart" > ${outputFolderName}/tmp_${outName}/Start_end_readmapping/long_read/fastas/${header}_EndAndStart.fasta
grep ">" -v ${outputFolderName}/tmp_${outName}/Start_end_overlap/fastas//${header}_end.fasta >> ${outputFolderName}/tmp_${outName}/Start_end_readmapping/long_read/fastas/${header}_EndAndStart.fasta
echo -e "NNNNNNNNNN" >> ${outputFolderName}/tmp_${outName}/Start_end_readmapping/long_read/fastas/${header}_EndAndStart.fasta
grep ">" -v ${outputFolderName}/tmp_${outName}/Start_end_overlap/fastas//${header}_start.fasta >> ${outputFolderName}/tmp_${outName}/Start_end_readmapping/long_read/fastas/${header}_EndAndStart.fasta


##make the assembly fit together nicely (60 nucleotides wide)
awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' ${outputFolderName}/tmp_${outName}/Start_end_readmapping/long_read/fastas/${header}_EndAndStart.fasta |sed 's/\t/\n/g' |awk -v FS= '/^>/{print;next}{for (i=0;i<=NF/60;i++) {for (j=1;j<=60;j++) printf "%s", $(i*60 +j); print ""}}'  > ${outputFolderName}/tmp_${outName}/Start_end_readmapping/long_read/fastas/${header}_EndAndStart_fitted.fasta

rm ${outputFolderName}/tmp_${outName}/Start_end_readmapping/long_read/fastas/${header}_EndAndStart.fasta


##-------------------------------------------- map to the new constructs----------------------------------

##map dnaA genes to the individual references and remove non-mapping or low quality mapping reads (mapq>50)
minimap2 -x map-ont -t ${threads} ${outputFolderName}/tmp_${outName}/Start_end_readmapping/long_read/fastas/${header}_EndAndStart_fitted.fasta ${longreads} > ${outputFolderName}/tmp_${outName}/Start_end_readmapping/long_read/mapping/${header}_mapping_unfiltered.paf 2> ${outputFolderName}/tmp_${outName}/Start_end_readmapping/long_read/mapping/${header}_mapping.log
awk -F "\t" '{OFS="\t"}{if($12>50) print $0}'  ${outputFolderName}/tmp_${outName}/Start_end_readmapping/long_read/mapping/${header}_mapping_unfiltered.paf >  ${outputFolderName}/tmp_${outName}/Start_end_readmapping/long_read/mapping/${header}_mapping.paf


##-------------------------------------------- set limitation that reads must map more then 1400bp----------------------------------


awk -F "\t" '{OFS="\t"}{if($10>1400) print $0}' ${outputFolderName}/tmp_${outName}/Start_end_readmapping/long_read/mapping/${header}_mapping.paf > ${outputFolderName}/tmp_${outName}/Start_end_readmapping/long_read/mapping/${header}_mapping_overlaping_reads.paf
overlaps=$(wc -l ${outputFolderName}/tmp_${outName}/Start_end_readmapping/long_read/mapping/${header}_mapping_overlaping_reads.paf |cut -d ' ' -f 1)


echo -e ${overlaps} | awk -v limitOverlapzz="$limitOverlaps" '{if ($1>limitOverlapzz)	print "Y" ; else 	print "N"}' >> ${outputFolderName}/tmp_${outName}/Start_end_readmapping/long_read/tmp.analysis
echo -e ${overlaps} | awk -v limitOverlapzz="$limitOverlaps"  -v contigName="$header" '{if ($1>limitOverlapzz)	print "           "contigName"......found "$1" mapping long reads" ; else 	print  "           "contigName"......No mapping long reads found"}'



done



##--------------------------------------------Merge with previous analysis----------------------------------
paste -d "\t" ${outputFolderName}/tmp_${outName}/Start_end_overlap/analysis_circularity_extended ${outputFolderName}/tmp_${outName}/Start_end_readmapping/long_read/tmp.analysis >  ${outputFolderName}/tmp_${outName}/Start_end_readmapping/long_read/analysis_circularity_extended


fi #finish of else if
#cat  ${outputFolderName}/tmp_${outName}/Start_end_readmapping/long_read/analysis_circularity_extended



echo

###########################################################
#3. Check circularity by aligning short reads
############################################################

echo -e "3. Checking circularity"
echo -e "       Aligning short reads...This may take some time!"


if [[ "$shortreads_1" == "" ]];
then
    echo  "           skip because no short reads supplied"

    echo -e "OverlapingLongReads" > ${outputFolderName}/tmp_${outName}/tmp_shortReads.analysis
    grep "^#" -v ${outputFolderName}/tmp_${outName}/Start_end_overlap/analysis_circularity_extended | awk '{print "NA"}' >>  ${outputFolderName}/tmp_${outName}/tmp_shortReads.analysis
    ##--------------------------------------------Merge with previous analysis----------------------------------
    paste -d "\t" ${outputFolderName}/tmp_${outName}/Start_end_readmapping/long_read/analysis_circularity_extended ${outputFolderName}/tmp_${outName}/tmp_shortReads.analysis >  ${outputFolderName}/tmp_${outName}/Start_end_readmapping/short_read/analysis_circularity_extended


  elif [ -f "${outputFolderName}/tmp_${outName}/Start_end_readmapping/short_read/tmp.analysis" ]; then
    echo  "           skip because output file already exists"
    else


    #rm -r ${outputFolderName}/tmp_${outName}/Start_end_readmapping/long_read/
    mkdir -p ${outputFolderName}/tmp_${outName}/Start_end_readmapping/short_read/{fastas,mapping}/
    echo -e "OverlapingShortReads" > ${outputFolderName}/tmp_${outName}/Start_end_readmapping/short_read/tmp.analysis

    for header in $(grep ">" ${outputFolderName}/tmp_${outName}/genome/tmp_wide_all.fasta |sed 's/>//g')
    do
    #echo -e ${header}
    #echo  "           "${header}

    ##--------------------------------------------merge start and end divided by NNNNs----------------------------------

    echo -e ">"${header}"_EndAndStart" > ${outputFolderName}/tmp_${outName}/Start_end_readmapping/short_read/fastas/${header}_EndAndStart.fasta
    grep ">" -v ${outputFolderName}/tmp_${outName}/Start_end_overlap/fastas//${header}_end.fasta >> ${outputFolderName}/tmp_${outName}/Start_end_readmapping/short_read/fastas/${header}_EndAndStart.fasta
    echo -e "NNNNNNNNNN" >> ${outputFolderName}/tmp_${outName}/Start_end_readmapping/short_read/fastas/${header}_EndAndStart.fasta
    grep ">" -v ${outputFolderName}/tmp_${outName}/Start_end_overlap/fastas//${header}_start.fasta >> ${outputFolderName}/tmp_${outName}/Start_end_readmapping/short_read/fastas/${header}_EndAndStart.fasta

    ##make the assembly fit together nicely (60 nucleotides wide)
    awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' ${outputFolderName}/tmp_${outName}/Start_end_readmapping/short_read/fastas/${header}_EndAndStart.fasta |sed 's/\t/\n/g' |awk -v FS= '/^>/{print;next}{for (i=0;i<=NF/60;i++) {for (j=1;j<=60;j++) printf "%s", $(i*60 +j); print ""}}'  > ${outputFolderName}/tmp_${outName}/Start_end_readmapping/short_read/fastas/${header}_EndAndStart_fitted.fasta

    rm ${outputFolderName}/tmp_${outName}/Start_end_readmapping/short_read/fastas/${header}_EndAndStart.fasta



    ##--------------------------------------------map to the new constructs---------------------------------

    #map dnaA genes to the individual references and remove non-mapping or low quality mapping reads (mapq>50)
   minimap2 -x sr -t ${threads} ${outputFolderName}/tmp_${outName}/Start_end_readmapping/short_read/fastas/${header}_EndAndStart_fitted.fasta ${shortreads_1} ${shortreads_2} > ${outputFolderName}/tmp_${outName}/Start_end_readmapping/short_read/mapping/${header}_mapping_unfiltered.paf 2> ${outputFolderName}/tmp_${outName}/Start_end_readmapping/short_read/mapping/${header}_mapping.paf.log
   awk -F "\t" '{OFS="\t"}{if($12>50) print $0}' ${outputFolderName}/tmp_${outName}/Start_end_readmapping/short_read/mapping/${header}_mapping_unfiltered.paf > ${outputFolderName}/tmp_${outName}/Start_end_readmapping/short_read/mapping/${header}_mapping.paf
   #wc -l ${outputFolderName}/tmp_${outName}/Start_end_readmapping/short_read/mapping/${header}_mapping.paf


    ##-------------------------------------------- reads that map longer than 50% of the read and map 10bp on either side----------------------------------

    reads_mapping=$(awk '{if($8<990 && $9>1020 && $10 > 0.5*$2) print $0}'  ${outputFolderName}/tmp_${outName}/Start_end_readmapping/short_read/mapping/${header}_mapping.paf |wc -l |cut -d ' ' -f 1)
    Read_counts=$(wc -l  ${outputFolderName}/tmp_${outName}/Start_end_readmapping/short_read/mapping/${header}_mapping.paf |cut -d ' ' -f 1)
    #echo -e ${header}" has in total "${reads_mapping}" supporting the ciruclarity with "${Read_counts}" mappig to the edges"

    echo -e ${reads_mapping} | awk -v limitOverlapzz="$limitOverlaps" '{if ($1>limitOverlapzz)	print "Y" ; else 	print "N"}' >> ${outputFolderName}/tmp_${outName}/Start_end_readmapping/short_read/tmp.analysis
    echo -e ${reads_mapping} | awk -v limitOverlapzz="$limitOverlaps"  -v contigName="$header" '{if ($1>limitOverlapzz)	print "           "contigName"......found "$1" mapping short reads" ; else 	print  "           "contigName"......No mapping short reads found"}'

    done


        ##--------------------------------------------Merge with previous analysis----------------------------------
        paste -d "\t" ${outputFolderName}/tmp_${outName}/Start_end_readmapping/long_read/analysis_circularity_extended ${outputFolderName}/tmp_${outName}/Start_end_readmapping/short_read/tmp.analysis >  ${outputFolderName}/tmp_${outName}/Start_end_readmapping/short_read/analysis_circularity_extended

    fi

  #  cat  ${outputFolderName}/tmp_${outName}/Start_end_readmapping/short_read/analysis_circularity_extended
echo


###########################################################
#Final verdict
############################################################

echo -e "Circularizability final verdict"

echo -e "Circularizable" > ${outputFolderName}/tmp_${outName}/Start_end_readmapping/tmp.circularizable

for header in $(grep ">" ${outputFolderName}/tmp_${outName}/genome/tmp_wide_all.fasta |sed 's/>//g')
do
#echo $header
awk -F "\t" -v contigNames="$header" '{OFS="\t"}{if($1==contigNames )print $0}' ${outputFolderName}/tmp_${outName}/Start_end_readmapping/short_read/analysis_circularity_extended  |awk -F "\t" -v contigNames="$contigName" '{OFS="\t"}{if($2=="dnaA-containing-contig"  && ($7=="Y"||$8=="Y"||$9=="Y"))print "Y"; else print "N"}' >> ${outputFolderName}/tmp_${outName}/Start_end_readmapping/tmp.circularizable
done

##--------------------------------------------Merge with previous analysis----------------------------------
paste -d "\t"  ${outputFolderName}/tmp_${outName}/Start_end_readmapping/short_read/analysis_circularity_extended ${outputFolderName}/tmp_${outName}/Start_end_readmapping/tmp.circularizable  >  ${outputFolderName}/analysis_circularity_extended.log
cat  ${outputFolderName}/analysis_circularity_extended.log

echo

###########################################################
#Circleries and start align
############################################################
echo -e "Circularising"

 mkdir -p ${outputFolderName}/tmp_${outName}/StartAlignedContigs/
 mkdir -p ${outputFolderName}/tmp_${outName}/tmp/

for contigName in $(grep "^#" -v  ${outputFolderName}/analysis_circularity_extended.log |awk -F "\t" '{OFS="\t"}{if($10=="Y")print $1}')
do

#echo -e "contig name: \t\t" ${contigName}
 orientation=$(awk -F "\t" -v contigzz="$contigName" '{OFS="\t"}{if($1==contigzz)print $4}' ${outputFolderName}/analysis_circularity_extended.log )
 contigStart=$(awk -F "\t" -v contigzz="$contigName" '{OFS="\t"}{if($1==contigzz)print $5}'   ${outputFolderName}/analysis_circularity_extended.log )
 contigLength=$(awk -F "\t" -v contigzz="$contigName" '{OFS="\t"}{if($1==contigzz)print $6}'  ${outputFolderName}/analysis_circularity_extended.log )

  if [ $orientation == "+" ]
    then


    #  $fasta_shift_path -i ${outputFolderName}/tmp_${outName}/genome/${contigName}.fasta -p ${contigStart} > \
    #          ${outputFolderName}/tmp_${outName}/StartAlignedContigs/${contigName}_startAligned_old.fasta

        cat ${outputFolderName}/tmp_${outName}/genome/${contigName}.fasta  | $seqkit_path restart -i ${contigStart} > \
        ${outputFolderName}/tmp_${outName}/StartAlignedContigs/${contigName}_startAligned.fasta

              echo -e "           forward oriented and start align at......"${contigStart}

  else

      #  echo "reverse oriented and start align at:  "${contigStart}

      #$fasta_shift_path -i ${outputFolderName}/tmp_${outName}/genome/${contigName}.fasta -p ${contigStart} > \
      #      ${outputFolderName}/tmp_${outName}/tmp/${contigName}_startAligned_wrongOrientation_old.fasta

      cat ${outputFolderName}/tmp_${outName}/genome/${contigName}.fasta  | $seqkit_path restart -i ${contigStart} > \
                  ${outputFolderName}/tmp_${outName}/tmp/${contigName}_startAligned_wrongOrientation.fasta

    #  revseq -sequence ${outputFolderName}/tmp_${outName}/tmp/${contigName}_startAligned_wrongOrientation.fasta -outseq \
    #      ${outputFolderName}/tmp_${outName}/StartAlignedContigs/${contigName}_startAligned.fasta -notag 2> ${outputFolderName}/tmp_${outName}/tmp/${contigName}_startAligned_wrongOrientation.log

      $seqkit_path seq -t dna ${outputFolderName}/tmp_${outName}/tmp/${contigName}_startAligned_wrongOrientation.fasta -r -p > \
      ${outputFolderName}/tmp_${outName}/StartAlignedContigs/${contigName}_startAligned.fasta 2> ${outputFolderName}/tmp_${outName}/tmp/${contigName}_startAligned_wrongOrientation.log


          echo -e "           reverse oriented and start align at......"${contigStart}

  fi
  #add startaligned tag
  sed -i '/^>/ s/$/_StartAligned/' ${outputFolderName}/tmp_${outName}/StartAlignedContigs/${contigName}_startAligned.fasta


done

echo

###########################################################
#Quality control
############################################################
echo -e "Quality Control"

###-------------------QC if DNA is really at first position

echo -e "#contigName\tContigOrigin\tNumberOfdnaAmappings\tOrientationOfMapping\tStartOfTargetMapping\tContigLength" >   ${outputFolderName}/tmp_${outName}/start_alignment_mapping/After_StartAlignment_contigs.minimap

for contigName in $(grep "^#" -v  ${outputFolderName}/analysis_circularity_extended.log |awk -F "\t" '{OFS="\t"}{if($10=="Y")print $1}')
do


  ##map dnaA genes to the individual references and remove non-mapping or low quality mapping reads (mapq>50)
  minimap2 ${outputFolderName}/tmp_${outName}/StartAlignedContigs/${contigName}_startAligned.fasta $startalining_genes  > ${outputFolderName}/tmp_${outName}/start_alignment_mapping/${contigName}_afterStartAlignment_unfiltered.minimap 2> ${outputFolderName}/tmp_${outName}/start_alignment_mapping/${contigName}_afterStartAlignment_unfiltered.log
  awk -F "\t" '{OFS="\t"}{if($12>50) print $0}' ${outputFolderName}/tmp_${outName}/start_alignment_mapping/${contigName}_afterStartAlignment_unfiltered.minimap > ${outputFolderName}/tmp_${outName}/start_alignment_mapping/${contigName}_afterStartAlignment.minimap

  ##--------------------------------------------minimap2 analysis------------------------------------


  #find out if there are mappings and in which orientation
  number_of_mappings=$(wc -l  ${outputFolderName}/tmp_${outName}/start_alignment_mapping/${contigName}_afterStartAlignment.minimap|cut -d ' ' -f 1)
  orientation=$(sort -k10 -n -r  ${outputFolderName}/tmp_${outName}/start_alignment_mapping/${contigName}_afterStartAlignment.minimap|head -1 |cut -f 5)


  # create a info table containing information about the presence of dnaA and if present the orientation

  if [ $number_of_mappings == "0" ]
  then
    echo -e "           "$contigName"......no dnaA map found"
        contig_length=$(grep -v ">" ${outputFolderName}/tmp_${outName}/genome/${contigName}.fasta|wc -c )

         echo -e ${contigName}"\tnon-dnaA-containing-contig\t"${number_of_mappings}"\tNA\tNA\t"${contig_length} >>  ${outputFolderName}/tmp_${outName}/start_alignment_mapping/After_StartAlignment_contigs.minimap
         echo "no dnaA map found"
         echo "SOMETHING IS WRONG BECAUSE dnaA cannot be found anymore...aborting"
         exit

  elif [ $orientation == "+" ]
  then

          start=$(sort -k8 -n ${outputFolderName}/tmp_${outName}/start_alignment_mapping/${contigName}_afterStartAlignment.minimap|head -1 |awk -F "\t" '{OFS="\t"}{print $8-$3-5}')
          contig_length=$(sort -k8 -n ${outputFolderName}/tmp_${outName}/start_alignment_mapping/${contigName}_afterStartAlignment.minimap|head -1 |cut -f 7)

          echo -e ${contigName}"\tdnaA-containing-contig\t"${number_of_mappings}"\t"${orientation}"\t"${start}"\t"${contig_length} >>  ${outputFolderName}/tmp_${outName}/start_alignment_mapping/After_StartAlignment_contigs.minimap
          #echo "Looks good!"
          echo -e "           "$contigName"......dnaA map found on "${orientation}" strand at position "${start}

  else
        echo -e "           "$contigName"......dnaA map found"


      start=$(sort -k9 -n -r ${outputFolderName}/tmp_${outName}/start_alignment_mapping/${contigName}_afterStartAlignment.minimap|head -1 |awk -F "\t" '{OFS="\t"}{print $9+$2-$4+5}')
          contig_length=$(sort -k9 -n -r  ${outputFolderName}/tmp_${outName}/start_alignment_mapping/${contigName}_afterStartAlignment.minimap|head -1 |cut -f 7)

          echo -e ${contigName}"\tdnaA-containing-contig\t"${number_of_mappings}"\t"${orientation}"\t"${start}"\t"${contig_length} >> ${outputFolderName}/tmp_${outName}/start_alignment_mapping/After_StartAlignment_contigs.minimap
          echo "SOMETHING IS WRONG BECAUSE dnaA IS STILL ON THE REVERSE STRAND...aborting"
          exit
  fi

  done

  echo

  ###########################################################
  #Bring together
  ############################################################
  echo -e "Bring contigs together"


  ##--------------------------------------------remove final fasta if present----------------------------------------

  if [ -f "${outputFolderName}/${outName}.fasta" ]; then
  rm  ${outputFolderName}/${outName}.fasta
fi


for header in $(grep ">" ${outputFolderName}/tmp_${outName}/genome/tmp_wide_all.fasta |sed 's/>//g')
  do

  circular=$(awk -F "\t" -v contigzz="$header" '{OFS="\t"}{if($1==contigzz)print $10} ' ${outputFolderName}/analysis_circularity_extended.log )


   if [ $circular == "Y" ];
      then
       echo -e "              Circularized contig name......" ${header}
       cat ${outputFolderName}/tmp_${outName}/StartAlignedContigs/${contigName}_startAligned.fasta  >> ${outputFolderName}/${outName}.fasta

    else
      echo -e "              Non-Circularized contig name......" ${header}

       cat  ${outputFolderName}/tmp_${outName}/genome/${header}.fasta  >> ${outputFolderName}/${outName}.fasta


    fi

  done


  ###########################################################
  #Final Infos
  ############################################################
  echo -e "Final Infos"


  numberContigs=$(grep "^#" -c ${outputFolderName}/analysis_circularity_extended.log)
  BacContigs=$(grep "^#" -v  ${outputFolderName}/analysis_circularity_extended.log |awk -F "\t" '{OFS="\t"}{if($2=="dnaA-containing-contig")print $0}'|wc -l)
  NonBacContigs=$(grep "^#" -v  ${outputFolderName}/analysis_circularity_extended.log |awk -F "\t" '{OFS="\t"}{if($2!="dnaA-containing-contig")print $0}' |wc -l)

  echo -e "            Total number of contigs......"${numberContigs}
  echo -e "            Number of startaligned contigs contigs......"${BacContigs}
  #grep "^#" -v  ${outputFolderName}/start_alignment_mapping/StartAlignment_contigs.minimap |awk -F "\t" '{OFS="\t"}{if($2=="bacterial_contig")print $0}'
  echo -e "            Number of non-startaligned contigs......"${NonBacContigs}
  #grep "^#" -v  ${outputFolderName}/start_alignment_mapping/StartAlignment_contigs.minimap |awk -F "\t" '{OFS="\t"}{if($2!="bacterial_contig")print $0}'
echo

  echo -e "The location of dnaA on the "${BacContigs}" Bacterial contigs..."
  #grep "^#"  ${outputFolderName}/start_alignment_mapping/After_StartAlignment_contigs.minimap
  #grep "^#" -v  ${outputFolderName}/start_alignment_mapping/After_StartAlignment_contigs.minimap|awk -F "\t" '{OFS="\t"}{if($2=="bacterial_contig")print $0}'
  cat ${outputFolderName}/tmp_${outName}/start_alignment_mapping/After_StartAlignment_contigs.minimap
echo
  echo -e "Output"

  echo -e "      Location of Startaligned contigs......"${outputFolderName}/${outName}".fasta"
  echo -e "      Location of log file......"${outputFolderName}"/analysis_circularity_extended.log"


  runtime=$((($(date +%s)-$start)/60))
  echo -e "      The script ran for......"${runtime}" minutes"

echo

##--------------------------------------------remove tmp files----------------------------------------

if [ "$tmpss" == "N" ]
then
    rm -r ${outputFolderName}/tmp_${outName}
    echo -e "      removing all tmps files......Done"
fi
echo -e "      See you soon!"
echo
