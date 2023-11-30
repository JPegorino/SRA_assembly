#!/bin/bash
# This is a combined job script to download reads, perform QC/adapter trimming and assemble genomes for SRA run accessions in EDDIE 
# The script takes a singled input file: comma-delimited table in the format output from the SRA run selector.
# The prererquisites for the script can be installed in a conda environment called 'sra' using the line below: 
# `mamba install -c bioconda -c conda-forge -c defaults fastqc=0.11.7 trimmomatic=0.36 spades=3.15.4 kraken2=2.0.9beta sra-tools=2.10.0 entrez-direct=13.3`

### Setting qsub parameters
#$ -cwd -l h_rt=05:00:00 -l h_vmem=16G
#$ -t 2-500:1 ### REPLACE 500 with the number of lines in the SraRunTable.csv file

### load software
source /etc/profile.d/modules.sh
module load anaconda/5.3.1
source activate sra # all software already loaded in this environment

# prepare miscellaneous input 
scratch_dir="/exports/eddie/scratch/$(whoami)"
kraken_db_path="/exports/cmvm/eddie/eb/groups/fitzgerald_grp/software/kraken_default_database/kraken_ddb" # location of a KRAKEN database
trimmomatic_adapter_path=$(find ${CONDA_PREFIX} -name "adapters" | xargs -I{} readlink -m {} | head -1) # path to trimmomatic adapters
cat ${trimmomatic_adapter_path}/*.fa > trimmomatic_adapters.fa # create a fasta file in cwd containing the combined adapter sequences included with installed trimmomatic version
params="ILLUMINACLIP:./trimmomatic_adapters.fa:2:30:10 HEADCROP:6 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:50" #trimmomatic parameters

### Setting command variables
infile=$1

### Main Script -- add the commands for the script you want to run in eddie ###################################################

# Troubleshooting options: print confirmation of key settings in output #
echo "This is the job ${JOB_NAME} with the unique ID: ${JOB_ID}"
echo "It is number ${SGE_TASK_ID} in an array job from ${SGE_TASK_FIRST} to ${SGE_TASK_LAST} by steps of ${SGE_TASK_STEPSIZE}"
echo "This script takes one input file, either a txt list or csv table with run information from NCBI. The file must have run accessions in column 1 (columns delimited by ',')"
if [ -n "${infile}" ] ; then echo "It was run for file ${infile} line ${SGE_TASK_ID} " ; else echo "ERROR: no input file" ; exit 1 ; fi

### Writing functions -- define any functions called later in the script
function qc_stats() { # function to calculate basic QC stats for an NCBI SRA fastq input file, genome size and sequencer
 local file_name=`basename "${1}"`
 local xsn=`basename "${1}" | grep -o "[SED]RR[0-9]*"`
 local readcount=`cat "${1}" | grep -E '^@[EDS]RR' | grep -Ec ' length='`
 local averagelength=`cat "${1}" | grep ' length=' | awk '{print $NF}' | cut -d '=' -f 2 | awk '{ total += $1 } END { print total/NR }'`
 local filesize=`du "${1}" | awk '{print $1}'`
 local depth=`awk -v v1=${averagelength} -v v2=${readcount} -v v3=${2} "BEGIN {printf \"%.2f\", v1 * v2 * v3}"`
 local approx_coverage=`echo "scale=2; ${depth} / 3000000" | bc -l`
 local filestats=`echo "$file_name" "$xsn" "$filesize" "$readcount" "$averagelength" "$depth" "$approx_coverage"`
 echo "$filestats" ; }
 # a header line for the qc_stats files
header_line=`echo {'FILE','RUN','FILE_SIZE_BYTES','READ_COUNT','MEAN_READ_LENGTH','SEQUENCING_DEPTH','COVERAGE_ESTIMATE'} | sed 's/ /,/g'` # make a header line for output stats tables

### Main Script (main)
xsn=$(head -"${SGE_TASK_ID}" ${infile} | tail -1 | cut -d, -f1) # the run accession
mkdir "./${xsn}" # make an empty directory for the accession
echo "downloading reads for ${xsn} ..."
  # Now, we'll keep attempting to download reads for each sequence until:
attempt_count=1 ; while [ $(ls -A "./${xsn}" | wc -l) -eq 0 ] && [ $attempt_count -lt 11 ] # either 10 downloads have been attempted OR output files exist
  do echo "download attempt $((attempt_count++)) ..."
  fasterq-dump -O "./${xsn}" -p -f --split-files "${xsn}" ; done # download reads as fastq files for accession - use -f to force it to overwrite empty existing directory for xsn
echo "organising data for ${xsn} reads ..."
# make variables (strings) formatted for each of the potentially expected output files
mate0=$(echo "./${xsn}/${xsn}.fastq") # fastq file for single-sequence
mate1=$(echo "./${xsn}/${xsn}_1.fastq") # fastq file for paired-sequence mate 1
mate2=$(echo "./${xsn}/${xsn}_2.fastq") # fastq file for paired-sequence mate 2
fastq_count=$(ls ./${xsn}/*.fastq | wc -l) # the number of downloaded fastq files (1 if single, 2 if paired)
# create the QC stats files and run basic qc for untrimmed reads
echo "${header_line}" > ${xsn}/untrimmed_QC_FASTQ_stats.csv 
echo "${header_line}" > ${xsn}/trimmed_QC_FASTQ_stats.csv
for fastq_file in ${xsn}/${xsn}*.fastq
  do echo "collecting stats for $xsn reads ..."
  qc_stats "${fastq_file}" ${fastq_count} | sed 's/ /,/g' >> ${xsn}/untrimmed_QC_FASTQ_stats.csv 
  echo "running fastqc for raw $xsn reads ..."
  fastqc --outdir="${xsn}" "${fastq_file}" ; done # add stats for each .fastq file to the .csv file
# check the formatting of any downloaded files
if [[ -f "${mate1}" && -f "${mate2}" ]] # if downloaded read files exist files exist for each mate in a pair (with formatting to indicate they are paired end)
  then echo "trimming $xsn reads ..."
  trimmomatic PE -threads 1 -phred33 "${mate1}" "${mate2}" "${xsn}/${xsn}_tP1.fastq" "${xsn}/${xsn}_tU1.fastq" "${xsn}/${xsn}_tP2.fastq" "${xsn}/${xsn}_tU2.fastq" ${params}
  cat "${xsn}/${xsn}_tU1.fastq" "${xsn}/${xsn}_tU2.fastq" > "${xsn}/${xsn}_tS.fastq" # combined orphaned reads for SPAdes
  for trimmed_fastq_file in ${xsn}/${xsn}_t* # collect stats and run fastqc for all the read files
    do echo "collecting trimmed stats for $xsn from $trimmed_fastq_file ..."
    qc_stats "${trimmed_fastq_file}" ${fastq_count} | sed 's/ /,/g' >> ${xsn}/trimmed_QC_FASTQ_stats.csv
    echo "running fastqc for trimmed reads in ${trimmed_fastq_file} ..."
    fastqc --outdir="${xsn}" "${trimmed_fastq_file}" ; done # Run FastQC on trimmed files
  # run trimmed reads against the specified kraken database
  echo "running kraken for $xsn reads against the database in $kraken_db_path"
  kraken2 --threads 1 --paired --db "${kraken_db_path}" --output "${xsn}/K_ddb_${xsn}.txt" --report "${xsn}/k_ddb_${xsn}_report.txt" "${mate1/_1/_tP1}" "${mate2/_2/_tP2}"
  echo "assembling $xsn reads with SPAdes ..."
  spades.py -t1 --careful --phred-offset 33 -o "${xsn}/assembly" -1 "${mate1/_1/_tP1}" -2 "${mate2/_2/_tP2}" -s "${mate2/_2/_tS}"
    if [ $(ls -A "./${xsn}/assembly" | wc -l) -eq 0 ] ; then
    echo "SPAdes was unsuccessful ..."
    final_outdir="assembly_failed"
  else
    mv -nv "${xsn}"/assembly/*.fasta* "${xsn}" && rm -r "${xsn}/assembly" # keep only fasta output files and assembly graph
    final_outdir="completed"
  fi
  rm -v "${mate1}" "${mate2}" # deletes raw read files
elif [[ -f "${mate0}" ]] # otherwise, if there is a single downloaded read file (with formatting to indicate single-end reads)
  then echo "trimming $xsn reads ..."
  trimmomatic SE -threads 1 -phred33 ${mate0} "${xsn}/${xsn}_t.fastq" ${params}
  echo "collecting trimmed stats for $xsn from ${xsn}_t.fastq ..."
  qc_stats "${xsn}/${xsn}_t.fastq" ${fastq_count} | sed 's/ /,/g' >> trimmed_QC_FASTQ_stats.csv
  echo "running post-trim fastqc for $out_readfile ..."
  fastqc --outdir="${xsn}" "${xsn}/${xsn}_t.fastq" # Run FastQC on trimmed files
  # run trimmed reads against the specified kraken database
  echo "running kraken for $xsn reads against the database in $kraken_db_path ..."
  kraken2 --threads 1 --db "${kraken_db_path}" --use-names --output "${xsn}/K_ddb_${xsn}.txt" --report "${xsn}/k_ddb_${xsn}_report.txt" "${xsn}/${xsn}_t.fastq"
  echo "assembling $xsn reads with SPAdes ..."
  spades.py -t1 --careful --phred-offset 33 -o "${xsn}/assembly" -s "${mate2/_2/_tS}"
  if [ $(ls -A "./${xsn}/assembly" | wc -l) -eq 0 ] ; then
    echo "SPAdes was unsuccessful ..."
    final_outdir="assembly_failed"
  else
    mv -nv "${xsn}"/assembly/*.fasta* "${xsn}" && rm -r "${xsn}/assembly" # keep only fasta output files and assembly graph
    final_outdir="completed" # indicate an output directory for single-end reads
  fi
  rm -v "${mate0}" # deletes raw read file
else echo "downloading and processing unsuccessful for ${xsn} ..." && final_outdir="./repeat_downloads" ; fi # in any other scenario, (i.e. if downloaded read files have unexpected formatting), flag for repeat

### Post-analysis clean-up
echo "cleaninig up ..."
gzip ./${xsn}/*.fastq # compress completed trimmed and untrimmed read files
cd ${xsn} # hop into the output directory and rename the SPAdes output
for fasta_file in *.fasta
  do mv -nv "${fasta_file}" "${xsn}_${fasta_file}"
done
rm -v *before_rr*.fasta # only the final assemblies are useful
cd ..
if [[ ! -e "${final_outdir}" ]] ; then mkdir "${final_outdir}" ; fi # make a directory for each platform (without one already)
mv -nv "${xsn}" "${final_outdir}" # move directories for successful downloads to a new location based on the output
rm ./SRA/sra/${xsn}* # remove any cache files (if this isn't done automatically)
############################################

### Results Summary
echo "analysis complete for ${xsn}, generating ${xsn}_README.txt:"
echo -e "$(ls ${final_outdir}/${xsn}/*[12].fastq.gz | wc -l) raw read files:\n$(ls ${final_outdir}/${xsn}/*[12].fastq.gz)" >> ${xsn}_README.txt
echo -e "$(ls ${final_outdir}/${xsn}/*t[SUP]*.fastq.gz | wc -l) trimmed read files:\n$(ls ${final_outdir}/${xsn}/*t[SUP]*.fastq.gz)" >> ${xsn}_README.txt
echo -e "$(ls ${final_outdir}/${xsn}/${xsn}_*.fasta | wc -l) assembly files:\n$(ls ${final_outdir}/${xsn}/${xsn}_*.fasta)" >> ${xsn}_README.txt
echo -e "$(ls ${final_outdir}/${xsn}/*fastqc* | wc -l) fastqc files:\n$(ls ${final_outdir}/${xsn}/*fastqc*)" >> ${xsn}_README.txt
echo -e "$(ls ${final_outdir}/${xsn}/*.csv | wc -l) read QC stats tables:\n$(ls ${final_outdir}/${xsn}/*.csv)" >> ${xsn}_README.txt
echo -e "$(ls ${final_outdir}/${xsn}/[Kk]_ddb_*.txt | wc -l) kraken output files:\nReport $(ls ${final_outdir}/${xsn}/K_ddb*.txt)\nSummary: $(ls ${final_outdir}/${xsn}/k_ddb*.txt)" >> ${xsn}_README.txt
mv ${xsn}_README.txt ${final_outdir}/${xsn}/README.txt

conda deactivate # exit conda environment - restore defaults
echo "script complete"
### end of script ########################################################################################################
