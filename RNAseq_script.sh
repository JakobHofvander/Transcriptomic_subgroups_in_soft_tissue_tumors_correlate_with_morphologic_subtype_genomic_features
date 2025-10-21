#!/bin/sh
#SBATCH -n 40
#SBATCH -N 1
#SBATCH -t 80:00:00
#SBATCH -A csens2024-3-4
#SBATCH -J STAR
#SBATCH -o STAR.%j.out
#SBATCH -e STARSingSing.%j.err
module purge

module load GCC/11.3.0
module load STAR/2.7.10b
module load GCCcore/11.3.0
module load parallel/20220722

export RUN="Path/to/Run/Directory"
export WD="Path/to/Working/Directory"
export FQ="Path/to/fastq/Directory"

parallel -v -j4 -k --progress -N2 "STAR --genomeDir /scale/gr02/shared/mertens-lab/lsens_fs4/References/STAR_overhang_75/ \
--readFilesCommand zcat \
      --readFilesIn {1} {2} \
      --twopassMode Basic \
      --chimSegmentMin 12 \
      --chimJunctionOverhangMin 12 \
      --alignSJDBoverhangMin 10 \
      --alignMatesGapMax 200000 \
      --alignIntronMax 200000 \
      --chimSegmentReadGapMax parameter 3 \
      --alignSJstitchMismatchNmax 5 -1 5 5 \
      --runThreadN 10 \
      --limitBAMsortRAM 31532137230 \
      --outSAMtype BAM SortedByCoordinate \
      --chimOutType WithinBAM \
      --outSAMstrandField intronMotif \
      --quantMode TranscriptomeSAM GeneCounts \
      --quantTranscriptomeBan IndelSoftclipSingleend \
--outFileNamePrefix $WD/$RUN/STAR-bam-hg38_CSENS/{1/.}_ \
--outFilterMultimapNmax 200" ::: $FQ/$RUN/*-*/*fastq.gz

##############################
### Strandness calculation ###
##############################
# Creating RSEM strandness and expression output directory in RUN folder
mkdir RSEM_strandness RSEM_expression RSEM_expression/Temp

#### strandedness calculation for RSEM
for file in $WD/$RUN/STAR-bam-hg38/*-*/*_ReadsPerGene.out.tab
do
fname=`echo $file | sed 's/ReadsPerGene.out.tab/Aligned.toTranscriptome.out.bam/g'`
samp=`basename $file | cut -d_ -f1,2`
bdir=`dirname $file | cut -d/ -f1,2,3,4,5,6`
Actual=`grep -v "N_" $file | awk -v OFS='\t' '{unst+=$2;forw+=$3;rev+=$4}END{print unst,forw,rev,100*(rev/(rev+forw))}'`
echo -e "$samp\t$bdir\t$fname\t$Actual" >> $WD/$RUN/RSEM_strandness/strandedness_temp.txt
done
#
awk -v OFS='\t' '{
if ($7 >=40 && $7 <= 55) print $0,"none";
else if ($7 >=0 && $7 <= 15) print $0,"forward";
else if ($7 >=85) print $0,"reverse";
else print $0,"error";
}' $WD/$RUN/RSEM_strandness/strandedness_temp.txt >> $WD/$RUN/RSEM_strandness/strandedness.txt
rm -rf $WD/$RUN/RSEM_strandness/strandedness_temp.txt

## check strandness file
cat $WD/$RUN/RSEM_strandness/strandedness.txt
grep none RSEM_strandness/strandedness.txt
grep reverse RSEM_strandness/strandedness.txt

###### RSEM
#!/bin/sh
#SBATCH -n 40
#SBATCH -N 1
#SBATCH -t 80:00:00
#SBATCH -A csens2024-3-4
#SBATCH -J RSEM
#SBATCH -o RSEM.%j.out
#SBATCH -e RSEMSingSing.%j.err

module purge
module load GCC/11.3.0 OpenMPI/4.1.4 RSEM/1.3.3
module load GCCcore/11.3.0 parallel/20220722


export RUN="Path/to/Run/Directory"
export WD="Path/to/Working/Directory"
export FQ="Path/to/fastq/Directory"

parallel -v -j4 -k "rsem-calculate-expression --bam \
                                                    --paired-end \
                                                    --num-threads 10 \
                                                    --no-bam-output \
                                                    --forward-prob 0 \
                                                    --time \
                                                    --estimate-rspd \
                                                    --temporary-folder $WD/$RUN/RSEM_expression/Temp/{1/.} \
                                                    {1} \
                                                    /scale/gr02/shared/mertens-lab/lsens_fs4/References/Ensembl/Homo_sapiens/hg38/RSEM_star/ \
                                                    $WD/$RUN/RSEM_expression/{1/.}" ::: $WD/$RUN/STAR-bam-hg38/*-*/*_Aligned.toTranscriptome.out.bam


####### run awk to generate seperate files for tpm,fpkm, counts
for i in *toTranscriptome.out.genes.results; do
awk '{print $1"\t"$5}' $i > $i.temp;
awk 'NR==1{$1="transcript_id""\t";$2="'${i%%.*}'"}1' $i.temp > $i.expected_count.table;
awk '{print $1"\t"$7}' $i > $i.temp2;
awk 'NR==1{$1="transcript_id""\t";$2="'${i%%.*}'"}1' $i.temp2 > $i.fpkm_count.table;
awk '{print $1"\t"$6}' $i > $i.temp3;
awk 'NR==1{$1="transcript_id""\t";$2="'${i%%.*}'"}1' $i.temp3 > $i.tpm_count.table;
rm $i.tem*;
done
