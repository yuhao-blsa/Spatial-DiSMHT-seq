#!/bin/bash
config_file=$1
source $config_file
mkdir -p ${outpath}/1_dorado_out
${dorado_path} basecaller hac,5mC_5hmC ${input_pod}  --reference ${reference_file} > ${outpath}/1_dorado_out/1_calls.all.mapped.${sample_id}.bam &&

${samtools_path} fastq ${outpath}/1_dorado_out/1_calls.all.mapped.${sample_id}.bam -@ $threads | gzip >${outpath}/1_dorado_out/1_calls.all.mapped.${sample_id}.fq.gz 

printf "%s dorado total: %s\n" "$sample_id" "$(${samtools_path} view -c -F 2304 ${outpath}/1_dorado_out/1_calls.all.mapped.${sample_id}.bam)" >> ${outpath}/Log.txt  # total 
printf "%s dorado mapped: %s\n" "$sample_id" "$(${samtools_path} view -c -F 2308 ${outpath}/1_dorado_out/1_calls.all.mapped.${sample_id}.bam)" >> ${outpath}/Log.txt # mapped

[ -s "${outpath}/1_dorado_out/1_calls.all.mapped.${sample_id}.bam" ] && [ -s "${outpath}/1_dorado_out/1_calls.all.mapped.${sample_id}.fq.gz" ]  && echo "Succeed" >${outpath}/log/1.dorado.log

