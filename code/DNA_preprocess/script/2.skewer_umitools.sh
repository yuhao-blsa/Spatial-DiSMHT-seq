#!/bin/bash
config_file=$1
source $config_file
mkdir -p ${outpath}/2_skewer_trim

${skewer_path} -f Sanger -t $threads -m head -l ${skewer_min_read_len} -x ${skewer_adapter} -z -o ${outpath}/2_skewer_trim/1_skewer_trim.${sample_id} ${outpath}/1_dorado_out/1_calls.all.mapped.${sample_id}.fq.gz &&

${umitools_path} extract --extract-method=regex --bc-pattern="${bcpattern}" -I ${outpath}/2_skewer_trim/1_skewer_trim.${sample_id}-trimmed.fastq.gz -S ${outpath}/2_skewer_trim/2_umi_extract.${sample_id}.fq.gz -L ${outpath}/2_skewer_trim/2_umi_extract.${sample_id}.log &&

zcat ${outpath}/2_skewer_trim/2_umi_extract.${sample_id}.fq.gz | ${seqkit_path} subseq -r 1:24 | gzip -c - > ${outpath}/2_skewer_trim/3.I7barcode2barcode1.${sample_id}.fq.gz 
zcat ${outpath}/2_skewer_trim/2_umi_extract.${sample_id}.fq.gz | ${seqkit_path} subseq -r 1:8 | gzip -c - > ${outpath}/2_skewer_trim/3.I7.${sample_id}.fq.gz 
zcat ${outpath}/2_skewer_trim/2_umi_extract.${sample_id}.fq.gz | ${seqkit_path} subseq -r 9:24 | gzip -c - > ${outpath}/2_skewer_trim/3.barcode2barcode1.${sample_id}.fq.gz 
zcat ${outpath}/2_skewer_trim/2_umi_extract.${sample_id}.fq.gz | ${seqkit_path} subseq -r 25:-1 | gzip -c - > ${outpath}/2_skewer_trim/3.read.${sample_id}.fq.gz 

[ -s "${outpath}/2_skewer_trim/3.I7barcode2barcode1.${sample_id}.fq.gz" ] &&\
       	[ -s "${outpath}/2_skewer_trim/1_skewer_trim.${sample_id}-trimmed.fastq.gz" ] &&\
       	[ -s "${outpath}/2_skewer_trim/2_umi_extract.${sample_id}.fq.gz" ] &&\
       	echo "Succeed" >${outpath}/log/2.skewer.trim.log
