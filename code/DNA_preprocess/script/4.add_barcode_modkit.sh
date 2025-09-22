#!/bin/bash
config_file=$1
source $config_file
mkdir -p ${outpath}/4_modkit_threshold
mkdir -p ${outpath}/5_Plot

${samtools_path} view -h -b -@ $threads -F 2308 ${outpath}/1_dorado_out/1_calls.all.mapped.${sample_id}.bam > ${outpath}/1_dorado_out/2_calls.all.mapped.${sample_id}_filtered.bam
${samtools_path} sort -@ $threads -o ${outpath}/1_dorado_out/2_calls.all.mapped.${sample_id}_filtered.sorted.bam ${outpath}/1_dorado_out/2_calls.all.mapped.${sample_id}_filtered.bam
${samtools_path} index ${outpath}/1_dorado_out/2_calls.all.mapped.${sample_id}_filtered.sorted.bam

awk '{print $1"\tCB\t"$5}' ${outpath}/3_corrected/3_out_id2coord.${sample_id}.tsv > ${outpath}/3_corrected/3_out_id2coord.sinto.${sample_id}.tsv
awk '{print $1"\tSP\t"$2}' ${outpath}/3_corrected/3_out_id2coord.${sample_id}.tsv >> ${outpath}/3_corrected/3_out_id2coord.sinto.${sample_id}.tsv
${sinto_path} addtags -m readname -b ${outpath}/1_dorado_out/2_calls.all.mapped.${sample_id}_filtered.sorted.bam -f ${outpath}/3_corrected/3_out_id2coord.sinto.${sample_id}.tsv -o ${outpath}/1_dorado_out/2_calls.all.mapped.${sample_id}_filtered.sorted.CB_tag.bam -p 20

${samtools_path} index ${outpath}/1_dorado_out/2_calls.all.mapped.${sample_id}_filtered.sorted.CB_tag.bam
printf "${sample_id} dorado barcode: %s\n" "$(${samtools_path} view -F 2308 ${outpath}/1_dorado_out/2_calls.all.mapped.${sample_id}_filtered.sorted.CB_tag.bam | grep "CB:Z:" | wc -l)" >> ${outpath}/Log.txt 

${bamCoverage_path} -b ${outpath}/1_dorado_out/2_calls.all.mapped.${sample_id}_filtered.sorted.CB_tag.bam -o ${outpath}/1_dorado_out/2_calls.all.mapped.${sample_id}_filtered.sorted.CB_tag.None.bw --binSize ${binsize} --normalizeUsing None --effectiveGenomeSize ${effectiveGenomeSize} --numberOfProcessors 20

# Split BAM file according to sample
awk '{print $2"\t"$2}' ${outpath}/3_corrected/3_out_id2coord.${sample_id}.tsv | sort | uniq > ${outpath}/3_corrected/3_out_id2coord.${sample_id}.sample.tsv
sinto filterbarcodes -b ${outpath}/1_dorado_out/2_calls.all.mapped.${sample_id}_filtered.sorted.CB_tag.bam -c ${outpath}/3_corrected/3_out_id2coord.${sample_id}.sample.tsv --barcodetag SP --outdir ${outpath}/1_dorado_out/3_split_sample_${sample_id} -p $threads
for sample_i in ${sample_ls}
do
${samtools_path} index ${outpath}/1_dorado_out/3_split_sample_${sample_id}/${sample_i}.bam
${mosdepth_path} -n -t $threads ${outpath}/1_dorado_out/3_split_sample_${sample_id}/1.${sample_i}.depth_out /${outpath}1_dorado_out/3_split_sample_${sample_id}/${sample_i}.bam
${mosdepth_path} -n -t $threads --by ${CpGbed} ${outpath}/1_dorado_out/3_split_sample_${sample_id}/2.${sample_i}.CpG.depth_out ${outpath}/1_dorado_out/3_split_sample_${sample_id}/${sample_i}.bam
done

# 2 repeat element
for file in ${SINE_file} ${LINE_file} ${LTR_file} ${Satellite_file} ${CGI_file}
do
cut -f-3 $file |${bedtools_path} coverage -a - -b ${outpath}/1_dorado_out/2_calls.all.mapped.${sample_id}_filtered.sorted.CB_tag.bam |awk -v file="$file" -v sample_id="$sample_id" '{total_covered += $5; total_length += $6} END {print sample_id, file, total_covered / total_length}' >> ${outpath}/Log.txt
done

# 3 modkit
## get a tab-separated file of threshold values for each modified base
${modkit_path} sample-probs ${outpath}/1_dorado_out/2_calls.all.mapped.${sample_id}_filtered.sorted.CB_tag.bam --hist --out-dir ${outpath}/4_modkit_threshold/${sample_id} -t $threads
ulimit -n 80000 # increase open file limit
mkdir -p ${outpath}/4_modkit_perpixel_CpG
# only CpG

${modkit_path} pileup ${outpath}/1_dorado_out/2_calls.all.mapped.${sample_id}_filtered.sorted.CB_tag.bam ${outpath}/4_modkit_perpixel_CpG/${sample_id}_soft_threshold --cpg --ref ${reference_file} \
--log-filepath ${outpath}/4_modkit_perpixel_CpG/${sample_id}_soft_threshold/pileup.log --partition-tag CB --prefix split --threads $threads --filter-threshold C:0.7 --mod-thresholds m:0.7 --mod-thresholds h:0.6
find ${outpath}/4_modkit_perpixel_CpG/${sample_id}_soft_threshold/ -type f ! -size 0 -name "split_i*.bed" | sort > ${outpath}/4_modkit_perpixel_CpG/pixel_id.${sample_id}_soft_threshold.txt

mkdir -p ${outpath}/4_modkit_perpixel_nonCpG
${modkit_path} pileup ${outpath}/1_dorado_out/2_calls.all.mapped.${sample_id}_filtered.sorted.CB_tag.bam ${outpath}/4_modkit_perpixel_nonCpG/${sample_id}_CHG_soft_threshold \
--motif CAG 0 --motif CCG 0 --motif CTG 0 --ref ${reference_file} \
--log-filepath ${outpath}/4_modkit_perpixel_nonCpG/${sample_id}_CHG_soft_threshold/pileup.log --partition-tag CB --prefix split --threads $threads \
--filter-threshold C:0.7 --mod-thresholds m:0.7 --mod-thresholds h:0.6
find ${outpath}/4_modkit_perpixel_nonCpG/${sample_id}_CHG_soft_threshold/ -type f ! -size 0 -name "split_i*.bed" | sort > ${outpath}/4_modkit_perpixel_nonCpG/pixel_id.${sample_id}_CHG_soft_threshold.txt

${modkit_path} pileup ${outpath}/1_dorado_out/2_calls.all.mapped.${sample_id}_filtered.sorted.CB_tag.bam ${outpath}/4_modkit_perpixel_nonCpG/${sample_id}_CHH_soft_threshold \
--motif CAA 0 --motif CCA 0 --motif CTA 0 --motif CAC 0 --motif CCC 0 --motif CTC 0 --motif CAT 0 --motif CCT 0 --motif CTT 0 --ref ${reference_file} \
--log-filepath ${outpath}/4_modkit_perpixel_nonCpG/${sample_id}_CHH_soft_threshold/pileup.log --partition-tag CB --prefix split --threads $threads \
--filter-threshold C:0.7 --mod-thresholds m:0.7 --mod-thresholds h:0.6 &&
find ${outpath}/4_modkit_perpixel_nonCpG/${sample_id}_CHH_soft_threshold/ -type f ! -size 0 -name "split_i*.bed" | sort > ${outpath}/4_modkit_perpixel_nonCpG/pixel_id.${sample_id}_CHH_soft_threshold.txt &&

# 4 calculate mapped ref length
${samtools_path} view ${outpath}/1_dorado_out/2_calls.all.mapped.${sample_id}_filtered.sorted.CB_tag.bam | awk '$0 ~ /CB:Z:/{
    match($0, /CB:Z:([^\t]+)/);
    cb_value = substr($0, RSTART+5, RLENGTH-5);
    cigar = $6;
    total_length = 0;
    while (match(cigar, /([0-9]+)([MIDSHN])/)) {
        leni = substr(cigar, RSTART, RLENGTH-1);
        type = substr(cigar, RSTART+RLENGTH-1, 1);
        if (type == "M" || type == "D" || type == "N") {
            total_length += leni;
        }
        cigar = substr(cigar, RSTART + RLENGTH);
    }
    print cb_value, total_length
}' > ${outpath}/5_Plot/1_mapped_ref_length.${sample_id}.tsv &&

${samtools_path} view ${outpath}/1_dorado_out/2_calls.all.mapped.${sample_id}_filtered.sorted.CB_tag.bam | awk '$0 ~ /CB:Z:/{
    match($0, /CB:Z:([^\t]+)/);
    cb_value = substr($0, RSTART+5, RLENGTH-5);
    match($0, /MN:i:([^\t]+)/);
    seq_len = substr($0, RSTART+5, RLENGTH-5);
    print cb_value, seq_len
}' > ${outpath}/5_Plot/1_sequencing_length.${sample_id}.tsv &&

awk '{a[$1]+=$2} END {for (i in a) print i, a[i]}' ${outpath}/5_Plot/1_mapped_ref_length.${sample_id}.tsv > ${outpath}/5_Plot/1_mapped_ref_length.summary.${sample_id}.tsv &&

[ -s "${outpath}/1_dorado_out/2_calls.all.mapped.${sample_id}_filtered.sorted.bam" ] &&\
        [ -s "${outpath}/4_modkit_perpixel_nonCpG/pixel_id.${sample_id}_CHH_soft_threshold.txt" ] &&\
	[ -s "${outpath}/4_modkit_perpixel_CpG/pixel_id.${sample_id}_soft_threshold.txt" ] &&\
        echo "Succeed" >${outpath}/log/4.add_barcode_modkit.log
