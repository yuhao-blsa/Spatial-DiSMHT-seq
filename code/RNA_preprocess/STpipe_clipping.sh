#!/bin/bash
# 1
core=8
mem=100
MAP=/path/to/STAR_index
ANN=/path/to/Homo_sapiens.GRCh38.113.chrname.transcript_as_exon.gtf
CONT=/path/to/STAR_index
gene_model=/path/to/hg38_GENCODE_V47.bed
R_path=/path/to/Rscript

# 2 Variable parameter
batch=exp4_250407
slice=sample
directory=/path/to/data
# 3 Variable parameter
ID_File=/path/to/combine_barcode.round2round1_index1_index2.192.txt
RNA_R1=/path/to/R1.fq.gz
RNA_R2=/path/to/R2.fq.gz


# Fixed combination
name=RNA_${batch}_${slice}
out_path=${directory}/scRNA_${slice}

sh ./1_extract_barcode.sh $name $RNA_R1 $RNA_R2 $out_path/fastq
$R_path ./2_check_barcode.R -b $ID_File -p $out_path/fastq -n ${name}_2.extract.barcode
rm $out_path/fastq/${name}_2.extract.barcode
sh ./3_st_pipeline.soft_clipping.sh $core $mem $name $out_path \
${out_path}/fastq/${name}_2.extract.fq.gz ${out_path}/fastq/${name}_1.extract.fq.gz $MAP $ANN $CONT $ID_File $gene_model
rm ${out_path}/1_stpipeline/tmp/R2_quality_trimmed.bam ${out_path}/1_stpipeline/tmp/contaminated_clean.bam 
