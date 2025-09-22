config=$1
source $config
mkdir -p ${outpath}/log

#1.dorado
if [ ! -f "${outpath}/log/1.dorado.log" ];then
        bash ./script/1.dorado.sh $config
fi
##2.skewer
if [ ! -f "${outpath}/log/2.skewer.trim.log" ]&&[ -f "${outpath}/log/1.dorado.log" ];then
	bash ./script/2.skewer_umitools.sh $config
fi
##3.correct_barcode
if [ ! -f "${outpath}/log/3.correct_barcode.log" ]&&[ -f "${outpath}/log/2.skewer.trim.log" ];then
	mkdir -p ${outpath}/3_corrected
	Rscript ./script/3.correct_i7_barcode.R $config
	Rscript ./script/3.Plot_raw_count.spatial_map.R $config
fi
##4
if [ ! -f "${outpath}/log/4.add_barcode_modkit.log" ]&&[ -f "${outpath}/log/3.correct_barcode.log" ];then
        bash ./script/4.add_barcode_modkit.sh $config
fi
##5
if [ -f "${outpath}/log/4.add_barcode_modkit.log" ]&&[ ! -f "${outpath}/log/5.summarize.log" ];then
	mkdir -p ${outpath}/5_Plot
	Rscript ./script/5.summarize_data.R $config
fi
###6
if [ -f "${outpath}/log/5.summarize.log" ]&&[ ! -f "${outpath}/log/6.merge_pixel_or_region.log" ];then
        Rscript ./script/6.merge_pixel_or_region.R $config
fi

