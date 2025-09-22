modkit extract calls 2_calls.all.mapped.HEM23_filtered.sorted.CB_tag.bam extract_calls.tsv -t 20 --mapped-only --cpg --reference hg38.p14.fa --filter-threshold C:0.7 --mod-thresholds m:0.7 --mod-thresholds h:0.6
awk -F'\t' 'NR==1 || $20 != "true"' extract_calls.tsv > extract_calls_filter.tsv

awk '
BEGIN {
    FS=OFS="\t";
    print "read_id", "5mC", "5hmC", "Unmodified";
}
NR > 1 {
    read_id = $1;
    call_code = $14;  # call_code 列：m=甲基化，h=羟甲基化，-=非甲基化

    if (call_code == "m") m[read_id]++;
    else if (call_code == "h") h[read_id]++;
    else if (call_code == "-") u[read_id]++;
}
END {
    # 输出所有 reads 的统计
    for (rid in m) print rid, m[rid], h[rid]+0, u[rid]+0;
    for (rid in h) if (!(rid in m)) print rid, 0, h[rid], u[rid]+0;
    for (rid in u) if (!(rid in m || rid in h)) print rid, 0, 0, u[rid];
}' extract_calls_filter.tsv > read_mod_counts_filter.tsv


awk 'BEGIN {FS=OFS="\t"} NR>1 {print $1, $3}' 3_out_id2coord.sinto.HEM23.tsv | gzip > readid_to_name.tsv.gz

awk '
BEGIN {
    FS = OFS = "\t";
    # 加载 read_id → name 的映射（无表头，直接读取）
    while ((getline < "readid_to_name.tsv.gz") > 0) {
        name[$1] = $2;  # 第一列是 read_id，第二列是 name
    }
    close("readid_to_name.tsv.gz");
}
# 处理 read_mod_counts.tsv（假设它有表头）
NR == 1 {
    print $0, "name";  # 为输出添加 name 列的表头
}
NR > 1 {
    print $0, (name[$1] ? name[$1] : "NA");  # 查找映射，找不到则填 NA
}
' read_mod_counts_filter.tsv  > read_mod_counts_with_name_filter.tsv

