# HAR_MPRA
# /net/shendure/vol9/seq/NEXTSEQ/170314_NS500488_0340_AHHJ52BGX2
# 
# 136K/mm2, 95% passing
# 
# 15/20/15
# 
# Read1: 15bp barcode 
# Index Read: 10bp molecular tag + TG + 1st 8bp of sample id = 20bp 
# Read2: 15bp barcode
# 
# WTC-NPC-1-DNA CGCTACGATT
# WTC-NPC-2-DNA AGGCAGTACC
# WTC-NPC-3-DNA GTTAATTCTA
# Pt2a-NPC-1-DNA CTAATATGAT
# Pt2a-NPC-2-DNA TGGCTGACCA
# Pt2a-NPC-3-DNA ATTGCATAGT
# WTC-NPC-1-RNA AGTCGCGTCG
# WTC-NPC-2-RNA GTTCCGTATC
# WTC-NPC-3-RNA ATCCGGATAA
# Pt2a-NPC-1-RNA GGTTGGCGGT
# Pt2a-NPC-2-RNA CTTCAGGCCA
# Pt2a-NPC-3-RNA GAATGAGGTT
# 
# loaded 3:1 RNA:DNA, but there was very little RNA for either of these, and very little DNA for WTC, so I just did what I could with them. I expect that some will have a lot more reads than others, and not at all even. 
# 
# HAR things look like this:
# AATGATACGGCGACCACCGAGATCTACACCTCGGCATGGACGAGCTGTACAAGTAGGAATTCNNNNNNNNNNNNNNNCATTGCGTGAACCGACAATTCGTCGAGGGACCTAATAACTTCGNNNNNNNNNNTG##########ATCTCGTATGCCGTCTTCTGCTTG

##################
## RUN INFO
##################

/net/shendure/vol9/seq/NEXTSEQ/170314_NS500488_0340_AHHJ52BGX2

less /net/shendure/vol9/seq/NEXTSEQ/170314_NS500488_0340_AHHJ52BGX2/RunInfo.xml
#       <Read Number="1" NumCycles="15" IsIndexedRead="N" />
#       <Read Number="2" NumCycles="20" IsIndexedRead="Y" />
#       <Read Number="3" NumCycles="15" IsIndexedRead="N" />

##################
## PROCESSING
##################

mkdir /net/shendure/vol10/nobackup/tmp/170314_NS500488_0340_AHHJ52BGX2
cd /net/shendure/vol10/nobackup/tmp/170314_NS500488_0340_AHHJ52BGX2

module load gmp/5.0.2 mpfr/3.1.0 mpc/0.8.2 gcc/4.9.1 bcl2fastq/latest

bcl2fastq --create-fastq-for-index-reads --sample-sheet NoSampleSheet.txt --use-bases-mask 'Y*,I*,Y*' -R /net/shendure/vol9/seq/NEXTSEQ/170314_NS500488_0340_AHHJ52BGX2/ -o ./ --minimum-trimmed-read-length 0 --mask-short-adapter-reads 0

zcat Undetermined_S0_L00*_I1_001.fastq.gz | wc -l | awk '{ print $1/4 }'
# 336669180

zcat Undetermined_S0_L00*_R1_001.fastq.gz | awk '{ count+=1; if (count == 2) { print }; if (count == 4) { count=0 } }' |  head -n 100000 | sort | uniq -c | sort -nr | head -n 10 
#    1481 GTCGAGGGACCTAAT
#      21 TAGAAGAATCGGACT
#      20 CCAAGACCTTCTGAA
#      20 CAGCTACCTTAGTAC
#      20 CAATTACCGGTCTCA
#      19 CTGGTTCTGCGCTCC
#      19 CCAAGCGTCCTGGAT
#      19 CAGAGGCTCCTTAGA
#      18 TTATCTGCCTTGGAT
#      18 CGGTTCATATGATTG
     
zcat Undetermined_S0_L00*_I1_001.fastq.gz | awk '{ count+=1; if (count == 2) { print }; if (count == 4) { count=0 } }' | head -n 100000 | sort | uniq -c | sort -nr | head -n 10 
#    2412 GGGGGGGGGGGGGGGGGGGG
#      86 CCCCCCCCCCTGGGTTGGCG
#      60 TATAGCATACATTATACGAA
#      27 AGGGGGGGGGGGGGGGGGGG
#      26 GGGGGGGGGGGGGGGGGGGT
#      17 CCCCCCCCCTGGGTTGGCGG
#      16 CGGGGGGGGGGGGGGGGGGG
#      14 TGGGGGGGGGGGGGGGGGGG
#      14 CCCCCCCCCGTGGGTTGGCG
#      13 GGGGGGTGGGGGGGGGGGGG
     
zcat Undetermined_S0_L00*_I1_001.fastq.gz | awk '{ count+=1; if (count == 2) { print }; if (count == 4) { count=0 } }' | head -n 100000 | cut -c 11-12 | sort | uniq -c | sort -nr | head -n 10 
#   93036 TG
#    4739 GG
#     449 CG
#     409 GC
#     334 GA
#     208 AG
#     171 TA
#     140 AT
#     127 TT
#     124 GT

zcat Undetermined_S0_L00*_I1_001.fastq.gz | awk '{ count+=1; if (count == 2) { print }; if (count == 4) { count=0 } }' | head -n 100000 | cut -c 13- | sort | uniq -c | sort -nr | head -n 20 
#   16761 CTTCAGGC
#   13832 AGTCGCGT
#   13750 GGTTGGCG
#    9462 GAATGAGG
#    8683 GTTCCGTA
#    4094 CGCTACGA
#    3942 AGGCAGTA
#    3762 TGGCTGAC
#    3476 ATCCGGAT
#    3317 CTAATATG
#    3153 ATTGCATA
#    2823 GGGGGGGG
#    2480 GTTAATTC
#     399 AATGAGGT
#     368 GTTGGCGG
#     341 TTCAGGCC
#     266 CAGGAATC
#     194 GTTCCGTT
#     194 GTCGCGTC
#     148 TCCGGATA

zcat Undetermined_S0_L00*_R2_001.fastq.gz | awk '{ count+=1; if (count == 2) { print }; if (count == 4) { count=0 } }' | head -n 100000 | sort | uniq -c | sort -nr | head -n 10 
#    2070 GGGGGGGGGGGGGGG
#      21 AGTCCGATTCTTCTA
#      20 TTCAGAAGGTCTTGG
#      20 TGAGACCGGTAATTG
#      20 GTACTAAGGTAGCTG
#      19 TCTAAGGAGCCTCTG
#      19 GGAGCGCAGAACCAG
#      19 ATCCAGGACGCTTGG
#      18 CAATCATATGAACCG
#      18 ATCCAAGGCAGATAA

##########################
## Basic read processing #
##########################

mkdir -p /net/shendure/vol1/home/mkircher/regulatory_tests/HARs/counts/NPC
cd /net/shendure/vol1/home/mkircher/regulatory_tests/HARs/counts/NPC

# CGCTACGATT	WTC-NPC-1-DNA
# AGGCAGTACC	WTC-NPC-2-DNA
# GTTAATTCTA	WTC-NPC-3-DNA
# CTAATATGAT	Pt2a-NPC-1-DNA
# TGGCTGACCA	Pt2a-NPC-2-DNA
# ATTGCATAGT	Pt2a-NPC-3-DNA
# AGTCGCGTCG	WTC-NPC-1-RNA
# GTTCCGTATC	WTC-NPC-2-RNA
# ATCCGGATAA	WTC-NPC-3-RNA
# GGTTGGCGGT	Pt2a-NPC-1-RNA
# CTTCAGGCCA	Pt2a-NPC-2-RNA
# GAATGAGGTT	Pt2a-NPC-3-RNA

vim index.lst
#--------------------
#barcode	sample
NNNNNNNNNNTGCGCTACGA	WTC-NPC-1-DNA
NNNNNNNNNNTGAGGCAGTA	WTC-NPC-2-DNA
NNNNNNNNNNTGGTTAATTC	WTC-NPC-3-DNA
NNNNNNNNNNTGCTAATATG	Pt2a-NPC-1-DNA
NNNNNNNNNNTGTGGCTGAC	Pt2a-NPC-2-DNA
NNNNNNNNNNTGATTGCATA	Pt2a-NPC-3-DNA
NNNNNNNNNNTGAGTCGCGT	WTC-NPC-1-RNA
NNNNNNNNNNTGGTTCCGTA	WTC-NPC-2-RNA
NNNNNNNNNNTGATCCGGAT	WTC-NPC-3-RNA
NNNNNNNNNNTGGGTTGGCG	Pt2a-NPC-1-RNA
NNNNNNNNNNTGCTTCAGGC	Pt2a-NPC-2-RNA
NNNNNNNNNNTGGAATGAGG	Pt2a-NPC-3-RNA
#--------------------

vim processLane.sh
#----------------------------
#!/bin/bash
( paste <( zcat /net/shendure/vol10/nobackup/tmp/170314_NS500488_0340_AHHJ52BGX2/Undetermined_S0_L00${1}_R1_001.fastq.gz ) <( zcat /net/shendure/vol10/nobackup/tmp/170314_NS500488_0340_AHHJ52BGX2/Undetermined_S0_L00${1}_I1_001.fastq.gz ) <( zcat /net/shendure/vol10/nobackup/tmp/170314_NS500488_0340_AHHJ52BGX2/Undetermined_S0_L00${1}_R2_001.fastq.gz ) | awk '{ count+=1; if ((count == 1) || (count == 3)) { print $1 } else { print $1$2$3 }; if (count == 4) { count=0 } }' | ~mkircher/bin/pipeline2.0/SplitFastQdoubleIndexBAM.py -p -l 20 -m 0 -s 36 -i index.lst --setN=NNNNNNNNNNyyyyyyyyyy --remove --summary | ~mkircher/bin/pipeline2.0/MergeTrimReadsBAM.py -p --mergeoverlap -f CATTGCGTGAACCGACAATTCGTCGAGGGACCTAATAACTTCG -c GTCGAGGGACCTAATAACTTCGCCCCCCCCCCTG > sample_l${1}.bam ) 2> processing_stats_l${1}.log
#----------------------------
chmod +x processLane.sh

for i in {1..4}; do qsub -N processCounts -l virtual_free=3G,mfree=3G -m n ./processLane.sh ${i}; done

#Multiple lanes / NextSeq:
cat *.log  | grep clusters | awk '{ print $5,$2 }' | sort | ~/bin/src/sum_numbers_category.py | column -t 
# Pt2a-NPC-1-DNA  11010851
# Pt2a-NPC-1-RNA  45259038
# Pt2a-NPC-2-DNA  12291000
# Pt2a-NPC-2-RNA  54928840
# Pt2a-NPC-3-DNA  10927342
# Pt2a-NPC-3-RNA  31235753
# WTC-NPC-1-DNA   13446094
# WTC-NPC-1-RNA   45325289
# WTC-NPC-2-DNA   12787683
# WTC-NPC-2-RNA   28429008
# WTC-NPC-3-DNA   8006761
# WTC-NPC-3-RNA   11595699
# unknown         51425822
# 
# -> 15% unknown/barcodes not included
# 
# Pt2a-NPC-1-DNA	4%
# Pt2a-NPC-1-RNA	16%
# Pt2a-NPC-2-DNA	4%
# Pt2a-NPC-2-RNA	19%
# Pt2a-NPC-3-DNA	4%
# Pt2a-NPC-3-RNA	11%
# WTC-NPC-1-DNA	5%
# WTC-NPC-1-RNA	16%
# WTC-NPC-2-DNA	4%
# WTC-NPC-2-RNA	10%
# WTC-NPC-3-DNA	3%
# WTC-NPC-3-RNA	4%

cat *.log  | grep Merged | sed s/';'//g | awk '{ total += $2; merged += $5+$8; paired += $11 }END{ printf("Total: %d; Merged %d (%.2f%%); Paired %d (%.2f%%)\n",total,merged,merged/total*100,paired,paired/total*100) }'
# Total: 285243358; Merged 276776814 (97.03%); Paired 8425521 (2.95%)

###################
## Extract counts #
###################

vim extractRawCounts.sh
#---------------------------
#!/bin/bash
mkdir -p raw_counts
echo ${1}
samtools merge -u sample_l{1..4}.bam | samtools view -M 20 -F pq -r ${1} - | awk 'BEGIN{ OFS= "\t" }{ for (i=12; i<=NF; i++) { if ($i ~ /^XI:Z:/) print $10,substr($i,6,10) }}' | sort | uniq -c | awk 'BEGIN{ OFS="\t" }{ print $2,$3,$1 }' | gzip -c > raw_counts/${1}.tsv.gz
#---------------------------
chmod +x extractRawCounts.sh

for i in $(tail -n +2 index.lst | cut -f 2); do echo ./extractRawCounts.sh $i; done > extractBarcodes_cmds.lst
~/bin/commands2arrayjob.sh extractBarcodes_cmds.lst -N extrBarcodes -l virtual_free=2G,mfree=2G -m n 

for i in raw_counts/*.tsv.gz; do echo $(basename $i) $(zcat $i | cut -f 3 | sort | uniq -c | sort -nr | head -n 5); done 
# Pt2a-NPC-1-DNA.tsv.gz 10239462 1 199580 2 2434 3 26 4
# Pt2a-NPC-1-RNA.tsv.gz 34173817 1 3958175 2 502038 3 67096 4 11560 5
# Pt2a-NPC-2-DNA.tsv.gz 11566088 1 199369 2 1914 3 14 4
# Pt2a-NPC-2-RNA.tsv.gz 43371509 1 4168771 2 437881 3 46462 4 5024 5
# Pt2a-NPC-3-DNA.tsv.gz 10130443 1 207294 2 2736 3 31 4 1 5
# Pt2a-NPC-3-RNA.tsv.gz 21428980 1 3369522 2 567855 3 89516 4 13644 5
# WTC-NPC-1-DNA.tsv.gz 11957504 1 473794 2 39998 3 5498 4 740 5
# WTC-NPC-1-RNA.tsv.gz 36352675 1 3292768 2 314093 3 29892 4 3006 5
# WTC-NPC-2-DNA.tsv.gz 11856637 1 269695 2 5516 3 206 4 32 5
# WTC-NPC-2-RNA.tsv.gz 22010913 1 2309355 2 267496 3 32185 4 4295 5
# WTC-NPC-3-DNA.tsv.gz 7338813 1 179153 2 3501 3 113 4 5 5
# WTC-NPC-3-RNA.tsv.gz 9114960 1 905423 2 97276 3 10433 4 1073 5

for i in raw_counts/*.tsv.gz; do echo $(basename $i) $(zcat $i | awk '{ lines+=1; counts+=$NF }END{ print lines,counts }'); done | column -t
# Pt2a-NPC-1-DNA.tsv.gz  10441502  10646028
# Pt2a-NPC-1-RNA.tsv.gz  38722084  44001651
# Pt2a-NPC-2-DNA.tsv.gz  11767385  11970624
# Pt2a-NPC-2-RNA.tsv.gz  48030491  53239009
# Pt2a-NPC-3-DNA.tsv.gz  10340505  10553368
# Pt2a-NPC-3-RNA.tsv.gz  25471970  30313241
# WTC-NPC-1-DNA.tsv.gz   12477705  13051926
# WTC-NPC-1-RNA.tsv.gz   39993125  44019904
# WTC-NPC-2-DNA.tsv.gz   12132097  12413628
# WTC-NPC-2-RNA.tsv.gz   24625834  27593834
# WTC-NPC-3-DNA.tsv.gz   7521585   7708099
# WTC-NPC-3-RNA.tsv.gz   10129290  11265502

vim filterCounts.sh
#----------------------------
#!/bin/bash
mkdir -p filtered_counts
zcat ${1} | awk 'BEGIN{ OFS="\t" }{ if (length($1) == 15) { print } }' | gzip -c > filtered_counts/$(basename ${1})
#----------------------------
chmod +x filterCounts.sh

for i in $(tail -n +2 index.lst  | cut -f 2); do echo ./filterCounts.sh raw_counts/${i}.tsv.gz; done > filter_cmds.lst
~mkircher/bin/commands2arrayjob.sh filter_cmds.lst -N filtCounts -l virtual_free=2G,mfree=2G -m n -hold_jid 331709

for i in filtered_counts/*.tsv.gz; do echo $(basename $i) $(zcat $i | awk 'BEGIN{ last="" }{ if($1 != last) { barcodes+=1; last=$1 }; lines+=1; counts+=$NF }END{ print barcodes,lines,counts }'); done | column -t
# Pt2a-NPC-1-DNA.tsv.gz  736260   10023824  10220346
# Pt2a-NPC-1-RNA.tsv.gz  1709009  37143585  42223068
# Pt2a-NPC-2-DNA.tsv.gz  792724   11296181  11491480
# Pt2a-NPC-2-RNA.tsv.gz  2003190  46065554  51071218
# Pt2a-NPC-3-DNA.tsv.gz  743638   9925926   10130358
# Pt2a-NPC-3-RNA.tsv.gz  1390839  24423757  29076050
# WTC-NPC-1-DNA.tsv.gz   856676   11975272  12527166
# WTC-NPC-1-RNA.tsv.gz   1727156  38372449  42242766
# WTC-NPC-2-DNA.tsv.gz   826631   11646371  11917482
# WTC-NPC-2-RNA.tsv.gz   1311881  23623676  26476891
# WTC-NPC-3-DNA.tsv.gz   618591   7220036   7399470
# WTC-NPC-3-RNA.tsv.gz   743239   9713890   10805480

# Frequent UMIs?
for i in filtered_counts/*.tsv.gz; do echo $(basename $i); zcat $i | cut -f 2 | sort | uniq -c | sort -nr | head; echo; done > freqUMIs.txt

for i in filtered_counts/*.tsv.gz; do echo $(basename $i | sed s/.tsv.gz//g) $(zcat $i | awk '{ total+=$3; if ($2 ~ /(AAAAAAAAAA|.AAAAAAAAA|A.AAAAAAAA|AA.AAAAAAA|AAA.AAAAAA|AAAA.AAAAA|AAAAA.AAAA|AAAAAA.AAA|AAAAAAA.AA|AAAAAAAA.A|AAAAAAAAA.)/) count+=$3 }END{ printf("%d/%d (%.2f%%)",count,total,count/total*100) }'); done
# Pt2a-NPC-1-DNA 2779/10220346 (0.03%)
# Pt2a-NPC-1-RNA 36618/42223068 (0.09%)
# Pt2a-NPC-2-DNA 11931/11491480 (0.10%)
# Pt2a-NPC-2-RNA 47290/51071218 (0.09%)
# Pt2a-NPC-3-DNA 251/10130358 (0.00%)
# Pt2a-NPC-3-RNA 29232/29076050 (0.10%)
# WTC-NPC-1-DNA 10777/12527166 (0.09%)
# WTC-NPC-1-RNA 862/42242766 (0.00%)
# WTC-NPC-2-DNA 18611/11917482 (0.16%)
# WTC-NPC-2-RNA 513/26476891 (0.00%)
# WTC-NPC-3-DNA 4839/7399470 (0.07%)
# WTC-NPC-3-RNA 398/10805480 (0.00%)

for i in filtered_counts/*.tsv.gz; do echo $(basename $i | sed s/.tsv.gz//g) $(zcat $i | awk '{ total+=$3; if ($2 ~ /(CCCCCCCCCC|.CCCCCCCCC|C.CCCCCCCC|CC.CCCCCCC|CCC.CCCCCC|CCCC.CCCCC|CCCCC.CCCC|CCCCCC.CCC|CCCCCCC.CC|CCCCCCCC.C|CCCCCCCCC.)/) count+=$3 }END{ printf("%d/%d (%.2f%%)",count,total,count/total*100) }'); done
# Pt2a-NPC-1-DNA 28459/10220346 (0.28%)
# Pt2a-NPC-1-RNA 810025/42223068 (1.92%)
# Pt2a-NPC-2-DNA 11127/11491480 (0.10%)
# Pt2a-NPC-2-RNA 232528/51071218 (0.46%)
# Pt2a-NPC-3-DNA 53945/10130358 (0.53%)
# Pt2a-NPC-3-RNA 101732/29076050 (0.35%)
# WTC-NPC-1-DNA 141578/12527166 (1.13%)
# WTC-NPC-1-RNA 81821/42242766 (0.19%)
# WTC-NPC-2-DNA 288723/11917482 (2.42%)
# WTC-NPC-2-RNA 35295/26476891 (0.13%)
# WTC-NPC-3-DNA 186182/7399470 (2.52%)
# WTC-NPC-3-RNA 10313/10805480 (0.10%)

for i in filtered_counts/*.tsv.gz; do echo $(basename $i | sed s/.tsv.gz//g) $(zcat $i | awk '{ total+=$3; if ($2 ~ /(GGGGGGGGGG|.GGGGGGGGG|G.GGGGGGGG|GG.GGGGGGG|GGG.GGGGGG|GGGG.GGGGG|GGGGG.GGGG|GGGGGG.GGG|GGGGGGG.GG|GGGGGGGG.G|GGGGGGGGG.)/) count+=$3 }END{ printf("%d/%d (%.2f%%)",count,total,count/total*100) }'); done
# Pt2a-NPC-1-DNA 116/10220346 (0.00%)
# Pt2a-NPC-1-RNA 4106/42223068 (0.01%)
# Pt2a-NPC-2-DNA 220/11491480 (0.00%)
# Pt2a-NPC-2-RNA 1021/51071218 (0.00%)
# Pt2a-NPC-3-DNA 805/10130358 (0.01%)
# Pt2a-NPC-3-RNA 285/29076050 (0.00%)
# WTC-NPC-1-DNA 280/12527166 (0.00%)
# WTC-NPC-1-RNA 2843/42242766 (0.01%)
# WTC-NPC-2-DNA 230/11917482 (0.00%)
# WTC-NPC-2-RNA 1153/26476891 (0.00%)
# WTC-NPC-3-DNA 123/7399470 (0.00%)
# WTC-NPC-3-RNA 410/10805480 (0.00%)

for i in filtered_counts/*.tsv.gz; do echo $(basename $i | sed s/.tsv.gz//g) $(zcat $i | awk '{ total+=$3; if ($2 ~ /(TTTTTTTTTT|.TTTTTTTTT|T.TTTTTTTT|TT.TTTTTTT|TTT.TTTTTT|TTTT.TTTTT|TTTTT.TTTT|TTTTTT.TTT|TTTTTTT.TT|TTTTTTTT.T|TTTTTTTTT.)/) count+=$3 }END{ printf("%d/%d (%.2f%%)",count,total,count/total*100) }'); done
# Pt2a-NPC-1-DNA 59/10220346 (0.00%)
# Pt2a-NPC-1-RNA 3133/42223068 (0.01%)
# Pt2a-NPC-2-DNA 82/11491480 (0.00%)
# Pt2a-NPC-2-RNA 3174/51071218 (0.01%)
# Pt2a-NPC-3-DNA 745/10130358 (0.01%)
# Pt2a-NPC-3-RNA 1627/29076050 (0.01%)
# WTC-NPC-1-DNA 1308/12527166 (0.01%)
# WTC-NPC-1-RNA 1711/42242766 (0.00%)
# WTC-NPC-2-DNA 1273/11917482 (0.01%)
# WTC-NPC-2-RNA 909/26476891 (0.00%)
# WTC-NPC-3-DNA 611/7399470 (0.01%)
# WTC-NPC-3-RNA 394/10805480 (0.00%)

#########################################
## Match to HARperm and HAR assignments
#########################################

mkdir assignment_counts_HARperm_UMIs

for i in filtered_counts/*.tsv.gz; do echo $(basename $i); join -1 1 -2 2 -t"$(echo -e '\t')" <( zcat $i | cut -f 1 | sort | uniq -c | awk 'BEGIN{ OFS="\t" }{ print $2,$1 }' ) <( zcat ~/regulatory_tests/HARs/design/probes_7.22.2015_HAR_permutation_library_1049_seqs.inclBarcode.tsv.gz | grep -v '^#' | cut -f -2 | sort -k2,2 ) | gzip -c > assignment_counts_HARperm_UMIs/$(basename $i); done

for i in assignment_counts_HARperm_UMIs/*.tsv.gz; do echo $(basename $i) $(zcat $i | awk '{ lines+=1; counts+=$2 }END{ print lines,counts }'); done | column -t
# Pt2a-NPC-1-DNA.tsv.gz  263  1180
# Pt2a-NPC-1-RNA.tsv.gz  621  4609
# Pt2a-NPC-2-DNA.tsv.gz  277  1358
# Pt2a-NPC-2-RNA.tsv.gz  750  5780
# Pt2a-NPC-3-DNA.tsv.gz  248  1276
# Pt2a-NPC-3-RNA.tsv.gz  450  3015
# WTC-NPC-1-DNA.tsv.gz   600  1811
# WTC-NPC-1-RNA.tsv.gz   594  4777
# WTC-NPC-2-DNA.tsv.gz   320  1474
# WTC-NPC-2-RNA.tsv.gz   433  2951
# WTC-NPC-3-DNA.tsv.gz   206  953
# WTC-NPC-3-RNA.tsv.gz   209  1153

for i in assignment_counts_HARperm_UMIs/*.tsv.gz; do echo $(basename $i) $(zcat $i | awk '{ if ($2 > 2) { lines+=1; counts+=$2 } }END{ print lines,counts }'); done | column -t
# Pt2a-NPC-1-DNA.tsv.gz  50  960
# Pt2a-NPC-1-RNA.tsv.gz  65  4028
# Pt2a-NPC-2-DNA.tsv.gz  51  1126
# Pt2a-NPC-2-RNA.tsv.gz  70  5074
# Pt2a-NPC-3-DNA.tsv.gz  52  1074
# Pt2a-NPC-3-RNA.tsv.gz  62  2611
# WTC-NPC-1-DNA.tsv.gz   53  1252
# WTC-NPC-1-RNA.tsv.gz   81  4231
# WTC-NPC-2-DNA.tsv.gz   53  1199
# WTC-NPC-2-RNA.tsv.gz   59  2557
# WTC-NPC-3-DNA.tsv.gz   48  790
# WTC-NPC-3-RNA.tsv.gz   50  986

mkdir assignment_counts_HAR_UMIs

for i in filtered_counts/*.tsv.gz; do echo $(basename $i); join -1 1 -2 2 -t"$(echo -e '\t')" <( zcat $i | cut -f 1 | sort | uniq -c | awk 'BEGIN{ OFS="\t" }{ print $2,$1 }' ) <( zcat /net/shendure/vol1/home/mkircher/regulatory_tests/HARs/design/design_9.2.2015_HAR_library_2440_seqs.inclBarcode.tsv.gz | grep -v '^#' | cut -f -2 | sort -k2,2 ) | gzip -c > assignment_counts_HAR_UMIs/$(basename $i); done

for i in assignment_counts_HAR_UMIs/*.tsv.gz; do echo $(basename $i) $(zcat $i | awk '{ lines+=1; counts+=$2 }END{ print lines,counts }'); done | column -t
# Pt2a-NPC-1-DNA.tsv.gz  198988  9184625
# Pt2a-NPC-1-RNA.tsv.gz  201696  33880522
# Pt2a-NPC-2-DNA.tsv.gz  199537  10350837
# Pt2a-NPC-2-RNA.tsv.gz  201946  41867492
# Pt2a-NPC-3-DNA.tsv.gz  198822  9078893
# Pt2a-NPC-3-RNA.tsv.gz  201013  22134610
# WTC-NPC-1-DNA.tsv.gz   199619  10925070
# WTC-NPC-1-RNA.tsv.gz   201576  35019545
# WTC-NPC-2-DNA.tsv.gz   199646  10646637
# WTC-NPC-2-RNA.tsv.gz   200785  21511141
# WTC-NPC-3-DNA.tsv.gz   197041  6595730
# WTC-NPC-3-RNA.tsv.gz   197591  8856236

for i in assignment_counts_HAR_UMIs/*.tsv.gz; do echo $(basename $i) $(zcat $i | awk '{ if ($2 > 2) { lines+=1; counts+=$2 } }END{ print lines,counts }'); done | column -t
# Pt2a-NPC-1-DNA.tsv.gz  190097  9171093
# Pt2a-NPC-1-RNA.tsv.gz  199185  33876908
# Pt2a-NPC-2-DNA.tsv.gz  192124  10339541
# Pt2a-NPC-2-RNA.tsv.gz  199817  41864452
# Pt2a-NPC-3-DNA.tsv.gz  189816  9065130
# Pt2a-NPC-3-RNA.tsv.gz  197160  22128902
# WTC-NPC-1-DNA.tsv.gz   192234  10913813
# WTC-NPC-1-RNA.tsv.gz   199047  35015822
# WTC-NPC-2-DNA.tsv.gz   192124  10635178
# WTC-NPC-2-RNA.tsv.gz   196659  21504918
# WTC-NPC-3-DNA.tsv.gz   183934  6575868
# WTC-NPC-3-RNA.tsv.gz   186447  8839376

# Pt2a-NPC-1-DNA	99.87%	99.99%
# Pt2a-NPC-1-RNA	99.69%	99.99%
# Pt2a-NPC-2-DNA	99.86%	99.99%
# Pt2a-NPC-2-RNA	99.63%	99.99%
# Pt2a-NPC-3-DNA	99.88%	99.99%
# Pt2a-NPC-3-RNA	99.78%	99.99%
# WTC-NPC-1-DNA	99.70%	99.98%
# WTC-NPC-1-RNA	99.71%	99.99%
# WTC-NPC-2-DNA	99.84%	99.99%
# WTC-NPC-2-RNA	99.78%	99.99%
# WTC-NPC-3-DNA	99.90%	99.99%
# WTC-NPC-3-RNA	99.89%	99.99%
# Average      	99.79%	99.99%
# 
# ==> 99.99% pure HAR libraries

######################
## Summary           #
######################

# .               Barcodes  UMIcounts  Counts    Unique  min1    .         Barcodes  Used   min2    .         Barcodes  Used
# Pt2a-NPC-1-DNA  736260    10023824   10220346  98.1%   198988  9184625   81.6%     89.9%  190097  9171093   77.9%     89.7%
# Pt2a-NPC-1-RNA  1709009   37143585   42223068  88.0%   201696  33880522  82.7%     80.2%  199185  33876908  81.6%     80.2%
# Pt2a-NPC-2-DNA  792724    11296181   11491480  98.3%   199537  10350837  81.8%     90.1%  192124  10339541  78.7%     90.0%
# Pt2a-NPC-2-RNA  2003190   46065554   51071218  90.2%   201946  41867492  82.8%     82.0%  199817  41864452  81.9%     82.0%
# Pt2a-NPC-3-DNA  743638    9925926    10130358  98.0%   198822  9078893   81.5%     89.6%  189816  9065130   77.8%     89.5%
# Pt2a-NPC-3-RNA  1390839   24423757   29076050  84.0%   201013  22134610  82.4%     76.1%  197160  22128902  80.8%     76.1%
# WTC-NPC-1-DNA   856676    11975272   12527166  95.6%   199619  10925070  81.8%     87.2%  192234  10913813  78.8%     87.1%
# WTC-NPC-1-RNA   1727156   38372449   42242766  90.8%   201576  35019545  82.6%     82.9%  199047  35015822  81.6%     82.9%
# WTC-NPC-2-DNA   826631    11646371   11917482  97.7%   199646  10646637  81.8%     89.3%  192124  10635178  78.7%     89.2%
# WTC-NPC-2-RNA   1311881   23623676   26476891  89.2%   200785  21511141  82.3%     81.2%  196659  21504918  80.6%     81.2%
# WTC-NPC-3-DNA   618591    7220036    7399470   97.6%   197041  6595730   80.8%     89.1%  183934  6575868   75.4%     88.9%
# WTC-NPC-3-RNA   743239    9713890    10805480  89.9%   197591  8856236   81.0%     82.0%  186447  8839376   76.4%     81.8%

######################
## Correlation check #
######################

Rscript barcodeCounts_correlation.R assignment_counts_HAR_UMIs
# [1] "WTC-NPC"
# [1] "DNA"
# [1] 0.9690922
# [1] 0.9594674
# [1] 0.959998
# [1] "RNA"
# [1] 0.9792575
# [1] 0.9638171
# [1] 0.9649053
# [1] "-----"
# [1] "Pt2a-NPC"
# [1] "DNA"
# [1] 0.9661346
# [1] 0.9658731
# [1] 0.9642839
# [1] "RNA"
# [1] 0.9887701
# [1] 0.986711
# [1] 0.9856516
# [1] "-----"


R --vanilla --quiet
#-------------------------------------------

for (folder in c('assignment_counts_HAR_UMIs')) 
{
  for (sample in c('WTC-NPC','Pt2a-NPC')) {
  RNA1 <- read.table(sprintf("%s/%s-1-RNA.tsv.gz",folder,sample),as.is=T)
  RNA2 <- read.table(sprintf("%s/%s-2-RNA.tsv.gz",folder,sample),as.is=T)
  RNA3 <- read.table(sprintf("%s/%s-3-RNA.tsv.gz",folder,sample),as.is=T)
  DNA1 <- read.table(sprintf("%s/%s-1-DNA.tsv.gz",folder,sample),as.is=T)
  DNA2 <- read.table(sprintf("%s/%s-2-DNA.tsv.gz",folder,sample),as.is=T)
  DNA3 <- read.table(sprintf("%s/%s-3-DNA.tsv.gz",folder,sample),as.is=T)

  data1 <- merge(RNA1,DNA1,by="V1")[,-3]
  colnames(data1) <- c("barcodes","X","Y","name")
  data1$RNA <- data1$X/sum(data1$X)*10^6
  data1$DNA <- data1$Y/sum(data1$Y)*10^6
  data1$ratio <- data1$RNA/data1$DNA

  data2 <- merge(RNA2,DNA2,by="V1")[,-3]
  colnames(data2) <- c("barcodes","X","Y","name")
  data2$RNA <- data2$X/sum(data2$X)*10^6
  data2$DNA <- data2$Y/sum(data2$Y)*10^6
  data2$ratio <- data2$RNA/data2$DNA

  data3 <- merge(RNA3,DNA3,by="V1")[,-3]
  colnames(data3) <- c("barcodes","X","Y","name")
  data3$RNA <- data3$X/sum(data3$X)*10^6
  data3$DNA <- data3$Y/sum(data3$Y)*10^6
  data3$ratio <- data3$RNA/data3$DNA

  mean(c(cor(merge(data1,data2,by="name")[,c('ratio.x','ratio.y')],method="spearman")[2],cor(merge(data1,data3,by="name")[,c('ratio.x','ratio.y')],method="spearman")[2],cor(merge(data2,data3,by="name")[,c('ratio.x','ratio.y')],method="spearman")[2]))
  
  png(sprintf("%s/%s-allTags.png",folder,sample),width=800,height=800,type="cairo")
  par(mfrow=c(3,3))

  selSamples <- c(sprintf('%s-1',sample),sprintf('%s-2',sample))
  res1 <- merge(data1,data2,by="name")
  plot(res1$DNA.x,res1$DNA.y, main=sprintf("DNA %s vs %s\n(p-rho %.2f, s-rho %.2f)",selSamples[1],selSamples[2], cor(res1$DNA.x,res1$DNA.y,method="pearson"), cor(res1$DNA.x,res1$DNA.y,method="spearman")), xlab=sprintf("%s normalized DNA count per insert",selSamples[1]), ylab=sprintf("%s normalized DNA count per insert",selSamples[2]), pch=19, cex=0.3, col = "darkblue", las=1, cex.axis=1.3, cex.lab=1.3, cex.main=1.2)
  abline(0,1,lty=2)
  plot(res1$RNA.x,res1$RNA.y, main=sprintf("RNA %s vs %s\n(p-rho %.2f, s-rho %.2f)",selSamples[1],selSamples[2], cor(res1$RNA.x,res1$RNA.y,method="pearson"), cor(res1$RNA.x,res1$RNA.y,method="spearman")), xlab=sprintf("%s normalized RNA count per insert",selSamples[1]), ylab=sprintf("%s normalized RNA count per insert",selSamples[2]), pch=19, cex=0.3, col = "darkgreen", las=1, cex.axis=1.3, cex.lab=1.3, cex.main=1.2)
  abline(0,1,lty=2)
  plot(log2(res1$ratio.x),log2(res1$ratio.y), main=sprintf("RNA/DNA %s vs %s\n(p-rho %.2f, s-rho %.2f)",selSamples[1],selSamples[2], cor(log2(res1$ratio.x),log2(res1$ratio.y),method="pearson"), cor(res1$ratio.x,res1$ratio.y,method="spearman")), xlab=sprintf("%s log2 RNA/DNA per insert",selSamples[1]), ylab=sprintf("%s log2 RNA/DNA per insert",selSamples[2]), pch=19, cex=0.3, col = "darkred", las=1, cex.axis=1.3, cex.lab=1.3, cex.main=1.2)
  abline(0,1,lty=2)

  selSamples <- c(sprintf('%s-1',sample),sprintf('%s-3',sample))
  res1 <- merge(data1,data3,by="name")
  plot(res1$DNA.x,res1$DNA.y, main=sprintf("DNA %s vs %s\n(p-rho %.2f, s-rho %.2f)",selSamples[1],selSamples[2], cor(res1$DNA.x,res1$DNA.y,method="pearson"), cor(res1$DNA.x,res1$DNA.y,method="spearman")), xlab=sprintf("%s normalized DNA count per insert",selSamples[1]), ylab=sprintf("%s normalized DNA count per insert",selSamples[2]), pch=19, cex=0.3, col = "darkblue", las=1, cex.axis=1.3, cex.lab=1.3, cex.main=1.2)
  abline(0,1,lty=2)
  plot(res1$RNA.x,res1$RNA.y, main=sprintf("RNA %s vs %s\n(p-rho %.2f, s-rho %.2f)",selSamples[1],selSamples[2], cor(res1$RNA.x,res1$RNA.y,method="pearson"), cor(res1$RNA.x,res1$RNA.y,method="spearman")), xlab=sprintf("%s normalized RNA count per insert",selSamples[1]), ylab=sprintf("%s normalized RNA count per insert",selSamples[2]), pch=19, cex=0.3, col = "darkgreen", las=1, cex.axis=1.3, cex.lab=1.3, cex.main=1.2)
  abline(0,1,lty=2)
  plot(log2(res1$ratio.x),log2(res1$ratio.y), main=sprintf("RNA/DNA %s vs %s\n(p-rho %.2f, s-rho %.2f)",selSamples[1],selSamples[2], cor(log2(res1$ratio.x),log2(res1$ratio.y),method="pearson"), cor(res1$ratio.x,res1$ratio.y,method="spearman")), xlab=sprintf("%s log2 RNA/DNA per insert",selSamples[1]), ylab=sprintf("%s log2 RNA/DNA per insert",selSamples[2]), pch=19, cex=0.3, col = "darkred", las=1, cex.axis=1.3, cex.lab=1.3, cex.main=1.2)
  abline(0,1,lty=2)

  selSamples <- c(sprintf('%s-2',sample),sprintf('%s-3',sample))
  res1 <- merge(data2,data3,by="name")
  plot(res1$DNA.x,res1$DNA.y, main=sprintf("DNA %s vs %s\n(p-rho %.2f, s-rho %.2f)",selSamples[1],selSamples[2], cor(res1$DNA.x,res1$DNA.y,method="pearson"), cor(res1$DNA.x,res1$DNA.y,method="spearman")), xlab=sprintf("%s normalized DNA count per insert",selSamples[1]), ylab=sprintf("%s normalized DNA count per insert",selSamples[2]), pch=19, cex=0.3, col = "darkblue", las=1, cex.axis=1.3, cex.lab=1.3, cex.main=1.2)
  abline(0,1,lty=2)
  plot(res1$RNA.x,res1$RNA.y, main=sprintf("RNA %s vs %s\n(p-rho %.2f, s-rho %.2f)",selSamples[1],selSamples[2], cor(res1$RNA.x,res1$RNA.y,method="pearson"), cor(res1$RNA.x,res1$RNA.y,method="spearman")), xlab=sprintf("%s normalized RNA count per insert",selSamples[1]), ylab=sprintf("%s normalized RNA count per insert",selSamples[2]), pch=19, cex=0.3, col = "darkgreen", las=1, cex.axis=1.3, cex.lab=1.3, cex.main=1.2)
  abline(0,1,lty=2)
  plot(log2(res1$ratio.x),log2(res1$ratio.y), main=sprintf("RNA/DNA %s vs %s\n(p-rho %.2f, s-rho %.2f)",selSamples[1],selSamples[2], cor(log2(res1$ratio.x),log2(res1$ratio.y),method="pearson"), cor(res1$ratio.x,res1$ratio.y,method="spearman")), xlab=sprintf("%s log2 RNA/DNA per insert",selSamples[1]), ylab=sprintf("%s log2 RNA/DNA per insert",selSamples[2]), pch=19, cex=0.3, col = "darkred", las=1, cex.axis=1.3, cex.lab=1.3, cex.main=1.2)
  abline(0,1,lty=2)

  dev.off()
 }
}

#-------------------------------------------

Rscript perInsertCounts_correlation.R assignment_counts_HAR_UMIs
# [1] "WTC-NPC"
# [1] "RNA"
# [1] 0.9982328
# [1] 0.9967485
# [1] 0.9966782
# [1] "DNA"
# [1] 0.9963769
# [1] 0.9944288
# [1] 0.9963273
# [1] "ratio"
# [1] 0.9038816
# [1] 0.8158461
# [1] 0.8323594
# [1] "NormSymmetry"
# [1] 144
# [1] "-----"
# [1] "Pt2a-NPC"
# [1] "RNA"
# [1] 0.9986555
# [1] 0.9986124
# [1] 0.9983381
# [1] "DNA"
# [1] 0.9953366
# [1] 0.9949083
# [1] 0.9974046
# [1] "ratio"
# [1] 0.9185729
# [1] 0.9052699
# [1] 0.9094181
# [1] "NormSymmetry"
# [1] 78
# [1] "-----"

R --vanilla --quiet
#-----------------------------------------

for (folder in c('assignment_counts_HAR_UMIs')) 
{
 for (sample in c('WTC-NPC','Pt2a-NPC')) {
  samples <- c(sprintf('%s-1',sample),sprintf('%s-2',sample),sprintf('%s-3',sample))

  png(sprintf("%s/%s-byInsert.png",folder,sample),width=800,height=800,type="cairo")
  par(mfrow=c(3,3))

  for (selSamples in combn(samples,2,simplify=F))
  {
    data1 <- read.table(sprintf("%s/%s-byInsert.tsv",folder,selSamples[1]),as.is=T,sep="\t",header=T)
    data2 <- read.table(sprintf("%s/%s-byInsert.tsv",folder,selSamples[2]),as.is=T,sep="\t",header=T)
    res1 <- merge(data1,data2,by="name")
    plot(res1$DNA.x,res1$DNA.y, main=sprintf("DNA %s vs %s\n(p-rho %.2f, s-rho %.2f)",selSamples[1],selSamples[2], cor(res1$DNA.x,res1$DNA.y,method="pearson"), cor(res1$DNA.x,res1$DNA.y,method="spearman")), xlab=sprintf("%s normalized DNA count per insert",selSamples[1]), ylab=sprintf("%s normalized DNA count per insert",selSamples[2]), pch=19, cex=0.3, col = "darkblue", las=1, cex.axis=1.3, cex.lab=1.3, cex.main=1.2)
    abline(0,1,lty=2)
    plot(res1$RNA.x,res1$RNA.y, main=sprintf("RNA %s vs %s\n(p-rho %.2f, s-rho %.2f)",selSamples[1],selSamples[2], cor(res1$RNA.x,res1$RNA.y,method="pearson"), cor(res1$RNA.x,res1$RNA.y,method="spearman")), xlab=sprintf("%s normalized RNA count per insert",selSamples[1]), ylab=sprintf("%s normalized RNA count per insert",selSamples[2]), pch=19, cex=0.3, col = "darkgreen", las=1, cex.axis=1.3, cex.lab=1.3, cex.main=1.2)
    abline(0,1,lty=2)
    plot(log2(res1$ratio.x),log2(res1$ratio.y), main=sprintf("RNA/DNA %s vs %s\n(p-rho %.2f, s-rho %.2f)",selSamples[1],selSamples[2], cor(log2(res1$ratio.x),log2(res1$ratio.y),method="pearson"), cor(res1$ratio.x,res1$ratio.y,method="spearman")), xlab=sprintf("%s log2 RNA/DNA per insert",selSamples[1]), ylab=sprintf("%s log2 RNA/DNA per insert",selSamples[2]), pch=19, cex=0.3, col = "darkred", las=1, cex.axis=1.3, cex.lab=1.3, cex.main=1.2)
    abline(0,1,lty=2)
  }
  dev.off()
 }
}

#-----------------------------------------

R --vanilla --quiet
#-----------------------------------------

for (folder in c('assignment_counts_HAR_UMIs')) 
{
 for (sample in c('WTC-NPC','Pt2a-NPC')) {
  print(sample)
  png(sprintf("%s/%s-TagCountsbyInsert.png",folder,sample),width=900,height=350,type="cairo")
  par(mfrow=c(1,3))
  samples <- c(sprintf('%s-1',sample),sprintf('%s-2',sample),sprintf('%s-3',sample))
  for (selSample in samples)
  {
    data1 <- read.table(sprintf("%s/%s-byInsert.tsv",folder,selSample),as.is=T,sep="\t",header=T)
    print(summary(data1$tags))
    res <- as.data.frame(table(data1$tags),stringsAsFactors=F)
    res$Var1 <- as.numeric(res$Var1)
    plot(res,xlim=c(0,100),type='h',ylab="Frequency",xlab="Number of tags per insert",main=sprintf("%s: %d inserts",selSample,dim(data1)[1]), las=1, cex.axis=1.3, cex.lab=1.3, cex.main=1.2)
  }
  dev.off()
  print("-----------")
 }
}
# [1] "WTC-NPC"
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   17.00   76.00   82.00   81.59   88.00  100.00 
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   17.00   76.00   82.00   81.46   88.00  100.00 
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   17.00   73.00   80.00   79.63   86.00  100.00 
# [1] "-----------"
# [1] "Pt2a-NPC"
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   18.00   75.00   82.00   81.38   88.00  100.00 
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   17.00   76.00   82.00   81.64   88.00  100.00 
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   17.00   75.00   81.00   81.21   88.00  100.00 
# [1] "-----------"

#-----------------------------------------

##############################################
# What controls were included in the design? #
##############################################

zcat ~/regulatory_tests/HARs/design/design_9.2.2015_HAR_library_2440_seqs.tsv.gz | grep ctrl | sort | rev | cut -f 2- -d':' | rev | uniq -c
    100 brain_vista_pos_ctrl10__cut=1_of_1__length=171__range_inclusive=1_to_171__VALIDATED_STRONG_CONTROL__no_modifications__NO_EcoRI_RESTRICTION_SITES_MODIFIED__NO_SbfI_RESTRICTION_SITES_MODIFIED__
    100 ENCODE_neg_ctrl42__cut=1_of_1__length=171__range_inclusive=1_to_171__VALIDATED_STRONG_CONTROL__no_modifications__NO_EcoRI_RESTRICTION_SITES_MODIFIED__NO_SbfI_RESTRICTION_SITES_MODIFIED__
    100 ENCODE_pos_ctrl5__cut=1_of_1__length=171__range_inclusive=1_to_171__VALIDATED_STRONG_CONTROL__no_modifications__NO_EcoRI_RESTRICTION_SITES_MODIFIED__NO_SbfI_RESTRICTION_SITES_MODIFIED__
    100 N1_K27me3_neg_ctrl42__cut=1_of_2__length=171__range_inclusive=1_to_171__VALIDATED_STRONG_CONTROL__no_modifications__NO_EcoRI_RESTRICTION_SITES_MODIFIED__NO_SbfI_RESTRICTION_SITES_MODIFIED__
    
R --vanilla --quiet
#-----------------------------------------

VISTA_brain_pos <- "brain_vista_pos_ctrl10__cut=1_of_1__length=171__range_inclusive=1_to_171__VALIDATED_STRONG_CONTROL__no_modifications__NO_EcoRI_RESTRICTION_SITES_MODIFIED__NO_SbfI_RESTRICTION_SITES_MODIFIED__"
ENCODE_ctrl42_neg <- "ENCODE_neg_ctrl42__cut=1_of_1__length=171__range_inclusive=1_to_171__VALIDATED_STRONG_CONTROL__no_modifications__NO_EcoRI_RESTRICTION_SITES_MODIFIED__NO_SbfI_RESTRICTION_SITES_MODIFIED__"
ENCODE_ctrl5_pos <- "ENCODE_pos_ctrl5__cut=1_of_1__length=171__range_inclusive=1_to_171__VALIDATED_STRONG_CONTROL__no_modifications__NO_EcoRI_RESTRICTION_SITES_MODIFIED__NO_SbfI_RESTRICTION_SITES_MODIFIED__"
N1_K27me3_ctrl42_neg <- "N1_K27me3_neg_ctrl42__cut=1_of_2__length=171__range_inclusive=1_to_171__VALIDATED_STRONG_CONTROL__no_modifications__NO_EcoRI_RESTRICTION_SITES_MODIFIED__NO_SbfI_RESTRICTION_SITES_MODIFIED__"

for (folder in c('assignment_counts_HAR_UMIs'))
{
  for (sample in c('WTC-NPC','Pt2a-NPC')) {
  png(sprintf("%s/%s-controls.png",folder,sample),width=700,height=1200,type="cairo")
  par(mfrow=c(3,1))

  for (i in 1:3)
  {
    RNA <- read.table(sprintf("%s/%s-%d-RNA.tsv.gz",folder,sample,i),as.is=T)
    DNA <- read.table(sprintf("%s/%s-%d-DNA.tsv.gz",folder,sample,i),as.is=T)

    data <- merge(RNA,DNA,by="V1")[,-3]
    colnames(data) <- c("barcodes","X","Y","name")
    data$category <- sapply(strsplit(data$name,":",fixed=T),FUN=function(x){ paste(x[1:(length(x)-1)],sep="",collapse=":") })
    data$ratio <- (data$X/(sum(data$X)*10^6))/(data$Y/(sum(data$Y)*10^6))

    res <- as.data.frame(t(sapply(unique(data$category),FUN=function(x) { sel <- which(data$category == x); c((sum(data$X[sel])/length(sel))/sum(data$X)*10^6,(sum(data$Y[sel])/length(sel))/sum(data$Y)*10^6,length(sel)) } )))
    res$name <- rownames(res)
    res$ratio <- res$V1/res$V2

    boxplot(list(VISTApos=log2(data$ratio[grep(VISTA_brain_pos,data$category)]),ENCODEpos=log2(data$ratio[grep(ENCODE_ctrl5_pos,data$category)]),ENCODEneg=log2(data$ratio[grep(ENCODE_ctrl42_neg,data$category)]),N1neg=log2(data$ratio[grep(N1_K27me3_ctrl42_neg,data$category)])),horizontal=T,main=sprintf("%s %d",sample,i),xlab="Log2 RNA/DNA ratio")

    points(log2(res$ratio[grep(VISTA_brain_pos,res$name)]),1,pch=19,cex=2,col="red")
    points(log2(res$ratio[grep(ENCODE_ctrl5_pos,res$name)]),2,pch=19,cex=2,col="red")
    points(log2(res$ratio[grep(ENCODE_ctrl42_neg,res$name)]),3,pch=19,cex=2,col="red")
    points(log2(res$ratio[grep(N1_K27me3_ctrl42_neg,res$name)]),4,pch=19,cex=2,col="red")
  }
  dev.off()
 }
}

#-----------------------------------------
