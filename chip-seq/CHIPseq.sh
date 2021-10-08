cd ~/bioinfo_tools/linuxnd
./lnd login -u X101SC21043714-Z01-F001 -p ba4c38ct
./lnd cp -d oss://CP2021042900040/H101SC21043714/RSCS0500/X101SC21043714-Z01/X101SC21043714-Z01-F001/2.cleandata ~/download/chip_rxx

cd /home/zhushu/zhouty19990625/download/chip_rxx/2.cleandata/15_BKDL210023553-1a
bowtie2 -p 8 -3 5 --local -x ~/genome/human/index_bowtie/hs -1 15_BKDL210023553-1a_1.clean.fq.gz -2 15_BKDL210023553-1a_2.clean.fq.gz -S ~/download/chip_rxx/align/chip15.sam

cd /home/zhushu/zhouty19990625/download/chip_rxx/2.cleandata/16_BKDL210023554-1a
bowtie2 -p 8 -3 5 --local -x ~/genome/human/index_bowtie/hs -1 16_BKDL210023554-1a_1.clean.fq.gz -2 16_BKDL210023554-1a_2.clean.fq.gz -S ~/download/chip_rxx/align/chip16.sam

cd /home/zhushu/zhouty19990625/download/chip_rxx/align/
ls *.sam|while read id; do samtools view -S ${id} -b > ${id}.bam;samtools sort ${id}.bam -o ${id}.sorted.bam; samtools index ${id}.sorted.bam; done


cd ~/download/chip_rxx/align/
macs3 callpeak -t chip15.sam.sorted.bam -f BAM -g hs -n chip15 -B -q 0.05
macs3 callpeak -t chip16.sam.sorted.bam -f BAM -g hs -n chip16 -B -q 0.05








####chipseq_version2.0
#datadownload
cd ~/bioinfo_tools/linuxnd
./lnd login -u X101SC21043714-Z01-J003 -p g06h1a5e
./lnd cp -d oss://CP2021042900040/H101SC21043714/RSCS0500/X101SC21043714-Z01/X101SC21043714-Z01-J003/2.cleandata ~/download/chip_rxx

#align
##index
cd /home/zhushu/zhouty19990625/genome/mouse/dna/
bowtie2-build Mus_musculus.GRCm38.dna_rm.primary_assembly.fa mm
ls *.bt2
cd ~/download/chip_hkx/2.cleandata/
ls|while read id;do cd ${id};bowtie2 -p 8 -3 5 --local -x ~/genome/human/index_bowtie/hs -1 ${id}_1.clean.fq.gz -2 ${id}_2.clean.fq.gz -S ~/download/chip_hkx/align/${id}.sam;cd ~/download/chip_hkx/2.cleandata;done
cd /home/zhushu/zhouty19990625/download/chip_hkx/align/
ls *.sam|while read id; do samtools view -S ${id} -b > ${id}.bam;samtools sort ${id}.bam -o ${id}.sorted.bam; samtools index ${id}.sorted.bam; done

#call peak
conda activate rna-seq


cd ~/download/chip_hkx/align/

macs3 callpeak -t 2-n1_BKDL210030916-1a.sam.sorted.bam \
               -c 1-ev_BKDL210030915-1a.sam.sorted.bam \
               -f BAM -g hs -n 2-n1_BKDL210030916-1a -B -q 0.01 --outdir ~/download/chip_hkx/macs

macs3 callpeak -t 3-132_BKDL210030917-1a.sam.sorted.bam \
               -c 1-ev_BKDL210030915-1a.sam.sorted.bam \
               -f BAM -g hs -n 3-132_BKDL210030917-1a -B -q 0.01 --outdir ~/download/chip_hkx/macs

macs3 callpeak -t 4-fl_BKDL210030918-1a.sam.sorted.bam \
               -c 1-ev_BKDL210030915-1a.sam.sorted.bam \
               -f BAM -g hs -n 4-fl_BKDL210030918-1a -B -q 0.01 --outdir ~/download/chip_hkx/macs

cd ~/download/chip_hkx/macs
ls *.bdg|while read id; do sort -k1,1 -k2,2n ${id} > ${id}.sort.bdg;done


cd /home/zhushu/zhouty19990625/genome/human/dna/
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
cut -f 1,2 Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai > GRCh38.chrom.sizes

cd ~/download/chip_hkx/macs
ls *.bdg|while read id; do bedGraphToBigWig ${id} /home/zhushu/zhouty19990625/genome/human/dna/GRCh38.chrom.sizes ${id}.bw;done
cd ~/download/chip_rxx/macs
ls *control_lambda.bdg|while read id; do sort -k1,1 -k2,2n ${id} > ${id}.sort.bdg;done
ls *control_lambda.bdg.sort.bdg|while read id; do bedGraphToBigWig ${id} /home/zhushu/zhouty19990625/genome/mouse/dna/mm10.chrom.sizes ${id}.bw;done

#vis
track type=bigWig name="2-n1" description="2-n1" visibility=2 color=255,0,0 bigDataUrl=https://data.cyverse.org/dav-anon/iplant/home/zhouty19990625/2-n1_BKDL210030916-1a_treat_pileup.bdg.sort.bdg.bw


track type=bigWig name="3-132" description="3-132" visibility=2 color=255,205,0 bigDataUrl=https://data.cyverse.org/dav-anon/iplant/home/zhouty19990625/3-132_BKDL210030917-1a_treat_pileup.bdg.sort.bdg.bw



track type=bigWig name="4-fl" description="4-fl" visibility=2 color=0,176,80 bigDataUrl=https://data.cyverse.org/dav-anon/iplant/home/zhouty19990625/4-fl_BKDL210030918-1a_treat_pileup.bdg.sort.bdg.bw



track type=bigWig name="2-n1_c" description="2-n1_c" visibility=2 color=255,255,255 bigDataUrl=https://data.cyverse.org/dav-anon/iplant/home/zhouty19990625/2-n1_BKDL210030916-1a_control_lambda.bdg.sort.bdg.bw



track type=bigWig name="3-132_c" description="3-132_c" visibility=2 color=255,255,255 bigDataUrl=https://data.cyverse.org/dav-anon/iplant/home/zhouty19990625/3-132_BKDL210030917-1a_control_lambda.bdg.sort.bdg.bw





track type=bigWig name="4-fl_c" description="4-fl_c" visibility=2 color=255,255,255 bigDataUrl=https://data.cyverse.org/dav-anon/iplant/home/zhouty19990625/4-fl_BKDL210030918-1a_control_lambda.bdg.sort.bdg.bw

#end of the run 


#以下部分展示使用chipseeker展示chipseq结果
#首先我会把bed文件下载到本地

awk 'OFS="\t" {if (NR >0 )$1="chr"$1; print}' 2D9_BKDL210030912-1a_summits.bed > D9.txt
$r












