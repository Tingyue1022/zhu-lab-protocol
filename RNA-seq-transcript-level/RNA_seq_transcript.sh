##质控
#cd /home/zhushu/zhouty19990625/download/rxx/
#ls | while read id; do fastqc -t 8 -o /home/zhushu/zhouty19990625/download/rxx_data/fastqc1/ $id; done
#cd /home/zhushu/zhouty19990625/download/rxx_data/fastqc1/
#multiqc *fastqc.zip --pdf -o  /home/zhushu/zhouty19990625/download/rxx_data/output
#cd /home/zhushu/zhouty19990625/download/rxx

##trimmomatic修剪
cd /home/zhushu/zhouty19990625/download/rxx
cat ~/test.txt | while read id; do java -jar /home/zhushu/zhouty19990625/bioinfo_tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE "$id".R1.fastq.gz "$id".R2.fastq.gz -baseout /home/zhushu/zhouty19990625/download/rxx_data/rxx_after_trimming_fastq/"$id".fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:true SLIDINGWINDOW:4:15 MINLEN:50; done

#cd /home/zhushu/zhouty19990625/download/rxx_data/rxx_after_trimming_fastq/
#ls | while read id; do fastqc -t 8 -o /home/zhushu/zhouty19990625/download/rxx_data/fastqc2/ $id; done
#cd /home/zhushu/zhouty19990625/download/rxx_data/fastqc2/
#multiqc *fastqc.zip --pdf -o  /home/zhushu/zhouty19990625/download/rxx_data/output
cd ~
##使用star建立索引，这里我用测序数据的读长是149bp，所以 --sjdbOverhang这个选项设置成148
#STAR  --runMode genomeGenerate --genomeDir ~/genome/mouse/index/ --runThreadN 8 --genomeFastaFiles ~/genome/mouse/dna/Mus_musculus.GRCm38.dna_rm.primary_assembly.fa --sjdbGTFfile ~/genome/mouse/ann/Mus_musculus.GRCm38.102.gtf --sjdbOverhang 148



cd /home/zhushu/zhouty19990625/Duerkop_2018/fastq/
ls | while read id; do fastqc -t 8 -o /home/zhushu/zhouty19990625/Duerkop_2018/fastqc/ $id; done
cd /home/zhushu/zhouty19990625/Duerkop_2018/fastqc/
multiqc *fastqc.zip --pdf -o ~/download/


#cd ~/download/rxx_data/rxx_after_trimming_fastq/
##我想把trimmomatic质控之后的pair和unpair的两个.fastq文件合成成一个然后直接去跑star比对，但是发现不行，见“trim.txt里面记录的报错，于是直接用原始数据跑star。
#cat ~/test_rxx.txt| while read id; do cat ${id}_1P.fastq ${id}_1U.fastq > ~/download/rxx_data/rxx_after_trimming_merge/${id}_1M.fastq; cat ${id}_2P.fastq ${id}_2U.fastq > ~/download/rxx_data/rxx_after_trimming_merge/${id}_2M.fastq; done


##run star
max_intron_size=1000000
star_p=" --outFilterType BySJout --outSAMattributes NH HI AS NM MD \
       --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
       --alignIntronMin 20 --alignIntronMax ${max_intron_size} \
                --alignMatesGapMax ${max_intron_size} \
                --outFilterMatchNminOverLread 0.66 --outFilterScoreMinOverLread 0.66 \
                --winAnchorMultimapNmax 70 --seedSearchStartLmax 45 \
                --outSAMattrIHstart 0 --outSAMstrandField intronMotif \
                --genomeLoad NoSharedMemory --outReadsUnmapped Fastx \
                --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts"

cat ~/test_rxx.txt |while read id; do
cd ~/download/rxx/;
              
STAR --runMode alignReads --runThreadN 8 \
         --readFilesIn ${id}.R1.fastq.gz ${id}.R2.fastq.gz \
         --genomeDir /home/zhushu/zhouty19990625/genome/mouse/index \
         --outFileNamePrefix ~/download/rxx_data/${id}/${id}.\
         --outBAMsortingThreadN 8 \
         --readFilesCommand zcat \
         ${star_p};
done

在设置好的输出文件夹里面会有这些文件：
#(base) [zhouty19990625@mgt S215WT_L2_337X37]$ cd ..; cd S316WT_L2_338X38; ls -lh
#总用量 8.3G
#-rw-r--r-- 1 zhouty19990625 zhushu  2.3G 2月  10 10:46 S316WT_L2_338X38.Aligned.sortedByCoord.out.bam
#-rw-r--r-- 1 zhouty19990625 zhushu  4.1G 2月  10 10:45 S316WT_L2_338X38.Aligned.toTranscriptome.out.bam
#-rw-r--r-- 1 zhouty19990625 zhushu  2.0K 2月  10 10:46 S316WT_L2_338X38.Log.final.out
#-rw-r--r-- 1 zhouty19990625 zhushu   17K 2月  10 10:46 S316WT_L2_338X38.Log.out
#-rw-r--r-- 1 zhouty19990625 zhushu  2.0K 2月  10 10:46 S316WT_L2_338X38.Log.progress.out
#-rw-r--r-- 1 zhouty19990625 zhushu  1.4M 2月  10 10:45 S316WT_L2_338X38.ReadsPerGene.out.tab
#-rw-r--r-- 1 zhouty19990625 zhushu  6.3M 2月  10 10:45 S316WT_L2_338X38.SJ.out.tab
#drwx------ 3 zhouty19990625 zhushu  4.0K 2月  10 10:46 S316WT_L2_338X38._STARtmp
#-rw-r--r-- 1 zhouty19990625 zhushu 1023M 2月  10 10:45 S316WT_L2_338X38.Unmapped.out.mate1
#-rw-r--r-- 1 zhouty19990625 zhushu 1023M 2月  10 10:45 S316WT_L2_338X38.Unmapped.out.mate2

##参见徐州更的教程，最常用的samtools三板斧就是格式转换，排序，索引：
#for i in `seq 56 58`
#do
#    samtools view -S SRR35899${i}.sam -b > SRR35899${i}.bam
#    samtools sort SRR35899${i}.bam -o SRR35899${i}_sorted.bam
#    samtools index SRR35899${i}_sorted.bam
#done
因为我用的star中有参数“--outSAMtype BAM SortedByCoordinate”，所以就不再使用samtool来完成格式转换和排序。使用samtools index选项之后，在设置好的输出文件夹里面会多个.bai文件，变成这样：
(base) [zhouty19990625@mgt rxx_data]$ cd S215WT_L2_337X37; ls -lh
#总用量 7.9G
#-rw-r--r-- 1 zhouty19990625 zhushu 2.2G 2月  10 10:29 S215WT_L2_337X37.Aligned.sortedByCoord.out.bam
#-rw-r--r-- 1 zhouty19990625 zhushu 2.0M 2月  10 17:11 S215WT_L2_337X37.Aligned.sortedByCoord.out.bam.bai
#-rw-r--r-- 1 zhouty19990625 zhushu 3.9G 2月  10 10:28 S215WT_L2_337X37.Aligned.toTranscriptome.out.bam
#-rw-r--r-- 1 zhouty19990625 zhushu 2.0K 2月  10 10:29 S215WT_L2_337X37.Log.final.out
#-rw-r--r-- 1 zhouty19990625 zhushu  17K 2月  10 10:29 S215WT_L2_337X37.Log.out
#-rw-r--r-- 1 zhouty19990625 zhushu 2.0K 2月  10 10:29 S215WT_L2_337X37.Log.progress.out
#-rw-r--r-- 1 zhouty19990625 zhushu 1.4M 2月  10 10:28 S215WT_L2_337X37.ReadsPerGene.out.tab
#-rw-r--r-- 1 zhouty19990625 zhushu 6.3M 2月  10 10:28 S215WT_L2_337X37.SJ.out.tab
#drwx------ 3 zhouty19990625 zhushu 4.0K 2月  10 10:29 S215WT_L2_337X37._STARtmp
#-rw-r--r-- 1 zhouty19990625 zhushu 967M 2月  10 10:28 S215WT_L2_337X37.Unmapped.out.mate1
#-rw-r--r-- 1 zhouty19990625 zhushu 967M 2月  10 10:28 S215WT_L2_337X37.Unmapped.out.mate2

##RSeQC进行质控

cd genome/mouse; mkdir ann_bed
wget  https://sourceforge.net/projects/rseqc/files/BED/Mouse_Mus_musculus/mm10_RefSeq.bed.gz
cat bam_path.txt
#/home/zhushu/zhouty19990625/download/rxx_data/S215WT_L2_337X37/S215WT_L2_337X37.Aligned.sortedByCoord.out.bam
#/home/zhushu/zhouty19990625/download/rxx_data/S316WT_L2_338X38/S316WT_L2_338X38.Aligned.sortedByCoord.out.bam
#/home/zhushu/zhouty19990625/download/rxx_data/S418KO_L2_339X39/S418KO_L2_339X39.Aligned.sortedByCoord.out.bam
#/home/zhushu/zhouty19990625/download/rxx_data/S519KO_L2_340X40/S519KO_L2_340X40.Aligned.sortedByCoord.out.bam
#/home/zhushu/zhouty19990625/download/rxx_data/S620KO_L2_341X41/S620KO_L2_341X41.Aligned.sortedByCoord.out.bam
#/home/zhushu/zhouty19990625/download/rxx_data/S723WT_L2_342X42/S723WT_L2_342X42.Aligned.sortedByCoord.out.bam

##mm10_RefSeq.bed从官网上下载的不行，bed文件和bam文件的格式不对应，还是应该gtfToGenePred，genePredToBed两步走把自己的gtf文件转化成bed
geneBody_coverage.py -r /home/zhushu/zhouty19990625/genome/mouse/ann_bed/mm10_RefSeq.bed -i ~/download/rxx_data/bam_path.txt  -o ~/download/rxx_data/output/

gtfToGenePred -ignoreGroupsWithoutExons ~/genome/mouse/ann/Mus_musculus.GRCm38.102.gtf ~/genome/mouse/ann_bed/Mus_musculus.GRCm38.102.gtf.50505050.pred
genePredToBed ~/genome/mouse/ann_bed/Mus_musculus.GRCm38.102.gtf.50505050.pred ~/genome/mouse/ann_bed/Mus_musculus.GRCm38.102.gtf.bed12

##这个Mus_musculus.GRCm38.102.gtf.bed12有17M，这么跑速度非常慢，过了一个多小时连一个bam文件都还没跑完。-o后面要加个前缀，下面的代码是错误的示范，没有前缀开头就是“."开头
geneBody_coverage.py -r /home/zhushu/zhouty19990625/genome/mouse/ann_bed/Mus_musculus.GRCm38.model.gtf.bed12 -i ~/download/rxx_data/bam_path.txt  -o ~/download/rxx_data/output/
##把17M的bed文件这样处理一下之后还有1M左右，再运行geneBody_coverage.py差不多10min一个bam，我看之后运行read_distribution.py和其他模块都是用的17M那个
awk '$3-$2>1000 && $3-$2<2000' Mus_musculus.GRCm38.102.gtf.bed12 >Mus_musculus.GRCm38.model.gtf.bed12

cat ~/test_rxx.txt|while read id; do echo ${id}; cd ~/download/rxx_data/${id}/; bam_stat.py  -i ${id}.Aligned.sortedByCoord.out.bam; read_distribution.py -i ${id}.Aligned.sortedByCoord.out.bam -r /home/zhushu/zhouty19990625/genome/mouse/ann_bed/Mus_musculus.GRCm38.102.gtf.bed12; RPKM_saturation.py -i ${id}.Aligned.sortedByCoord.out.bam -r /home/zhushu/zhouty19990625/genome/mouse/ann_bed/Mus_musculus.GRCm38.102.gtf.bed12 -s 10 -q 0 -o ~/download/rxx_data/output_test/${id}.RPKM_saturation;done


##结果发现我的数据不是链特异性的数据集，所以RPKM_saturation.py里面的 -d 就使用默认的选项。### 测序饱和度评估（RPKM_saturation.py）-s: 采样频率，0-100之间的整数，类似于步长，-q: 过滤低质量比对
cat ~/test_rxx.txt|while read id; do echo ${id};cd ~/download/rxx_data/${id}/; infer_experiment.py -r /home/zhushu/zhouty19990625/genome/mouse/ann_bed/Mus_musculus.GRCm38.102.gtf.bed12 -i ${id}.Aligned.sortedByCoord.out.bam; done

##转录本组装
cd ~/download/rxx_data/
cat ~/test_rxx.txt|while read id; do stringtie ${id}/${id}.Aligned.sortedByCoord.out.bam -G /home/zhushu/zhouty19990625/genome/mouse/ann/Mus_musculus.GRCm38.102.gtf -l ${id} -o ${id}/${id}.stringtie_first.gtf -f 0.01 -p 8; done
##合并
stringtie --merge -p 8 -G /home/zhushu/zhouty19990625/genome/mouse/ann/Mus_musculus.GRCm38.102.gtf -o stringtie_merged.gtf -l stringtie_merged mergelist.txt
#mergelist:cat ~/test_rxx.txt|while read id; do echo $(pwd)/${id}/${id}.stringtie_first.gtf;done

gffcompare -R -r /home/zhushu/zhouty19990625/genome/mouse/ann/Mus_musculus.GRCm38.102.gtf -o gffcmp ~/download/rxx_data/stringtie_merged.gtf

cat ~/lnc_dir.txt|while read id; do stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/${id}/${id}.gtf ~/download/linc1216/${id}/${id}.Aligned.sortedByCoord.out.bam;done


##r
library('ballgown')
pheno_data<-read.table('~/mhd_dir.txt',header = T)
AS<-ballgown(dataDir = "/Users/zhulab/Desktop/download/ballgown",samplePattern = "CD11b",pData = pheno_data)
save.image('/home/zhushu/zhouty19990625/download/mhd_seq/ballgown.Rdata')

##rmats
run_rmats --b1 ~/download/rxx_data/path_WT.txt --b2 ~/download/rxx_data/path_KO.txt --gtf /home/zhushu/zhouty19990625/genome/mouse/ann/Mus_musculus.GRCm38.102.gtf -t paired --readLength 150 --nthread 8 --od ~/download/rxx_data/WT_KO --tmp ~/download/rxx_data/tmp_WT_KO --tstat 8 --cstat 0.0001 --libType fr-unstranded

#awk 'NR == 1 || NR==6 ||NR==2' ~/test_rxx.txt|while read id; do echo $(pwd)/${id}/${id}.Aligned.sortedByCoord.out.bam;done>path_WT.txt
#awk 'NR == 3 || NR==4 ||NR==5' ~/test_rxx.txt|while read id; do echo $(pwd)/${id}/${id}.Aligned.sortedByCoord.out.bam;done>path_KO.txt
##上面的不对，跑出来报错，这软件要求.txt文件中bam路径写在一行上面，用逗号隔开

awk 'FNR==1 || $20<0.2' SE.MATS.JC.txt >SE.MATS.JC.sig.txt

rmats2sashimiplot --b1 /home/zhushu/zhouty19990625/download/rxx_data/S215WT_L2_337X37/S215WT_L2_337X37.Aligned.sortedByCoord.out.bam,/home/zhushu/zhouty19990625/download/rxx_data/S316WT_L2_338X38/S316WT_L2_338X38.Aligned.sortedByCoord.out.bam,/home/zhushu/zhouty19990625/download/rxx_data/S723WT_L2_342X42/S723WT_L2_342X42.Aligned.sortedByCoord.out.bam --b2 /home/zhushu/zhouty19990625/download/rxx_data/S418KO_L2_339X39/S418KO_L2_339X39.Aligned.sortedByCoord.out.bam,/home/zhushu/zhouty19990625/download/rxx_data/S519KO_L2_340X40/S519KO_L2_340X40.Aligned.sortedByCoord.out.bam,/home/zhushu/zhouty19990625/download/rxx_data/S620KO_L2_341X41/S620KO_L2_341X41.Aligned.sortedByCoord.out.bam -t SE -e ~/download/rxx_data/WT_KO/SE.MATS.JC.test.txt --exon_s 1 --intron_s 5 --l1 WT --l2 KO -o ~/download/rxx_data/WT_KO_plots

##bowtie alignment
bowtie2-build ~/genome/EMCV/EMCV_complete_genome.fasta ~/genome/EMCV/EMCV_complete_genome
cat ~/wam_test.txt|while read id; do echo ${id}; bowtie2 -x ~/genome/EMCV/EMCV_complete_genome -1 ~/download/wam/${id}_combined_R1.fastq.gz -2 ~/download/wam/${id}_combined_R2.fastq.gz -S ~/download/wam/${id}.sam; done




cat ~/test_wdc.txt|while read id; do cp -r /home/zhushu/zhouty19990625/download/wdc_data/${id}/${id}.Aligned.sortedByCoord.out.bam.bai out_bai;done
