##质控

ls | while read id; do fastqc -t 8 -o /home/zhushu/zhouty19990625/download/rxx_data/fastqc1/ $id; done
multiqc *fastqc.zip --pdf -o  /home/zhushu/zhouty19990625/download/rxx_data/output

##修剪

cat ~/test.txt | while read id; do java -jar /home/zhushu/zhouty19990625/bioinfo_tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE "$id".R1.fastq.gz "$id".R2.fastq.gz -baseout /home/zhushu/zhouty19990625/download/rxx_data/rxx_after_trimming_fastq/"$id".fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:true SLIDINGWINDOW:4:15 MINLEN:50; done

##质控

ls | while read id; do fastqc -t 8 -o /home/zhushu/zhouty19990625/download/rxx_data/fastqc2/ $id; done
multiqc *fastqc.zip --pdf -o  /home/zhushu/zhouty19990625/download/rxx_data/output

##比对

STAR  --runMode genomeGenerate --genomeDir ~/genome/mouse/index/ --runThreadN 8 --genomeFastaFiles ~/genome/mouse/dna/Mus_musculus.GRCm38.dna_rm.primary_assembly.fa --sjdbGTFfile ~/genome/mouse/ann/Mus_musculus.GRCm38.102.gtf --sjdbOverhang 148

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

##质控（用RSeQC）

geneBody_coverage.py -r /home/zhushu/zhouty19990625/genome/mouse/ann_bed/mm10_RefSeq.bed -i ~/download/rxx_data/bam_path.txt  -o ~/download/rxx_data/output/

gtfToGenePred -ignoreGroupsWithoutExons ~/genome/mouse/ann/Mus_musculus.GRCm38.102.gtf ~/genome/mouse/ann_bed/Mus_musculus.GRCm38.102.gtf.50505050.pred
genePredToBed ~/genome/mouse/ann_bed/Mus_musculus.GRCm38.102.gtf.50505050.pred ~/genome/mouse/ann_bed/Mus_musculus.GRCm38.102.gtf.bed12

awk '$3-$2>1000 && $3-$2<2000' Mus_musculus.GRCm38.102.gtf.bed12 >Mus_musculus.GRCm38.model.gtf.bed12

cat ~/test_rxx.txt|while read id; do echo ${id}; cd ~/download/rxx_data/${id}/; bam_stat.py  -i ${id}.Aligned.sortedByCoord.out.bam; read_distribution.py -i ${id}.Aligned.sortedByCoord.out.bam -r /home/zhushu/zhouty19990625/genome/mouse/ann_bed/Mus_musculus.GRCm38.102.gtf.bed12; RPKM_saturation.py -i ${id}.Aligned.sortedByCoord.out.bam -r /home/zhushu/zhouty19990625/genome/mouse/ann_bed/Mus_musculus.GRCm38.102.gtf.bed12 -s 10 -q 0 -o ~/download/rxx_data/output_test/${id}.RPKM_saturation;done

##整理表达矩阵

mkdir count

cat ~/zwq_dir.txt |while read id;do
cd ${id};
awk '{print $1,$3}' ${id}.ReadsPerGene.out.tab > ~/download/zwq/count/${id}.count;
cd /home/zhushu/zhouty19990625/download/zwq;
done


perl -lne 'if ($ARGV=~/(.*).count/){print "$1\t$_"}' *.count >matrix.count
