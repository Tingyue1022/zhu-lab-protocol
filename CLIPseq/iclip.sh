cat ~/n6clip.txt|while read id; do cd ~/download/N6clip/${id}; fastqc -t 8 -o ~/download/N6clip/fastqc1/ ${id}_1.fq.gz;done 
cat ~/n6clip.txt|while read id; do cd ~/download/N6clip/${id}; fastqc -t 8 -o ~/download/N6clip/fastqc1/ ${id}_2.fq.gz;done
cd ~/download/N6clip/fastqc1/
multiqc *fastqc.zip --pdf -o  /home/zhushu/zhouty19990625/download/N6clip/output/


cat ~/n6clip.txt | while read id;do cd ~/download/N6clip/${id}; java -jar /home/zhushu/zhouty19990625/bioinfo_tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE ${id}_1.fq.gz ${id}_2.fq.gz -baseout /home/zhushu/zhouty19990625/download/N6clip/iclip_after_trimming_fastq/${id}.fastq ILLUMINACLIP:/home/zhushu/zhouty19990625/bioinfo_tools/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:10 MINLEN:18; done

cat ~/n6clip.txt|while read id; do echo ${id}; bowtie2 -p 8 -x ~/genome/EMCV/EMCV_complete_genome -1 /home/zhushu/zhouty19990625/download/N6clip/iclip_after_trimming_fastq/${id}_1P.fastq -2 /home/zhushu/zhouty19990625/download/N6clip/iclip_after_trimming_fastq/${id}_2P.fastq -S ~/download/N6clip/align/${id}.sam; done

cd /home/zhushu/zhouty19990625/download/N6clip/align/
ls *.sam|while read id; do samtools view -S ${id} -b > ${id}.bam;samtools sort ${id}.bam -o ${id}.sorted.bam; samtools index ${id}.sorted.bam; done

Piranha N6-CLIP-emcvEV-1_FKDL192533959-1a-13.sam.sorted.bam