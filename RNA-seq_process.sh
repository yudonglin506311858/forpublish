#人类
for i in {hs_proerythroblast_1,hs_proerythroblast_2,hs_proerythroblast_3,hs_early_basophilic_1,hs_early_basophilic_2,hs_early_basophilic_3,hs_late_basophilic_1,hs_late_basophilic_2,hs_late_basophilic_3,hs_polychromatic_1,hs_polychromatic_2,hs_polychromatic_3,hs_orthochromatic_1,hs_orthochromatic_2,hs_orthochromatic_3};
do 
cutadapt --minimum-length 15 --max-n 8 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o ${i}_rmadp.fastq.gz ${i}.fastq.gz >>filter.txt 2>&1
java -jar /data/yudonglin/software/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 10 -phred33 ${i}_rmadp.fastq.gz ${i}_fliter.fq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
rm ${i}_rmadp.fastq.gz
#hisat2比对到基因组：
hisat2 --thread 10 -x /data/yudonglin/reference/hg19/genome -U ${i}_fliter.fq.gz -S ${i}.sam >>hisat2.txt 2>&1
rm ${i}_fliter.fq.gz
#SAM to BAM
samtools sort -@ 30 -o ${i}.bam ${i}.sam >>samtools.txt 2>&1
rm ${i}.sam
#基因表达可以用featureCounts
featureCounts -T 10 -t exon -g gene_id -a /data/yudonglin/reference/hg19/Homo_sapiens.GRCh37.75.gtf -o ${i}.count ${i}.bam >>count.txt 2>&1
samtools sort ${i}.bam -o ${i}_sorted.bam;
samtools index ${i}_sorted.bam;
conda activate atac-seq
bamCoverage --bam ${i}_sorted.bam -o ${i}_sorted.bam.bw  --binSize 10 ;
rm ${i}_sorted.bam;
conda deactivate
#rm ${i}.bam;
done

#小鼠
for i in {mm_proerythroblast_1,mm_proerythroblast_2,mm_proerythroblast_3,mm_basophilic_1,mm_basophilic_2,mm_basophilic_3,mm_polychromatic_1,mm_polychromatic_2,mm_polychromatic_3,mm_orthochromatic_1,mm_orthochromatic_2,mm_orthochromatic_3};
do 
cutadapt --minimum-length 15 --max-n 8 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o ${i}_rmadp.fastq.gz ${i}.fastq.gz >>filter.txt 2>&1
java -jar /data/yudonglin/software/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 10 -phred33 ${i}_rmadp.fastq.gz ${i}_fliter.fq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
rm ${i}_rmadp.fastq.gz
#hisat2比对到基因组：
hisat2 --thread 10 -x /data/yudonglin/reference/mm10/genome -U ${i}_fliter.fq.gz -S ${i}.sam >>hisat2.txt 2>&1
rm ${i}_fliter.fq.gz
#SAM to BAM
samtools sort -@ 30 -o ${i}.bam ${i}.sam >>samtools.txt 2>&1
rm ${i}.sam
#基因表达可以用featureCounts
featureCounts -T 10 -t exon -g gene_id -a /data/yudonglin/reference/Mus_musculus.GRCm38.102.gtf -o ${i}.count ${i}.bam >>count.txt 2>&1
samtools sort ${i}.bam -o ${i}_sorted.bam;
samtools index ${i}_sorted.bam;
conda activate atac-seq
bamCoverage --bam ${i}_sorted.bam -o ${i}_sorted.bam.bw  --binSize 10 ;
rm ${i}_sorted.bam;
conda deactivate
#rm ${i}.bam;
done
















































#第二篇cell report的RNA-seq数据

#RNA-seq
for((i=71;i<=98;i++));
do /data/yudonglin/software/ena-fast-download/ena-fast-download.py SRR72954${i};
done

#重新命名
mv SRR7295471_1.fastq.gz Donor1_P1_Rep1_RNA_1.fastq.gz
mv SRR7295471_2.fastq.gz Donor1_P1_Rep1_RNA_2.fastq.gz

mv SRR7295472_1.fastq.gz Donor1_P1_Rep2_RNA_1.fastq.gz
mv SRR7295472_2.fastq.gz Donor1_P1_Rep2_RNA_2.fastq.gz

mv SRR7295473_1.fastq.gz Donor1_P1_Rep3_RNA_1.fastq.gz
mv SRR7295473_2.fastq.gz Donor1_P1_Rep3_RNA_2.fastq.gz

mv SRR7295474_1.fastq.gz Donor1_P1_Rep4_RNA_1.fastq.gz
mv SRR7295474_2.fastq.gz Donor1_P1_Rep4_RNA_2.fastq.gz

mv SRR7295475_1.fastq.gz Donor1_P1_Rep5_RNA_1.fastq.gz
mv SRR7295475_2.fastq.gz Donor1_P1_Rep5_RNA_2.fastq.gz

mv SRR7295476_1.fastq.gz Donor1_P1_Rep6_RNA_1.fastq.gz
mv SRR7295476_2.fastq.gz Donor1_P1_Rep6_RNA_2.fastq.gz

mv SRR7295477_1.fastq.gz Donor1_P1_Rep7_RNA_1.fastq.gz
mv SRR7295477_2.fastq.gz Donor1_P1_Rep7_RNA_2.fastq.gz

mv SRR7295478_1.fastq.gz Donor1_P1_Rep8_RNA_1.fastq.gz
mv SRR7295478_2.fastq.gz Donor1_P1_Rep8_RNA_2.fastq.gz

mv SRR7295479_1.fastq.gz Donor1_P5_Rep2_RNA_1.fastq.gz
mv SRR7295479_2.fastq.gz Donor1_P5_Rep2_RNA_2.fastq.gz

mv SRR7295480_1.fastq.gz Donor1_P6_Rep2_RNA_1.fastq.gz
mv SRR7295480_2.fastq.gz Donor1_P6_Rep2_RNA_2.fastq.gz

mv SRR7295481_1.fastq.gz Donor1_P7_Rep2_RNA_1.fastq.gz
mv SRR7295481_2.fastq.gz Donor1_P7_Rep2_RNA_2.fastq.gz

mv SRR7295482_1.fastq.gz Donor1_P8_Rep2_RNA_1.fastq.gz
mv SRR7295482_2.fastq.gz Donor1_P8_Rep2_RNA_2.fastq.gz

mv SRR7295483_1.fastq.gz Donor2_P1_Rep2_RNA_1.fastq.gz
mv SRR7295483_2.fastq.gz Donor2_P1_Rep2_RNA_2.fastq.gz

mv SRR7295484_1.fastq.gz Donor2_P2_Rep2_RNA_1.fastq.gz
mv SRR7295484_2.fastq.gz Donor2_P2_Rep2_RNA_2.fastq.gz

mv SRR7295485_1.fastq.gz Donor2_P3_Rep2_RNA_1.fastq.gz
mv SRR7295485_2.fastq.gz Donor2_P3_Rep2_RNA_2.fastq.gz

mv SRR7295486_1.fastq.gz Donor2_P4_Rep2_RNA_1.fastq.gz
mv SRR7295486_2.fastq.gz Donor2_P4_Rep2_RNA_2.fastq.gz

mv SRR7295487_1.fastq.gz Donor3_P1_Rep3_RNA_1.fastq.gz
mv SRR7295487_2.fastq.gz Donor3_P1_Rep3_RNA_2.fastq.gz

mv SRR7295488_1.fastq.gz Donor3_P2_Rep3_RNA_1.fastq.gz
mv SRR7295488_2.fastq.gz Donor3_P2_Rep3_RNA_2.fastq.gz

mv SRR7295489_1.fastq.gz Donor3_P3_Rep3_RNA_1.fastq.gz
mv SRR7295489_2.fastq.gz Donor3_P3_Rep3_RNA_2.fastq.gz

mv SRR7295490_1.fastq.gz Donor3_P4_Rep3_RNA_1.fastq.gz
mv SRR7295490_2.fastq.gz Donor3_P4_Rep3_RNA_2.fastq.gz

mv SRR7295491_1.fastq.gz Donor3_P5_Rep3_RNA_1.fastq.gz
mv SRR7295491_2.fastq.gz Donor3_P5_Rep3_RNA_2.fastq.gz

mv SRR7295492_1.fastq.gz Donor3_P5_Rep4_RNA_1.fastq.gz
mv SRR7295492_2.fastq.gz Donor3_P5_Rep4_RNA_2.fastq.gz

mv SRR7295493_1.fastq.gz Donor3_P6_Rep3_RNA_1.fastq.gz
mv SRR7295493_2.fastq.gz Donor3_P6_Rep3_RNA_2.fastq.gz

mv SRR7295494_1.fastq.gz Donor3_P6_Rep4_RNA_1.fastq.gz
mv SRR7295494_2.fastq.gz Donor3_P6_Rep4_RNA_2.fastq.gz

mv SRR7295495_1.fastq.gz Donor3_P7_Rep3_RNA_1.fastq.gz
mv SRR7295495_2.fastq.gz Donor3_P7_Rep3_RNA_2.fastq.gz

mv SRR7295496_1.fastq.gz Donor3_P7_Rep4_RNA_1.fastq.gz
mv SRR7295496_2.fastq.gz Donor3_P7_Rep4_RNA_2.fastq.gz

mv SRR7295497_1.fastq.gz Donor3_P8_Rep3_RNA_1.fastq.gz
mv SRR7295497_2.fastq.gz Donor3_P8_Rep3_RNA_2.fastq.gz

mv SRR7295498_1.fastq.gz Donor3_P8_Rep4_RNA_1.fastq.gz
mv SRR7295498_2.fastq.gz Donor3_P8_Rep4_RNA_2.fastq.gz



for i in {Donor1_P1_Rep1_RNA,Donor1_P1_Rep6_RNA,Donor2_P3_Rep2_RNA,Donor3_P4_Rep3_RNA,Donor3_P7_Rep3_RNA,Donor1_P1_Rep1_RNA,Donor1_P1_Rep6_RNA,Donor2_P3_Rep2_RNA,Donor3_P4_Rep3_RNA,Donor3_P7_Rep3_RNA,Donor1_P1_Rep2_RNA,Donor1_P1_Rep8_RNA,Donor2_P4_Rep2_RNA,Donor3_P5_Rep3_RNA,Donor3_P7_Rep4_RNA,Donor1_P1_Rep2_RNA,Donor1_P1_Rep8_RNA,Donor2_P4_Rep2_RNA,Donor3_P5_Rep3_RNA,Donor3_P7_Rep4_RNA,Donor1_P1_Rep3_RNA,Donor1_P5_Rep2_RNA,Donor3_P1_Rep3_RNA,Donor3_P5_Rep4_RNA,Donor3_P8_Rep3_RNA,Donor1_P1_Rep3_RNA,Donor1_P5_Rep2_RNA,Donor3_P1_Rep3_RNA,Donor3_P5_Rep4_RNA,Donor3_P8_Rep3_RNA,Donor1_P1_Rep4_RNA,Donor1_P8_Rep2_RNA,Donor3_P2_Rep3_RNA,Donor3_P6_Rep3_RNA,Donor3_P8_Rep4_RNA,Donor1_P1_Rep4_RNA,Donor1_P8_Rep2_RNA,Donor3_P2_Rep3_RNA,Donor3_P6_Rep3_RNA,Donor3_P8_Rep4_RNA,Donor1_P1_Rep5_RNA,Donor2_P2_Rep2_RNA,Donor3_P3_Rep3_RNA,Donor3_P6_Rep4_RNA,Donor1_P1_Rep5_RNA,Donor2_P2_Rep2_RNA,Donor3_P3_Rep3_RNA,Donor3_P6_Rep4_RNA};
do 
cutadapt --pair-filter=any --minimum-length 15 --max-n 8 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o ${i}_rmadp_1.fastq.gz -p ${i}_rmadp_2.fastq.gz ${i}_1.fastq.gz ${i}_2.fastq.gz >>filter.txt 2>&1
java -jar /data/yudonglin/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 20 -phred33 ${i}_rmadp_1.fastq.gz ${i}_rmadp_2.fastq.gz -baseout ${i}_fliter.fq.gz  AVGQUAL:20 SLIDINGWINDOW:4:15 MINLEN:15 1>>filter.txt 2>&1;
rm ${i}_rmadp_1.fastq.gz
rm ${i}_rmadp_2.fastq.gz
#hisat2比对到基因组：
hisat2 --thread 10 -x /data/yudonglin/reference/hg19/genome -U ${i}_fliter.fq.gz -S ${i}.sam >>hisat2.txt 2>&1
rm ${i}_fliter.fq.gz
#SAM to BAM
samtools sort -@ 30 -o ${i}.bam ${i}.sam >>samtools.txt 2>&1
rm ${i}.sam
#基因表达可以用featureCounts
featureCounts -T 10 -t exon -g gene_id -a /data/yudonglin/reference/Homo_sapiens.GRCh37.75.gtf -o ${i}.count ${i}.bam >>count.txt 2>&1
samtools sort ${i}.bam -o ${i}_sorted.bam;
samtools index ${i}_sorted.bam;
conda activate atac-seq
bamCoverage --bam ${i}_sorted.bam -o ${i}_sorted.bam.bw  --binSize 10 ;
rm ${i}_sorted.bam;
conda deactivate
#rm ${i}.bam;
done



