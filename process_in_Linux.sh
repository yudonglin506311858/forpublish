
mkdir fastqc
#质量控制
fastqc -t 30 -o ./fastqc/ *.gz


#第一批的两个


#过滤接头
for i in {XH-Z2B654-Baso-19-11-26_FKDL192547546-1a-38,XH-Z2B654-Baso-19-11-27_FKDL192547546-1a-46,XH-Z2B654-ortho-19-11-26_FKDL192547546-1a-40,XH-Z2B654-ortho-19-11-27_FKDL192547546-1a-48,XH-Z2B654-poly-19-11-26_FKDL192547546-1a-39,XH-Z2B654-poly-19-11-27_FKDL192547546-1a-47,XH-Z2B654-pro-19-11-26_FKDL192547546-1a-37,XH-Z2B654-pro-19-11-27_FKDL192547546-1a-45,XH-Z2wt-Baso-19-11-26_FKDL192547546-1a-34,XH-Z2wt-Baso-19-11-27_FKDL192547546-1a-42,XH-Z2wt-ortho-19-11-26_FKDL192547546-1a-36,XH-Z2wt-ortho-19-11-27_FKDL192547546-1a-44,XH-Z2wt-poly-19-11-26_FKDL192547546-1a-35,XH-Z2wt-poly-19-11-27_FKDL192547546-1a-43,XH-Z2wt-pro-19-11-26_FKDL192547546-1a-33,XH-Z2wt-pro-19-11-27_FKDL192547546-1a-41};
do cutadapt --pair-filter=any --minimum-length 15 --max-n 8 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o ${i}_rmadp_1.fastq.gz -p ${i}_rmadp_2.fastq.gz ${i}_1.fq.gz ${i}_2.fq.gz >>filter.txt 2>&1
java -jar /data/yudonglin/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 20 -phred33 ${i}_rmadp_1.fastq.gz ${i}_rmadp_2.fastq.gz -baseout ${i}_fliter.fq.gz  AVGQUAL:20 SLIDINGWINDOW:4:15 MINLEN:15 1>>filter.txt 2>&1;
rm ${i}_rmadp_1.fastq.gz;
rm ${i}_rmadp_2.fastq.gz;
done


#mkdir fastqc1
#去接头之后再进行一次质量控制
#fastqc -t 30 -o ./fastqc1/ *.gz


#比对
for i in {XH-Z2B654-Baso-19-11-26_FKDL192547546-1a-38,XH-Z2B654-Baso-19-11-27_FKDL192547546-1a-46,XH-Z2B654-ortho-19-11-26_FKDL192547546-1a-40,XH-Z2B654-ortho-19-11-27_FKDL192547546-1a-48,XH-Z2B654-poly-19-11-26_FKDL192547546-1a-39,XH-Z2B654-poly-19-11-27_FKDL192547546-1a-47,XH-Z2B654-pro-19-11-26_FKDL192547546-1a-37,XH-Z2B654-pro-19-11-27_FKDL192547546-1a-45,XH-Z2wt-Baso-19-11-26_FKDL192547546-1a-34,XH-Z2wt-Baso-19-11-27_FKDL192547546-1a-42,XH-Z2wt-ortho-19-11-26_FKDL192547546-1a-36,XH-Z2wt-ortho-19-11-27_FKDL192547546-1a-44,XH-Z2wt-poly-19-11-26_FKDL192547546-1a-35,XH-Z2wt-poly-19-11-27_FKDL192547546-1a-43,XH-Z2wt-pro-19-11-26_FKDL192547546-1a-33,XH-Z2wt-pro-19-11-27_FKDL192547546-1a-41};
do hisat2 --threads 35 -x /data/yudonglin/reference/mm10/genome -1 ${i}_fliter_1P.fq.gz -2 ${i}_fliter_2P.fq.gz -S ${i}.sam;
rm ${i}_fliter_1P.fq.gz;
rm ${i}_fliter_2P.fq.gz ;
done


#sam文件转换为bam文件
for i in {XH-Z2B654-Baso-19-11-26_FKDL192547546-1a-38,XH-Z2B654-Baso-19-11-27_FKDL192547546-1a-46,XH-Z2B654-ortho-19-11-26_FKDL192547546-1a-40,XH-Z2B654-ortho-19-11-27_FKDL192547546-1a-48,XH-Z2B654-poly-19-11-26_FKDL192547546-1a-39,XH-Z2B654-poly-19-11-27_FKDL192547546-1a-47,XH-Z2B654-pro-19-11-26_FKDL192547546-1a-37,XH-Z2B654-pro-19-11-27_FKDL192547546-1a-45,XH-Z2wt-Baso-19-11-26_FKDL192547546-1a-34,XH-Z2wt-Baso-19-11-27_FKDL192547546-1a-42,XH-Z2wt-ortho-19-11-26_FKDL192547546-1a-36,XH-Z2wt-ortho-19-11-27_FKDL192547546-1a-44,XH-Z2wt-poly-19-11-26_FKDL192547546-1a-35,XH-Z2wt-poly-19-11-27_FKDL192547546-1a-43,XH-Z2wt-pro-19-11-26_FKDL192547546-1a-33,XH-Z2wt-pro-19-11-27_FKDL192547546-1a-41};
do samtools view -S ${i}.sam -b > ${i}.bam;
rm ${i}.sam;
samtools sort ${i}.bam -o ${i}_sorted.bam;
samtools index ${i}_sorted.bam;
done

conda activate pyscenic
#将bam文件转换为bw文件
#bam文件转bw文件
for i in {XH-Z2B654-Baso-19-11-26_FKDL192547546-1a-38,XH-Z2B654-Baso-19-11-27_FKDL192547546-1a-46,XH-Z2B654-ortho-19-11-26_FKDL192547546-1a-40,XH-Z2B654-ortho-19-11-27_FKDL192547546-1a-48,XH-Z2B654-poly-19-11-26_FKDL192547546-1a-39,XH-Z2B654-poly-19-11-27_FKDL192547546-1a-47,XH-Z2B654-pro-19-11-26_FKDL192547546-1a-37,XH-Z2B654-pro-19-11-27_FKDL192547546-1a-45,XH-Z2wt-Baso-19-11-26_FKDL192547546-1a-34,XH-Z2wt-Baso-19-11-27_FKDL192547546-1a-42,XH-Z2wt-ortho-19-11-26_FKDL192547546-1a-36,XH-Z2wt-ortho-19-11-27_FKDL192547546-1a-44,XH-Z2wt-poly-19-11-26_FKDL192547546-1a-35,XH-Z2wt-poly-19-11-27_FKDL192547546-1a-43,XH-Z2wt-pro-19-11-26_FKDL192547546-1a-33,XH-Z2wt-pro-19-11-27_FKDL192547546-1a-41};
do bamCoverage --bam ${i}_sorted.bam -o ${i}_sorted.bam.bw  --binSize 10 ;
rm ${i}.bam;
done

conda deactivate
#基因表达定量
for i in {XH-Z2B654-Baso-19-11-26_FKDL192547546-1a-38,XH-Z2B654-Baso-19-11-27_FKDL192547546-1a-46,XH-Z2B654-ortho-19-11-26_FKDL192547546-1a-40,XH-Z2B654-ortho-19-11-27_FKDL192547546-1a-48,XH-Z2B654-poly-19-11-26_FKDL192547546-1a-39,XH-Z2B654-poly-19-11-27_FKDL192547546-1a-47,XH-Z2B654-pro-19-11-26_FKDL192547546-1a-37,XH-Z2B654-pro-19-11-27_FKDL192547546-1a-45,XH-Z2wt-Baso-19-11-26_FKDL192547546-1a-34,XH-Z2wt-Baso-19-11-27_FKDL192547546-1a-42,XH-Z2wt-ortho-19-11-26_FKDL192547546-1a-36,XH-Z2wt-ortho-19-11-27_FKDL192547546-1a-44,XH-Z2wt-poly-19-11-26_FKDL192547546-1a-35,XH-Z2wt-poly-19-11-27_FKDL192547546-1a-43,XH-Z2wt-pro-19-11-26_FKDL192547546-1a-33,XH-Z2wt-pro-19-11-27_FKDL192547546-1a-41};
do featureCounts -T 30 -t exon -g gene_id -a /data/yudonglin/reference/gencode.vM16.basic.annotation.gtf -o ${i}.count ${i}_sorted.bam >>~/count.txt 2>&1
rm ${i}_sorted.bam;
done




#第二批结果，四个时期分别各一个
for i in {wt-pro-19-12-24_FKDL202555457-1a-A1,wt-Baso-19-12-24_FKDL202555457-1a-A2,wt-poly-19-12-24_FKDL202555457-1a-A3,wt-ortho-19-12-24_FKDL202555457-1a-A4};
do cutadapt --pair-filter=any --minimum-length 15 --max-n 8 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o ${i}_rmadp_1.fastq.gz -p ${i}_rmadp_2.fastq.gz ${i}_1.fq.gz ${i}_2.fq.gz >>filter.txt 2>&1
java -jar /data/yudonglin/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 20 -phred33 ${i}_rmadp_1.fastq.gz ${i}_rmadp_2.fastq.gz -baseout ${i}_fliter.fq.gz  AVGQUAL:20 SLIDINGWINDOW:4:15 MINLEN:15 1>>filter.txt 2>&1;
done


#mkdir fastqc1
#去接头之后再进行一次质量控制
#fastqc -t 30 -o ./fastqc1/ *.gz


#比对
for i in {wt-pro-19-12-24_FKDL202555457-1a-A1,wt-Baso-19-12-24_FKDL202555457-1a-A2,wt-poly-19-12-24_FKDL202555457-1a-A3,wt-ortho-19-12-24_FKDL202555457-1a-A4};
do hisat2 --threads 35 -x /data/yudonglin/reference/mm10/genome -1 ${i}_fliter_1P.fq.gz -2 ${i}_fliter_2P.fq.gz -S ${i}.sam;
rm ${i}_rmadp_1.fastq.gz;
rm ${i}_rmadp_2.fastq.gz;
rm ${i}_fliter_1P.fq.gz;
rm ${i}_fliter_2P.fq.gz ;
done



#sam文件转换为bam文件
for i in {wt-pro-19-12-24_FKDL202555457-1a-A1,wt-Baso-19-12-24_FKDL202555457-1a-A2,wt-poly-19-12-24_FKDL202555457-1a-A3,wt-ortho-19-12-24_FKDL202555457-1a-A4};
do samtools view -S ${i}.sam -b > ${i}.bam
rm ${i}.sam
samtools sort ${i}.bam -o ${i}_sorted.bam
samtools index ${i}_sorted.bam;
rm ${i}.bam;
done

conda activate pyscenic
#将bam文件转换为bw文件
#bam文件转bw文件
for i in {wt-pro-19-12-24_FKDL202555457-1a-A1,wt-Baso-19-12-24_FKDL202555457-1a-A2,wt-poly-19-12-24_FKDL202555457-1a-A3,wt-ortho-19-12-24_FKDL202555457-1a-A4};
do bamCoverage --bam ${i}_sorted.bam -o ${i}_sorted.bam.bw  --binSize 10 
done

conda deactivate
#基因表达定量
for i in {wt-pro-19-12-24_FKDL202555457-1a-A1,wt-Baso-19-12-24_FKDL202555457-1a-A2,wt-poly-19-12-24_FKDL202555457-1a-A3,wt-ortho-19-12-24_FKDL202555457-1a-A4};
do featureCounts -T 30 -t exon -g gene_id -a /data/yudonglin/reference/gencode.vM16.basic.annotation.gtf -o ${i}.count ${i}_sorted.bam >>~/count.txt 2>&1
rm ${i}_sorted.bam;
done












#师兄的两个重复

for i in {XH-PRO-1_HL3NLCCXY_L3,XH-PRO-2_HL3NLCCXY_L3,XH-BASO-1_HL3NLCCXY_L3,XH-BASO-2_HL3NLCCXY_L3,XH-POLY-1_HL3NLCCXY_L3,XH-POLY-2_HL3NLCCXY_L3,XH-ORTHO-1_HL3NLCCXY_L3,XH-ORTHO-2_HL3NLCCXY_L3};
do cutadapt --pair-filter=any --minimum-length 15 --max-n 8 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o ${i}_rmadp_1.fastq.gz -p ${i}_rmadp_2.fastq.gz ${i}_1.fq.gz ${i}_2.fq.gz >>filter.txt 2>&1
java -jar /data/yudonglin/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 20 -phred33 ${i}_rmadp_1.fastq.gz ${i}_rmadp_2.fastq.gz -baseout ${i}_fliter.fq.gz  AVGQUAL:20 SLIDINGWINDOW:4:15 MINLEN:15 1>>filter.txt 2>&1;
done


#mkdir fastqc1
#去接头之后再进行一次质量控制
#fastqc -t 30 -o ./fastqc1/ *.gz


#比对
for i in {XH-PRO-1_HL3NLCCXY_L3,XH-PRO-2_HL3NLCCXY_L3,XH-BASO-1_HL3NLCCXY_L3,XH-BASO-2_HL3NLCCXY_L3,XH-POLY-1_HL3NLCCXY_L3,XH-POLY-2_HL3NLCCXY_L3,XH-ORTHO-1_HL3NLCCXY_L3,XH-ORTHO-2_HL3NLCCXY_L3};
do hisat2 --threads 35 -x /data/yudonglin/reference/mm10/genome -1 ${i}_fliter_1P.fq.gz -2 ${i}_fliter_2P.fq.gz -S ${i}.sam;
rm ${i}_rmadp_1.fastq.gz;
rm ${i}_rmadp_2.fastq.gz;
rm ${i}_fliter_1P.fq.gz;
rm ${i}_fliter_2P.fq.gz ;
done


#sam文件转换为bam文件
for i in {XH-PRO-1_HL3NLCCXY_L3,XH-PRO-2_HL3NLCCXY_L3,XH-BASO-1_HL3NLCCXY_L3,XH-BASO-2_HL3NLCCXY_L3,XH-POLY-1_HL3NLCCXY_L3,XH-POLY-2_HL3NLCCXY_L3,XH-ORTHO-1_HL3NLCCXY_L3,XH-ORTHO-2_HL3NLCCXY_L3};
do samtools view -S ${i}.sam -b > ${i}.bam;
rm ${i}.sam;
samtools sort ${i}.bam -o ${i}_sorted.bam;
samtools index ${i}_sorted.bam;
rm ${i}.bam;
done

conda activate pyscenic
#将bam文件转换为bw文件
#bam文件转bw文件
for i in {XH-PRO-1_HL3NLCCXY_L3,XH-PRO-2_HL3NLCCXY_L3,XH-BASO-1_HL3NLCCXY_L3,XH-BASO-2_HL3NLCCXY_L3,XH-POLY-1_HL3NLCCXY_L3,XH-POLY-2_HL3NLCCXY_L3,XH-ORTHO-1_HL3NLCCXY_L3,XH-ORTHO-2_HL3NLCCXY_L3};
do bamCoverage --bam ${i}_sorted.bam -o ${i}_sorted.bam.bw  --binSize 10 
done


conda deactivate
#基因表达定量
for i in {XH-PRO-1_HL3NLCCXY_L3,XH-PRO-2_HL3NLCCXY_L3,XH-BASO-1_HL3NLCCXY_L3,XH-BASO-2_HL3NLCCXY_L3,XH-POLY-1_HL3NLCCXY_L3,XH-POLY-2_HL3NLCCXY_L3,XH-ORTHO-1_HL3NLCCXY_L3,XH-ORTHO-2_HL3NLCCXY_L3};
do featureCounts -T 30 -t exon -g gene_id -a /data/yudonglin/reference/gencode.vM16.basic.annotation.gtf -o ${i}.count ${i}_sorted.bam >>~/count.txt 2>&1
rm ${i}_sorted.bam;
done


