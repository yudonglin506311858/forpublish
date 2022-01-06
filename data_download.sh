mkdir cellreport
cd cellreport/
#第一篇cell report

#SRX5516957: GSM3669548: ProE_ATAC_1; Homo sapiens; ATAC-seq
#SRX5516958: GSM3669549: ProE_ATAC_2; Homo sapiens; ATAC-seq
#SRX5516959: GSM3669550: EBaso_ATAC_1; Homo sapiens; ATAC-seq
#SRX5516960: GSM3669551: EBaso_ATAC_2; Homo sapiens; ATAC-seq
#SRX5516961: GSM3669552: LBaso_ATAC_1; Homo sapiens; ATAC-seq
#SRX5516962: GSM3669553: LBaso_ATAC_2; Homo sapiens; ATAC-seq
#SRX5516963: GSM3669554: Poly_ATAC_1; Homo sapiens; ATAC-seq
#SRX5516964: GSM3669555: Poly_ATAC_2; Homo sapiens; ATAC-seq
#SRX5516965: GSM3669556: Ortho_ATAC_1; Homo sapiens; ATAC-seq
#SRX5516966: GSM3669557: Ortho_ATAC_2; Homo sapiens; ATAC-seq
#ATAC-seq
for((i=0;i<=9;i++));
do /data/yudonglin/software/ena-fast-download/ena-fast-download.py SRR872373${i};
done

for((i=4;i<=9;i++));
do /data/yudonglin/software/ena-fast-download/ena-fast-download.py SRR872372${i};
done

#RNA-seq
#HSC-BFUE-CFUE
for((i=65;i<=73;i++));
do /data/yudonglin/software/ena-fast-download/ena-fast-download.py SRR87237${i};
done
#安秀丽
#SRX424826: GSM1304777: hs_proerythroblast_1; Homo sapiens; RNA-Seq SRR1106084
#SRX424852: GSM1304803: mm_orthochromatic_3; Mus musculus; RNA-Seq SRR1106110

#for((i=6084;i<=6110;i++));
#do /data/yudonglin/software/ena-fast-download/ena-fast-download.py SRR110${i};
#done

for((i=84;i<=99;i++));
do fastq-dump --gzip SRR11060${i};
done

for((i=0;i<=9;i++));
do fastq-dump --gzip SRR110610${i};
done
fastq-dump --gzip SRR1106110




#第二篇cell report
#ATAC-seq
for((i=61;i<=88;i++));
do /data/yudonglin/software/ena-fast-download/ena-fast-download.py SRR72952${i};
done

#RNA-seq
for((i=71;i<=98;i++));
do /data/yudonglin/software/ena-fast-download/ena-fast-download.py SRR72954${i};
done
