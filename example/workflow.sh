SAM1=nramp1
SAM2=nramp1_Mn
SAM3=mdg10nramp1
SAM4=mdg10nramp1_Mn
SAM1_1=nramp1_1
SAM1_2=nramp1_2
SAM1_3=nramp1_3
SAM2_1=nramp1_Mn_1
SAM2_2=nramp1_Mn_2
SAM2_3=nramp1_Mn_3
SAM3_1=mdg10nramp1_1
SAM3_2=mdg10nramp1_2
SAM3_3=mdg10nramp1_3
SAM4_1=mdg10nramp1_Mn_1
SAM4_2=mdg10nramp1_Mn_2
SAM4_3=mdg10nramp1_Mn_3


BAM_DIR=~/Project/CFH0003_20180828/RNA-Seq/bam

## Arabidopsis
GTF=/cluster/home/xfu/Gmatic7/gene/tair10/tair10.gtf
FASTA=/cluster/home/xfu/Gmatic7/genome/tair10/tair10.fa
AS=/cluster/home/xfu/Gmatic7/gene/tair10/tair10_sorted.gtf_astalavista.gtf
ANNO=/cluster/home/xfu/Gmatic7/gene/tair10/tair10_gene_anno.tsv

mkdir spanki_out
## Junction alignment evaluation
find $BAM_DIR/*.bam -printf "%f\n" | parallel --gnu "spankijunc -i $BAM_DIR/{} -o spanki_out/{}.out -g $GTF -f $FASTA"

## Generating a curated cross-sample junction set
cat spanki_out/*.bam.out/juncs.all|grep -v '^juncid'|cut -f1-7|sort|uniq > spanki_out/all_juncs.txt
cat spanki_out/*.bam.out/juncs.all|head -1|cut -f1-7 > spanki_out/all_juncs_an.txt
cat spanki_out/all_juncs.txt | awk '{if($4 == "an" || $4 == "annostatus")print}' >> spanki_out/all_juncs_an.txt
find $BAM_DIR/*.bam -printf "%f\n" | parallel --gnu "make_curated_jtab -o spanki_out/{}.out_curated -i $BAM_DIR/{} -jlist spanki_out/all_juncs_an.txt -jtab spanki_out/{}.out/juncs.all"

## Merging junction tables of replicates
mkdir spanki_out/reps_merged
merge_jtabs spanki_out/${SAM1_1}.bam.out_curated/juncs.all,spanki_out/${SAM1_2}.bam.out_curated/juncs.all,spanki_out/${SAM1_3}.bam.out_curated/juncs.all > spanki_out/reps_merged/${SAM1}.juncs
merge_jtabs spanki_out/${SAM2_1}.bam.out_curated/juncs.all,spanki_out/${SAM2_2}.bam.out_curated/juncs.all,spanki_out/${SAM2_3}.bam.out_curated/juncs.all > spanki_out/reps_merged/${SAM2}.juncs
merge_jtabs spanki_out/${SAM3_1}.bam.out_curated/juncs.all,spanki_out/${SAM3_2}.bam.out_curated/juncs.all,spanki_out/${SAM3_3}.bam.out_curated/juncs.all > spanki_out/reps_merged/${SAM3}.juncs
merge_jtabs spanki_out/${SAM4_1}.bam.out_curated/juncs.all,spanki_out/${SAM4_2}.bam.out_curated/juncs.all,spanki_out/${SAM4_3}.bam.out_curated/juncs.all > spanki_out/reps_merged/${SAM4}.juncs

## Comparison of junctions between samples
mkdir spanki_out/junccomp
junccomp -a spanki_out/reps_merged/${SAM1}.juncs -b spanki_out/reps_merged/${SAM2}.juncs -o spanki_out/junccomp/${SAM1}_vs_${SAM2}
junccomp -a spanki_out/reps_merged/${SAM3}.juncs -b spanki_out/reps_merged/${SAM4}.juncs -o spanki_out/junccomp/${SAM3}_vs_${SAM4}
junccomp -a spanki_out/reps_merged/${SAM3}.juncs -b spanki_out/reps_merged/${SAM1}.juncs -o spanki_out/junccomp/${SAM3}_vs_${SAM1}
junccomp -a spanki_out/reps_merged/${SAM4}.juncs -b spanki_out/reps_merged/${SAM2}.juncs -o spanki_out/junccomp/${SAM4}_vs_${SAM2}

## Quantification of inclusion / exclusion paths
find spanki_out/reps_merged/*.juncs| parallel --gnu "spankisplice -j {} -o {}_events -f $FASTA -g $GTF -a $AS"

## Comparison of splicing events between samples
mkdir spanki_out/splicecomp
splicecomp -a spanki_out/reps_merged/${SAM1}.juncs_events/events.out -b spanki_out/reps_merged/${SAM2}.juncs_events/events.out -o spanki_out/splicecomp/${SAM1}_vs_${SAM2}
splicecomp -a spanki_out/reps_merged/${SAM3}.juncs_events/events.out -b spanki_out/reps_merged/${SAM4}.juncs_events/events.out -o spanki_out/splicecomp/${SAM3}_vs_${SAM4}
splicecomp -a spanki_out/reps_merged/${SAM3}.juncs_events/events.out -b spanki_out/reps_merged/${SAM1}.juncs_events/events.out -o spanki_out/splicecomp/${SAM3}_vs_${SAM1}
splicecomp -a spanki_out/reps_merged/${SAM4}.juncs_events/events.out -b spanki_out/reps_merged/${SAM2}.juncs_events/events.out -o spanki_out/splicecomp/${SAM4}_vs_${SAM2}

## combine and filter AS events: p-value < 0.05
mkdir tables
cat spanki_out/splicecomp/${SAM1}_vs_${SAM2}/event_compare.out |grep -v 'Unclassified'|cut -f2,4,5,7,23,24,25|awk '{print $4"\t"$3"\t"$1"\t"$2"\t"$5"\t"$6"\t"$7}' > tables/${SAM1}_vs_${SAM2}_event_compare.out
cat spanki_out/splicecomp/${SAM3}_vs_${SAM4}/event_compare.out |grep -v 'Unclassified'|cut -f2,4,5,7,23,24,25|awk '{print $4"\t"$3"\t"$1"\t"$2"\t"$5"\t"$6"\t"$7}' > tables/${SAM3}_vs_${SAM4}_event_compare.out
cat spanki_out/splicecomp/${SAM3}_vs_${SAM1}/event_compare.out |grep -v 'Unclassified'|cut -f2,4,5,7,23,24,25|awk '{print $4"\t"$3"\t"$1"\t"$2"\t"$5"\t"$6"\t"$7}' > tables/${SAM3}_vs_${SAM1}_event_compare.out
cat spanki_out/splicecomp/${SAM4}_vs_${SAM2}/event_compare.out |grep -v 'Unclassified'|cut -f2,4,5,7,23,24,25|awk '{print $4"\t"$3"\t"$1"\t"$2"\t"$5"\t"$6"\t"$7}' > tables/${SAM4}_vs_${SAM2}_event_compare.out

Rscript script/merge_AS_out.R ${SAM1}_vs_${SAM2} ${SAM3}_vs_${SAM4} ${SAM3}_vs_${SAM1} ${SAM4}_vs_${SAM2}

## gene annotation
Rscript script/AS_gene_anno.R $ANNO tables/merge_AS_out
