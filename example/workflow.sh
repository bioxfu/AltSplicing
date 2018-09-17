SAM1=nramp1_-Mn
SAM2=nramp1_+Mn
SAM3=mdg10nramp1_-Mn
SAM4=mdg10nramp1_+Mn
SAM1_1=nramp1-1_-Mn
SAM1_2=nramp1-2_-Mn
SAM1_3=nramp1-3_-Mn
SAM2_1=nramp1-1_+Mn
SAM2_2=nramp1-2_+Mn
SAM2_3=nramp1-3_+Mn
SAM3_1=mdg10nramp1-1_-Mn
SAM3_2=mdg10nramp1-2_-Mn
SAM3_3=mdg10nramp1-3_-Mn
SAM4_1=mdg10nramp1-1_+Mn
SAM4_2=mdg10nramp1-2_+Mn
SAM4_3=mdg10nramp1-3_+Mn

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

## filter AS events: p-value < 0.05
mkdir tables
cat spanki_out/splicecomp/${SAM1}_vs_${SAM2}/event_compare.out |awk '{if($1=="event_id" || $24<0.05)print}'|grep -v 'Unclassified' > tables/${SAM1}_vs_${SAM2}_event_compare.out
cat spanki_out/splicecomp/${SAM3}_vs_${SAM4}/event_compare.out |awk '{if($1=="event_id" || $24<0.05)print}'|grep -v 'Unclassified' > tables/${SAM3}_vs_${SAM4}_event_compare.out
cat spanki_out/splicecomp/${SAM3}_vs_${SAM1}/event_compare.out |awk '{if($1=="event_id" || $24<0.05)print}'|grep -v 'Unclassified' > tables/${SAM3}_vs_${SAM1}_event_compare.out
cat spanki_out/splicecomp/${SAM4}_vs_${SAM2}/event_compare.out |awk '{if($1=="event_id" || $24<0.05)print}'|grep -v 'Unclassified' > tables/${SAM4}_vs_${SAM2}_event_compare.out

## gene annotation
Rscript script/AS_gene_anno.R $ANNO tables/${SAM1}_vs_${SAM2}_event_compare.out
Rscript script/AS_gene_anno.R $ANNO tables/${SAM3}_vs_${SAM4}_event_compare.out
Rscript script/AS_gene_anno.R $ANNO tables/${SAM3}_vs_${SAM1}_event_compare.out
Rscript script/AS_gene_anno.R $ANNO tables/${SAM4}_vs_${SAM2}_event_compare.out
