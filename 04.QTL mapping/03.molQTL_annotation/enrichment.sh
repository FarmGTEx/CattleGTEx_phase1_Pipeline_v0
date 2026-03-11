# to do enrichment analysis using QTL loci
QTL_design_file="/faststorage/project/farmgtex/pipeline/QTL_design"

# collect all the QTL loci in all the tissues
rm -f ./all_eQTL ./temp_file
cat $QTL_design_file | cut -f 1 |sort | uniq | sed "s/ /_/g"|while read tissue
  do cat /faststorage/project/farmgtex/pipeline/QTL/${tissue}/${tissue}.cis_qtl_pairs.significant.txt | cut -f 2  >> temp_file
done
  cat temp_file |grep -v "id" |sort -n | uniq > all_eQTL
  
  cat all_eQTL | sed "s/_/\t/g" > keep.txt
  bcftools view -s 1_NL_2_B_bqsr /faststorage/project/farmgtex/gtex/VCFfile/chr_final/recal/final_vcf/VCF_file/INFO_0.75_all.vcf.gz  -Ov -o ./test.vcf 
  snpEff   eff -csvStats variants.SnpEff.csv -s variants.SnpEff.html -c snpEff.config -v Cattle test.vcf >annotated.vcf
  
   
  
  
  
  mkdir annotation_file
  mkdir QTL_loci
  bcftools query -f '%CHROM\t%POS\t%ANN\n' annotated.vcf > annotation.txt
  less annotation.txt | sed "s/|/\t/g" | cut -f 1,2,4|awk '{print $1"_"$2"\t"$3}' > SNP_annotate
  less SNP_annotate | cut -f 1 >./QTL_loci/all_SNP.txt
  Rscript classify.r
  cat type_counts.txt |cut -f 1 | while read sample
  do mv ${sample}* ./annotation_file/
  done
  cp all_eQTL ./QTL_loci/all_eQTL.txt
  zcat /faststorage/project/farmgtex/pipeline/code/fine_mapping/susier.credible.sig.lead.gz |cut -f 3  | grep -v "SNP" | sort| uniq >./QTL_loci/lead_eQTL.txt
  
  
  ###tissue by tissue
  cat $QTL_design_file | cut -f 1 |sort | uniq | sed "s/ /_/g"|while read tissue
  do cat /faststorage/project/farmgtex/pipeline/QTL/${tissue}/${tissue}.cis_qtl_pairs.significant.txt | cut -f 2 | grep -v "id" >./QTL_loci/${tissue}.txt
  done
  
  
rm -f  ./temp_file
cat $QTL_design_file | cut -f 1 |sort | uniq | sed "s/ /_/g"|while read tissue
  do zcat /faststorage/project/farmgtex/QTL_result/sQTL/${tissue}/${tissue}.cis_qtl.txt.gz | cut -f 7  >> temp_file
done
  cat temp_file |grep -v "id" |sort -n | uniq > sQTL.txt
  
rm -f  ./temp_file
cat $QTL_design_file | cut -f 1 |sort | uniq | sed "s/ /_/g"|while read tissue
  do zcat /faststorage/project/farmgtex/QTL_result/isoform_QTL/${tissue}/${tissue}/${tissue}.cis_qtl.txt.gz | cut -f 6  >> temp_file
done
  cat temp_file |grep -v "id" |sort -n | uniq > isoQTL.txt
  
  
  rm -f  ./temp_file
cat $QTL_design_file | cut -f 1 |sort | uniq | sed "s/ /_/g"|while read tissue
  do zcat /faststorage/project/farmgtex/QTL_result/polyA_QTL/${tissue}/${tissue}.cis_qtl.txt.gz | cut -f 6  >> temp_file
done
  cat temp_file |grep -v "id" |sort -n | uniq > ./QTL_loci/3aQTL.txt