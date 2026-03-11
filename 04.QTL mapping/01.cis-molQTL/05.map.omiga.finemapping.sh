#!/bin/bash

basename -a /faststorage/project/farmgtex/QTL_result_new/*/ | grep "QTL" | grep -v "_" | while read type
do
    less "/faststorage/project/farmgtex/pipeline/eQTL_tissues" | while read tissue
      do
      qtl_work_dir="/faststorage/project/farmgtex/pipeline/${type}/${tissue}/"
      qtlresult_dir="/faststorage/project/farmgtex/QTL_result_new/${type}/${tissue}/"
      gene_list_dir="/faststorage/project/farmgtex/pipeline/susieR//${type}/gene_list/${tissue}/"
      result_dir="/faststorage/project/farmgtex/pipeline/susieR//${type}/${tissue}/"
      geno_dir="/faststorage/project/farmgtex/QTL_result_new/genotype/${tissue}/"
      log_dir="${result_dir}/log/"
      mkdir -p ${gene_list_dir}
      cd ${gene_list_dir}
      zcat ${qtlresult_dir}/${tissue}.cis_qtl.txt.sig.gz | awk 'NR==1{for(i=1;i<=NF;i++){if($i=="pheno_id"){col=i;break}};print $col} NR>1{print $col}' | grep -v "pheno" |split -dl 200
    done
done



basename -a /faststorage/project/farmgtex/QTL_result_new/*/ | grep "QTL" | grep -v "_"  | while read type
do
   less "/faststorage/project/farmgtex/pipeline/eQTL_tissues" | while read tissue
  do   
    geno_dir="/faststorage/project/farmgtex/QTL_result_new/genotype/${tissue}/"
    gene_list_dir="/faststorage/project/farmgtex/pipeline/susieR//${type}/gene_list/${tissue}/"
    qtl_work_dir="/faststorage/project/farmgtex/pipeline/${type}/${tissue}/"
    qtlresult_dir="/faststorage/project/farmgtex/QTL_result_new/${type}/${tissue}/"
    result_dir="/faststorage/project/farmgtex/pipeline/susieR//${type}/"
    log_dir="${result_dir}/log/"
    #cat ${geno_dir}/${tissue}.bim | sed "s/ /\t/g" > ${geno_dir}/${tissue}.bim_tab
    ls ${gene_list_dir}/x?? | xargs -n 1 basename|while read num
      do 
        cat susier.sh | sed "s|PATH|${log_dir}|g" > temp.sh
        sbatch temp.sh ${tissue} ${num} ${qtlresult_dir} ${qtl_work_dir} ${geno_dir} ${result_dir} ${gene_list_dir}
      done
  done
done       

basename -a /faststorage/project/farmgtex/QTL_result_new/*/ | grep "QTL" | grep -v "_"  | while read type
do     
  less "/faststorage/project/farmgtex/pipeline/eQTL_tissues" | while read tissue
  do
      qtlresult_dir="/faststorage/project/farmgtex/QTL_result_new/${type}/${tissue}/"
      result_dir="/faststorage/project/farmgtex/pipeline/susieR//${type}/"
      log_dir="${result_dir}/log/"
        cat susier_combine.sh | sed "s|PATH|${log_dir}|g" >temp.sh
           sbatch temp.sh ${tissue} ${qtlresult_dir} ${result_dir} 
  done       
done


basename -a /faststorage/project/farmgtex/QTL_result_new/*/ | grep "QTL" | grep -v "_"  | while read type
do
  result_dir="/faststorage/project/farmgtex/pipeline/susieR//${type}/"
  zcat ${result_dir}/*/result/*.susier.credible.gz | awk 'NR==1||$1!="tissue"' | gzip -c > ${result_dir}/susier.credible.gz
  zcat ${result_dir}/*/result/*.susier.credible.sig.gz | awk 'NR==1||$1!="tissue"' | gzip -c > ${result_dir}/susier.credible.sig.gz
  zcat ${result_dir}/*/result/*.susier.credible.sig.lead.gz | awk 'NR==1||$1!="tissue"' | gzip -c > ${result_dir}/susier.credible.sig.lead.gz
  cat ${result_dir}/*/result/*.crediblenum.txt > ${result_dir}/susier.crediblenum.txt
  cat ${result_dir}/*/result/*.causalnum.txt > ${result_dir}/susier.causalnum.txt
done
