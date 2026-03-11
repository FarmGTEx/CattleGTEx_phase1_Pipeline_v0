#!/usr/bin/env python3

"""
QTLEnrichV2 was written under the direction and supervision of Ayellet Segre
at Massachusetts Eye and Ear, Department of Ophthalmology and Ocular Genomics Institute

Author: Andrew Hamel

Code written in Python version 3.6

QTLEnrichV2 assesses enrichment of gwas variant associations amongst qtls in
a given tissue correcting for potential confounding factors 
(MAF, distance to TSS, local LD).
"""

import sys
import os
import os.path
import re
import argparse
import logging
import datetime
import pandas as pd
import numpy as np
import os_sys_qtlenrich
import qtl_gwas_parsing 
import rand_null_samp_permute_matching_confounders
import compute_fold_enrichments
import output_table
import GeneEnrich_inputs
from argparse import Namespace

if __name__ == "__main__":
    args = os_sys_qtlenrich.parse_args()
    #args = Namespace(**args)
    required_args = ['exp_label', 'qtl_type', 'gwas_file', 'trait_name', 'qtl_directory', 'file_name',
                     'confounders_table', 'null_table', 'qtl_q_value', 'gwas_p_value', 'num_quantiles']
    for arg in required_args:
        if not hasattr(args, arg):
            raise ValueError(f"Missing required argument: {arg}")

    exp_label = os_sys_qtlenrich.parse_exp_label(args.exp_label)

    gwas_trait = qtl_gwas_parsing.extract_gwas_trait_name(args.gwas_file, args.trait_name)

    output_directory = os_sys_qtlenrich.create_output_directory(args, datetime.date.today(), exp_label)
    output_files_directory = os_sys_qtlenrich.create_output_files_directory(output_directory)

    qtl_files = qtl_gwas_parsing.extract_qtl_files(args.qtl_directory, args.file_name)

    if hasattr(args, "subset_tissues") and args.subset_tissues:
        qtl_files = os_sys_qtlenrich.subset_qtl_files(qtl_files, args.file_name, args.tissue_names)

    print("Preparing GWAS data...")
    gwas = qtl_gwas_parsing.prepare_gwas_file(args.gwas_file, "ARS-UCD1.2")

    confounders_table = qtl_gwas_parsing.read_confounders_table(args.confounders_table)
    null_table = qtl_gwas_parsing.prepare_null_table(args.null_table, gwas)

    #gencode = qtl_gwas_parsing.declare_gencode(args)

    print("Parsing QTLs...")
    significant_qtl_dict, significant_qtl_gwas_dict, geneenrich_input_gwas_dict = qtl_gwas_parsing.parse_qtls(
        output_files_directory, datetime.date.today(), gwas_trait, qtl_files, args.qtl_directory,
        args.qtl_type, args.file_name, gwas, 0.05, False, False, None, False, False
    )

    if getattr(args, "GeneEnrich_input", False):
        print("Generating GeneEnrich input files...")
        if args.qtl_type in ["dapg_independent", "conditional_independent"]:
            os_sys_qtlenrich.check_independent_GeneEnrich_parameters(args)
            GeneEnrich_inputs.create_GeneEnrich_inputs_independent_qtls(
                args.eGenes_directory, args.eGene_file_name, gwas, gwas_trait, output_directory,
                datetime.date.today(), geneenrich_input_gwas_dict, args.qtl_type, p=args.gwas_p_value,
                q=args.qtl_q_value, independent_ranking=args.independent_ranking, gencode=gencode,
                subset_genes=args.subset_genes
            )
        else:
            GeneEnrich_inputs.create_GeneEnrich_inputs_best_eqtl_sqtl(
                gwas_trait, output_directory, datetime.date.today(), geneenrich_input_gwas_dict,
                args.qtl_type, exp_label, p=args.gwas_p_value, q=0.05
            )

    del geneenrich_input_gwas_dict, gwas

    print("Computing observed fold-enrichment...")
    trait_dict, original_length_dict, observed_fold, obs_dict, best_eqtl_variants_dict = compute_fold_enrichments.compute_qtl_statistics(
        gwas_trait, significant_qtl_dict, significant_qtl_gwas_dict, args.gwas_p_value
    )
    log_file="/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/QTLenrich/ct_seQTL/result/QTLEnrich_eQTL_ct.Dtr_Calv_Ease_ct.Dtr_Calv_Ease_best_eqtl_Feb20_2025_16:40:22.log"
    enrichment_p_value_dict, adjusted_fold_enrichment_dict, upper_bound_confidence_interval_dict, lower_bound_confidence_interval_dict, pi_1_dict, true_trait_dict = rand_null_samp_permute_matching_confounders.sample_match_null_confounders(
        log_file, datetime.date.today(), output_files_directory, 'best_eqtl', gwas_trait,
        confounders_table, null_table, significant_qtl_gwas_dict, observed_fold, 0.05
        #num_permutations=args.lower_bound_permutations,
        #upper_bound_permutations=args.upper_bound_permutations, lambda_factor=args.lambda_factor,
        #independent_ranking=args.independent_ranking, compute_tss_distance=args.compute_tss_distance,
        #keep_null_variants=args.keep_null_variants, keep_pvalue_matrix=args.keep_pvalue_matrix,
        #interaction_qtl_dict={}, num_quantiles=args.num_quantiles
    )

    print("Generating output table...")
    results_directory = os_sys_qtlenrich.create_results_directory(output_directory)
    output_filename = os_sys_qtlenrich.create_output_filename(
        results_directory, '2025', gwas_trait, 'best_eqtl', None, exp_label
    )
    
    def estimated_true_trait(output):
        output["Estimated_num_trait_associations_with_GWAS_p<0.05"] = (
            (output["Adjusted_Fold_Enrichment"] - 1) / output["Adjusted_Fold_Enrichment"]
        ) * output["Num_OBS_QTLs_with_GWAS_p<0.05"]
        output["Estimated_num_trait_associations_with_GWAS_p<0.05"].replace([np.inf, -np.inf], np.nan, inplace=True)
        output["Estimated_num_trait_associations_with_GWAS_p<0.05"].fillna(0, inplace=True)  
        output["Estimated_num_trait_associations_with_GWAS_p<0.05"] = np.nan_to_num(
          output["Estimated_num_trait_associations_with_GWAS_p<0.05"], nan=0, posinf=0, neginf=0
        )
        output["Estimated_num_trait_associations_with_GWAS_p<0.05"] = output["Estimated_num_trait_associations_with_GWAS_p<0.05"].astype(int)
        return(output)

    def expected_trait(output):
        """
        Compute expected num trait associations based on adjusted fold-enrichment

        equation: OBS * (1-((adjusted fold-enrichment - 1)/adjusted fold-enrichment))
        """
        output = output[output["Adjusted_Fold_Enrichment"] > 1e-6]
        output["Num_EXP_QTLs_with_GWAS_p<0.05"] = (output["Num_OBS_QTLs_with_GWAS_p<0.05"])*(1-((output["Adjusted_Fold_Enrichment"]-1)/output["Adjusted_Fold_Enrichment"]))
        output["Num_EXP_QTLs_with_GWAS_p<0.05"] = output["Num_EXP_QTLs_with_GWAS_p<0.05"].astype(int)
        return(output)

    def dict_to_df(dict,col_name):
        """
        Converts dictionary into a dataframe
        """
        return(pd.DataFrame(dict.items(),columns = ["Tissue_(QTL)",col_name]))

    def trait_dict_to_df(trait_dict):
        """
        Special case of converting dictionary to dataframe
        """
        return(pd.DataFrame(trait_dict).melt().rename(columns = {"variable":"Trait_(GWAS)","value":"Tissue_(QTL)"}))

    def sort_output_table(output):
        """
        Take QTLEnrich output table sort by enrichment p-value and adjusted fold-enrichment.

        Enrichment p-value is sorted in ascending order.
        Adjusted fold-enrichment is sorted in descending order.
        """
        output["intermediate_columns"] = output["Enrichment_P_value"].apply(lambda x: np.float32(x))
        output = output.sort_values(["intermediate_columns","Adjusted_Fold_Enrichment"],ascending=[True,False])
        output = output.drop(["intermediate_columns"],axis=1)
        return(output)

    def merge_tables(df_list):
        """
        Merges list of dataframes into one dataframe
        """
        output = df_list[0]
        for df_ in df_list[1:]:
            output = output.merge(df_, on="Tissue_(QTL)")
        output["Tissue_(QTL)"] = output["Tissue_(QTL)"].str.replace("Tissue_","")
        return(output)
    
    trait_dict_df = trait_dict_to_df(trait_dict)
    original_length_dict_df = dict_to_df(original_length_dict,col_name="Num_QTL_variants")
    best_QTL_variants_dict_df = dict_to_df(best_eqtl_variants_dict,col_name="Num_QTL_variants_in_GWAS")
    obs_dict_df = dict_to_df(obs_dict,col_name="Num_OBS_QTLs_with_GWAS_p<0.05")
    observed_fold_enrichment_df = dict_to_df(observed_fold,col_name="Fold-Enrichment_(OBS/EXP)")
    adjusted_fold_enrichment_dict_df = dict_to_df(adjusted_fold_enrichment_dict,col_name="Adjusted_Fold_Enrichment")
    lower_bound_confidence_interval_dict_df = dict_to_df(lower_bound_confidence_interval_dict,col_name="Lower_bound_95%_CI")
    upper_bound_confidence_interval_dict_df = dict_to_df(upper_bound_confidence_interval_dict,col_name="Upper_bound_95%_CI")
    enrichment_p_value_dict_df = dict_to_df(enrichment_p_value_dict,col_name="Enrichment_P_value")
    pi1_dict_df = dict_to_df(pi_1_dict,col_name="Empirical_Pi1")
    true_trait_dict_df = dict_to_df(true_trait_dict,col_name="Estimated_total_num_trait_associations")
    df_list = [trait_dict_df,original_length_dict_df,best_QTL_variants_dict_df,
               obs_dict_df,observed_fold_enrichment_df,adjusted_fold_enrichment_dict_df,
               lower_bound_confidence_interval_dict_df,upper_bound_confidence_interval_dict_df,enrichment_p_value_dict_df,pi1_dict_df,true_trait_dict_df]
    output_table = merge_tables(df_list)
    output_table = estimated_true_trait(output_table)
    output_table = expected_trait(output_table)
    output_table = output_table[["Trait_(GWAS)","Tissue_(QTL)","Num_QTL_variants","Num_QTL_variants_in_GWAS",
                                 "Num_OBS_QTLs_with_GWAS_p<0.05","Num_EXP_QTLs_with_GWAS_p<0.05","Estimated_num_trait_associations_with_GWAS_p<0.05",
                                 "Empirical_Pi1","Estimated_total_num_trait_associations","Fold-Enrichment_(OBS/EXP)",
                                 "Adjusted_Fold_Enrichment","Lower_bound_95%_CI","Upper_bound_95%_CI","Enrichment_P_value"]].copy()
    output_table = sort_output_table(output_table)
    #output_df = output_table.create_output_table(
        #trait_dict, original_length_dict, best_eqtl_variants_dict, observed_fold, obs_dict,
        #enrichment_p_value_dict, adjusted_fold_enrichment_dict, upper_bound_confidence_interval_dict,
        #lower_bound_confidence_interval_dict, pi_1_dict, true_trait_dict
    #)

    output_table.to_csv(output_filename, index=None, sep="\t")
    print("QTLEnrich analysis complete. Results saved to:", output_filename)
