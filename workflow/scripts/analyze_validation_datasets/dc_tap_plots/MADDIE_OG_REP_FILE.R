#!/usr/bin/env Rscript
#for encode wtc11 tap random screeen
#Evelyn Jagoda

#in this version pooled all reps and lanes together
args = commandArgs(trailingOnly=TRUE)

#devtools::install_github("katsevich-lab/sceptre")

library(dplyr)
library(ggplot2)
library(cowplot)
library(Seurat)
library(stringr)
library(Matrix)
library(sceptre)


setwd("/oak/stanford/groups/engreitz/Users/ejagoda/230327_Encode_K562_Tap_seq_full_seq")

mtx_name = "/oak/stanford/groups/engreitz/Users/ejagoda/231114-WTC11-ENCODE-DC-TAP-MOI5-Frozen-16-Lanes/count/"
perturb_status_name = "x"
gene_gRNA_group_pairs_name = "/oak/stanford/groups/engreitz/Users/ejagoda/230901_WTC11_ENCODE_DC_TAP_MOI5_FreshvsFrozen/new_merged_wtc11_gene_targets_table.txt"
gRNA_groups_table_name = "/oak/stanford/groups/engreitz/Users/ejagoda/230901_WTC11_ENCODE_DC_TAP_MOI5_FreshvsFrozen/new_merged_wtc11_grna_groups_table.txt" #"gRNA_groups_table_ej_use.txt"
sample_name = "231114-WTC11-ENCODE-DC-TAP-MOI5-Frozen-16-Lanes_new_sceptre_merged" #"moi5_sample1_only_testing_for_his_threshold_server"
sample_number = "all" #should be either 1 or all, will do the different options based on that"
guides_pre_assigned = "no" #should be yes or no, if yes use perturb, if no use sceptre
moi = 'high'
side = "both"
grna_integration_strategy = "union"
guide_assignment_method =  "thresholding"
outdir = "/oak/stanford/groups/engreitz/Users/ejagoda/231114-WTC11-ENCODE-DC-TAP-MOI5-Frozen-16-Lanes/"
do_pos = "yes"
formula_object_mod = "x"


mtx_name = args[1]
perturb_status_name = args[2] #moi5_sample_1_get_perturbation_status_transposed_from_sh.txt fixed_moi5_test_pertubation_matrix_u.txt
gene_gRNA_group_pairs_name = args[3] #gene_gRNA_group_pairs_use.txt
gRNA_groups_table_name = args[4] #"gRNA_groups_table_ej_use.txt"
sample_name = args[5] #"moi5_sample1_only_testing_for_his_threshold_server"
sample_number = args[6] #should be either 1 or all, will do the different options based on that"
guides_pre_assigned = args[7] #should be yes or no, if yes use perturb, if no use sceptre
moi = args[8]
side = args[9]
grna_integration_strategy = args[10]
guide_assignment_method = args[11]
outdir = args[12]
do_pos = args[13]

setwd(outdir)

sample_name = paste0(sample_name,"_side_",side,"_grna_integration_",grna_integration_strategy,"_guide_assignment_",guide_assignment_method,"_formula_mod_",formula_object_mod)

mtx = Read10X(mtx_name)

gRNA_matrix <- mtx[["CRISPR Guide Capture"]]
gene_matrix <- mtx[["Gene Expression"]]

#split them up
batch_tab = data.frame(cbind(colnames(gene_matrix)))

for (i in 1:nrow(batch_tab)){
  batch = str_split(batch_tab[i,1],"-")[[1]][2] ##or - depending
  batch_tab$batch[i] = batch
}
row.names(batch_tab) = batch_tab$cbind.colnames.gene_matrix..
lanes = unique(batch_tab$batch)
num_lanes = length(lanes)
sample1 = sample(1:num_lanes, size=floor(0.5*num_lanes))
sample2 = setdiff(lanes, sample1)
print(paste0("batch1 lanes: ", sample1))
print(paste0("batch2 lanes: ", sample2))
batch1 <- batch_tab[batch_tab$batch %in% sample1, ]
batch2 <- batch_tab[batch_tab$batch %in% sample2, ]
colnames(batch1) <- colnames(batch_tab)
colnames(batch2) <- colnames(batch_tab)
batch1$cbind.colnames.gene_matrix.. = NULL
batch2$cbind.colnames.gene_matrix.. = NULL

batches <- list(batch1=batch1, batch2=batch2)
for (batch_name in names(batches)){
    batch_tab <- batches[[batch_name]]
    gRNA_matrix_sub <- gRNA_matrix[, colnames(gRNA_matrix) %in% row.names(batch_tab)]
    gene_matrix_sub <- gene_matrix[, colnames(gene_matrix) %in% row.names(batch_tab)]
    #need to add to server version
    saveRDS(gRNA_matrix_sub,paste0(outdir,sample_name,"_gRNA_matrix_", batch_name, ".RDS"))
    saveRDS(gene_matrix_sub,paste0(outdir,sample_name,"_gene_matrix_", batch_name, ".RDS"))

    gene_gRNA_group_pairs = read.table(gene_gRNA_group_pairs_name,header=T,sep = '\t')

    gRNA_groups_table_ej = read.table(gRNA_groups_table_name,header=T,sep = '\t')

    extra_covariates = batch_tab

    if (guides_pre_assigned == "yes"){
        gRNA_groups_table_ej_use = gRNA_groups_table_ej[gRNA_groups_table_ej$grna_id %in% row.names(perturbation_matrix),]
    }

    if (guides_pre_assigned == "no"){
        gRNA_groups_table_ej_use = gRNA_groups_table_ej[gRNA_groups_table_ej$grna_id %in% row.names(gRNA_matrix_sub),]
    }

    gene_gRNA_group_pairs1 = gene_gRNA_group_pairs[gene_gRNA_group_pairs$response_id %in% row.names(gene_matrix_sub),]

    response_names = row.names(gene_matrix)

    #colnames(gRNA_groups_table_ej_use)[2] = "grna_target"

    #takeout extra covariates if it's 1 sample - but then have to take out any zeros ahead of time

    sceptre_object <- import_data(
        response_matrix = gene_matrix_sub,
        grna_matrix = gRNA_matrix_sub,
        grna_target_data_frame =gRNA_groups_table_ej_use,
        moi = moi,
        extra_covariates = extra_covariates,
        response_names = response_names
    )

    #check for zeros and remove

    ###new thing to update, make positive controls tab something you premake
    #make like this


    #need to have this be more comlicated, there's a few that are like TSS_1, TSS_2
    colnames(gene_gRNA_group_pairs1)[2] = "grna_target"
    gene_gRNA_group_pairs1$type = "x"

    #for (i in 1:nrow(gene_gRNA_group_pairs1)){
    #  if (gene_gRNA_group_pairs1$target_type[i] == "TSS" & (gene_gRNA_group_pairs1$grna_target[i] == paste0(gene_gRNA_group_pairs1$response_id[i],"_TSS"))){
    #   gene_gRNA_group_pairs1$type[i] = "pos"
    #  }
    #}

    #for (i in 1:nrow(gene_gRNA_group_pairs1)){
    #  if (gene_gRNA_group_pairs1$target_type[i] == "TSS" & (grepl(gene_gRNA_group_pairs1$grna_target[i],pattern = paste0(gene_gRNA_group_pairs1$response_id[i],"_TSS")))){
    #    gene_gRNA_group_pairs1$type[i] = "pos"
    #  }
    #}

    for (i in 1:nrow(gene_gRNA_group_pairs1)){
      if (gene_gRNA_group_pairs1$target_type[i] == "pos"){
        gene_gRNA_group_pairs1$type[i] = "pos"
      }
    }

    positive_control_pairs_all =  unique(gene_gRNA_group_pairs1[gene_gRNA_group_pairs1$type == "pos",c("grna_target", "response_id")])
    discovery_pairs = unique(gene_gRNA_group_pairs1[gene_gRNA_group_pairs1$type != "pos",c("grna_target", "response_id")])
    positive_control_pairs = positive_control_pairs_all[positive_control_pairs_all$grna_target %in% gRNA_groups_table_ej_use$grna_target, ]


    sceptre_object <- set_analysis_parameters(
        sceptre_object = sceptre_object,
        discovery_pairs = unique(discovery_pairs),
        positive_control_pairs = positive_control_pairs,
        side = side,
        grna_integration_strategy = grna_integration_strategy
    )


    print(sceptre_object)
    #log(response_n_nonzero) + log(response_n_umis) + log(grna_n_nonzero + 1) + log(grna_n_umis + 1) + batch is the default? wanna put it back to what I had

    png(paste0(sample_name,"_", batch_name, "_plot_grna_count_distributions.png"))
    p = plot_grna_count_distributions(
        sceptre_object = sceptre_object,
        n_grnas_to_plot = 9)
    print(p)
    dev.off()

    #use the default
    sceptre_object_guides_assigned <- assign_grnas(
        sceptre_object = sceptre_object,
        method = guide_assignment_method, parallel = F
    )

    png(paste0(sample_name,"_", batch_name, "_plot_grna_cell_assignment_method",guide_assignment_method,".png"))
    p = plot(sceptre_object_guides_assigned,n_grnas_to_plot =5)
    print(p)
    dev.off()

    sceptre_object = sceptre_object_guides_assigned

    print(sceptre_object)

    png(paste0(sample_name,"_", batch_name, "_plot_covariates.png"))
    p = plot_covariates(sceptre_object)
    print(p)
    dev.off()

    sceptre_object = run_qc(sceptre_object = sceptre_object)

    png(paste0(sample_name,"_", batch_name, "_plot_qc.png"))
    p = plot(sceptre_object)
    print(p)
    dev.off()

    print(sceptre_object)

    sceptre_object <- run_calibration_check(
        sceptre_object = sceptre_object,
        parallel = FALSE #set to T on comp
    )

    png(paste0(sample_name,"_", batch_name, "_calibration_check.png"))
    p = plot(sceptre_object)
    print(p)
    dev.off()

    calibration_result <- get_result(
        sceptre_object = sceptre_object,
        analysis = "run_calibration_check"
    )

    write.table(calibration_result,paste0(sample_name,"_", batch_name, "_calibration_result.txt"),quote = F,sep = '\t',row.names = F)

    print(sceptre_object)

    if (do_pos == "yes"){
        sceptre_object <- run_power_check(
            sceptre_object = sceptre_object,
            parallel = FALSE
        )
        png(paste0(sample_name,"_", batch_name, "_power_check.png"))
        p = plot(sceptre_object)
        print(p)
        dev.off()
    }



    sceptre_object <- run_discovery_analysis(
        sceptre_object = sceptre_object,
        parallel = FALSE
    )

    png(paste0(sample_name,"_", batch_name, "_discovery_result.png"))
    p = plot(sceptre_object)
    print(p)
    dev.off()

    discovery_result <- get_result(
        sceptre_object = sceptre_object,
        analysis = "run_discovery_analysis"
    )
    head(discovery_result)

    write.table(discovery_result,paste0(sample_name,"_", batch_name, "_discovery_result.txt"),quote = F,sep = '\t',row.names = F)

    print(sceptre_object)

    write_outputs_to_directory(sceptre_object = sceptre_object,directory = paste0(getwd(),"/",sample_name,"_", batch_name, "secptre_output_all"))

    saveRDS(sceptre_object,paste0(sample_name,"_", batch_name, "final_sceptre_object.rds"))

}
