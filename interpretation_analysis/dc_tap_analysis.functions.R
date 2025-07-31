### --- DEFINE CHROMATIN CATEGORIES ---
# get thresholds based on genome-wide elements
get_category_thresholds <- function(enh, quantiles) {
    enh_ctcf <- enh %>% filter(CTCF_peak_overlap == 1) %>% 
        select(cell_type, CTCF.H3K27ac.ratio.CTCF_peak = CTCF.H3K27ac.ratio) %>% 
        pivot_longer(cols = -cell_type, names_to = "feature", values_to = "value") %>% 
        group_by(cell_type, feature) %>% 
        reframe(quantile = quantiles, value = quantile(value, probs = quantiles, na.rm = TRUE))
    
    enh_h3k27ac <- enh %>% filter(H3K27ac_peak_overlap == 1) %>% 
        select(cell_type, H3K27ac.RPM.H3K27ac_peak = H3K27ac.RPM) %>% 
        pivot_longer(cols = -cell_type, names_to = "feature", values_to = "value") %>% 
        group_by(cell_type, feature) %>% 
        reframe(quantile = quantiles, value = quantile(value, probs = quantiles, na.rm = TRUE))

    enh_h3k27me3 <- enh %>% filter(H3K27me3_peak_overlap == 1) %>% 
        select(cell_type, H3K27me3.RPM.expandedRegion.H3K27me3_peak = H3K27me3.RPM.expandedRegion) %>% 
        pivot_longer(cols = -cell_type, names_to = "feature", values_to = "value") %>% 
        group_by(cell_type, feature) %>% 
        reframe(quantile = quantiles, value = quantile(value, probs = quantiles, na.rm = TRUE))

    enh_other <- enh %>% 
        select(cell_type, CTCF.RPM, H3K27ac.RPM, H3K27ac.RPM.expandedRegion, DHS.RPM) %>%
        pivot_longer(cols = -cell_type, names_to = "feature", values_to = "value") %>% 
        group_by(cell_type, feature) %>% 
        reframe(quantile = quantiles, value = quantile(value, probs = quantiles, na.rm = TRUE))

    res <- rbind(enh_ctcf, enh_h3k27ac, enh_h3k27me3, enh_other)
    
    return(res)
}

# get table of quantile values from genome-wide elements
get_threshold_key <- function(thresholds, feature_col, quantile_this) {
    filt <- thresholds %>% filter(feature == feature_col, quantile == quantile_this)
    key <- setNames(filt$value, filt$cell_type)
    return(key)
}

# categorize elements
categorize_elements <- function(enh, thresholds, H3K27ac_q_high = 0.9, H3K27ac_q_low = 0.5) {
    ### CATEGORIZATION LOGIC ###
    # if element overlaps H3K27ac peak:
        # if H3K27ac.RPM.expandedRegion > 90% --> High H3K27ac
        # if H3K27ac.RPM.expandedRegion < 90% --> H3K27ac
    # else:
        # if H3K27ac.RPM.expandedRegion > 90% --> High H3K27ac
        # if H3K27ac.RPM.expandedRegion > 50% --> H3K27ac
        # if element overlaps CTCF peak --> CTCF element
        # if element overlaps H3K27me3 peak --> H3K27me3 element
        # else: No H3K27ac
        

    key_high <- get_threshold_key(thresholds, "H3K27ac.RPM.expandedRegion", H3K27ac_q_high)
    key_low <- get_threshold_key(thresholds, "H3K27ac.RPM.expandedRegion", H3K27ac_q_low)

    enh <- enh %>%
         mutate(element_category = case_when(
                H3K27ac_peak_overlap == 1 & H3K27ac.RPM.expandedRegion >= key_high[cell_type] ~ "High H3K27ac",
                H3K27ac_peak_overlap == 1 ~ "H3K27ac",
                H3K27ac.RPM.expandedRegion >= key_high[cell_type] ~ "High H3K27ac",
                H3K27ac.RPM.expandedRegion >= key_low[cell_type] ~ "H3K27ac",
                CTCF_peak_overlap == 1 ~ "CTCF element",
                H3K27me3_peak_overlap == 1 ~ "H3K27me3 element",
                TRUE ~ "No H3K27ac"))
    return(enh)
}

# summarize properties of genome-wide element-gene pairs
annotate_genomewide_pairs <- function(enh, e2g_files, cell_types, remove_promoters, distance_threshold) {
    res_list <- vector("list", length(cell_types))
    res_list_genes <- vector("list", length(cell_types))

    for (i in seq_along(cell_types)) {
        ct <- cell_types[i]
        print(ct)

        enh_ct <- filter(enh, cell_type == ct)
        pairs_file <- e2g_files[ct]; print(pairs_file)

        pairs <- fread(pairs_file, sep = "\t") %>% 
            select(chr, start, end, class, TargetGene, ubiquitousExpressedGene, distanceToTSS, CellType) %>% 
            filter(distanceToTSS < distance_threshold) %>%
            mutate(distance_category = case_when(distanceToTSS < 10e3 ~  "0-10 kb",
                                             distanceToTSS < 100e3 ~ "10-100 kb",
                                             distanceToTSS < 250e3 ~ "100-250 kb",
                                             distanceToTSS < 1000e3 ~ "250 kb-1 Mb",
                                             distanceToTSS < 2000e3 ~ "1 Mb-2 Mb",
                                             TRUE ~ ">2 Mb"),
                    ubiq_category = ifelse(ubiquitousExpressedGene %in% c("True", TRUE), "Ubiq. expr. gene", "Other gene"))

        if (remove_promoters) {
            pairs <- filter(pairs, class != "promoter")
        }

        res_list[[i]] <- left_join(pairs, enh_ct, by = c("chr", "start", "end")) %>% 
            group_by(cell_type, element_category, distance_category, ubiq_category) %>% 
            summarize(n_pairs = n())

        res_list_genes[[i]] <- pairs %>% select(cell_type = CellType, TargetGene, ubiq_category) %>% distinct() %>% 
            group_by(cell_type, ubiq_category) %>%
            summarize(n_e2g_genes = n())
    }

    res <- rbindlist(res_list) %>% as.data.frame()
    res_genes <- rbindlist(res_list_genes) %>% as.data.frame()

    return(list(res, res_genes))
}

### --- PLOT PROPORTIONS OF ELEMENT-GENE PAIRS ---
# helper function to format numbers nicely
format_number <- function(values) {
  vapply(values, function(value) {
    if (is.na(value)) return(NA_character_)
    
    if (value >= 100000) {
      # Scientific notation
      formatC(value, digits = 1, format = "E") %>% trimws()
    } else if (value %% 1 == 0) {
      # Integer < 100k
      format(value, big.mark = ",") %>% trimws()
    } else {
      # Float < 100k
      formatC(round(value, 1), format = "f", big.mark = ",", digits = 1)
    }
  }, character(1))
}

# plot proportions 
plot_category_proportion_by_dataset_sceptre <- function(enh, crispr, category_col = "element_category",
    direct_effect_weighted = FALSE, direct_effect_threshold = 0.5) {

    # plotting params
    if (category_col == "element_category") {
        category_names <- c("H3K27me3 element", "CTCF element", "High H3K27ac", "H3K27ac", "No H3K27ac")
        cp <- c("#429130", "#49bcbc", "#c5373d", "#d9694a", "#c5cad7")
        names(cp) <- category_names
     } else if (category_col == "distance_category") {
        cp <- c("#002359", "#00488d", "#006eae", "#5496ce", "#9bcae9", "#9bcae9")
        names(cp) <- c("0-10 kb", "10-100 kb", "100-250 kb", "250 kb-1 Mb", "> 1 Mb", "1 Mb-2 Mb")
    } else if (category_col == "ubiq_category") {
        cp <- c("#792374", "#b778b3")
        names(cp) <- c("Ubiq. expr. gene", "Other gene")
    } 

    ct_include <- c("K562_WTC11", unique(crispr$cell_type))
    ds_include <- c("Gasperini2019", "Nasser2021", "Schraivogel2020", "K562_DC_TAP", "WTC11_DC_TAP", "DC_TAP")

    cat_order <- c(gw = "Genome-wide\nE-G pairs", crispr_tested = "All CRISPR\ntested pairs", crispr_powered = "Well-powered\nCRISPR pairs",
        crispr_pos_downreg = "CRISPR+\ndownregulated", crispr_pos_downreg_weighted = "Weighted CRISPR+\ndownregulated",
        crispr_pos_upreg ="CRISPR+\nupregulated", crispr_pos_upreg_weighted ="Weighted CRISPR+\nupregulated")
    ct_order <- c(K562_WTC11 = "K562 + WTC11", K562_WTC11_validation = "K562 + WTC11 (DC-TAP)",
        K562 = "K562", K562_training = "K562 (previous)", K562_validation = "K562 (DC-TAP)",
        WTC11 = "WTC11",  WTC11_validation = "WTC11 (DC-TAP)")

    text_angle <- 45

    # prepare genome-wide prediciton data
    enh <- enh %>%
        filter(cell_type %in% ct_include) %>% 
        mutate(cell_type_label = cell_type, category_use = !!sym(category_col)) %>% 
        group_by(cell_type_label, category_use) %>% 
        summarize(n_pairs = sum(n_pairs), .groups = "drop") %>%
        mutate(plot_category = cat_order["gw"])

    # all crispr
    crispr_all <- crispr %>% 
        mutate(cell_type_label = paste0(cell_type, "_", data_category), category_use = !!sym(category_col)) %>%
        select(cell_type_label, category_use) %>% 
        mutate(plot_category = cat_order["crispr_tested"]) %>% 
        group_by(cell_type_label, category_use, plot_category) %>% 
        summarize(n_pairs = n(), .groups = "drop")
    
    # crispr well-powered
    crispr_wp <- crispr %>% 
        filter(WellPowered == TRUE) %>% 
        mutate(cell_type_label = paste0(cell_type, "_", data_category), category_use = !!sym(category_col)) %>%
        select(cell_type_label, category_use) %>% 
        mutate(plot_category = cat_order["crispr_powered"]) %>% 
        group_by(cell_type_label, category_use, plot_category) %>% 
        summarize(n_pairs = n(), .groups = "drop")

    # crispr pos
    crispr_pos_downreg <- crispr %>% 
        mutate(cell_type_label = paste0(cell_type, "_", data_category), category_use = !!sym(category_col)) %>%
        filter(WellPowered == TRUE, Regulated == 1, direct_vs_indirect > direct_effect_threshold) %>%
        select(cell_type_label, category_use) %>% 
        mutate(plot_category = cat_order["crispr_pos_downreg"]) %>% 
        group_by(cell_type_label, category_use, plot_category) %>% 
        summarize(n_pairs = n(), .groups = "drop")

    crispr_pos_downreg_weighted <- crispr %>% 
        mutate(cell_type_label = paste0(cell_type, "_", data_category), category_use = !!sym(category_col)) %>%
        filter(WellPowered == TRUE, Regulated == 1) %>%
        select(cell_type_label, category_use, direct_vs_indirect, EffectSize) %>% 
        mutate(plot_category = cat_order["crispr_pos_downreg_weighted"]) %>% 
        group_by(cell_type_label, category_use, plot_category) %>% 
        summarize(n_pairs = sum(direct_vs_indirect), .groups = "drop")
    
    crispr_pos_upreg <- crispr %>% 
        mutate(cell_type_label = paste0(cell_type, "_", data_category), category_use = !!sym(category_col)) %>%
        filter(WellPowered == TRUE, Significant == TRUE, EffectSize > 0, direct_vs_indirect > direct_effect_threshold) %>%
        select(cell_type_label, category_use) %>% 
        mutate(plot_category = cat_order["crispr_pos_upreg"]) %>% 
        group_by(cell_type_label, category_use, plot_category) %>% 
        summarize(n_pairs = n(), .groups = "drop")
    
    crispr_pos_upreg_weighted <- crispr %>% 
        mutate(cell_type_label = paste0(cell_type, "_", data_category), category_use = !!sym(category_col)) %>%
        filter(WellPowered == TRUE, Significant == TRUE, EffectSize > 0) %>%
        select(cell_type_label, category_use, direct_vs_indirect, EffectSize) %>% 
        mutate(plot_category = cat_order["crispr_pos_upreg_weighted"]) %>% 
        group_by(cell_type_label, category_use, plot_category) %>% 
        summarize(n_pairs = sum(direct_vs_indirect), .groups = "drop")
    
    if (direct_effect_weighted) {
        crispr_pos_downreg <- crispr_pos_downreg_weighted
        crispr_pos_upreg <- crispr_pos_upreg_weighted
    }

    res <- rbind(enh, crispr_all, crispr_wp, crispr_pos_downreg) %>% #, crispr_pos_upreg) %>% 
        mutate(category_use = factor(category_use, levels = names(cp), ordered = TRUE),
            plot_category = factor(plot_category, levels = cat_order, ordered = TRUE),
            cell_type_label = ct_order[cell_type_label],
            cell_type_label = factor(cell_type_label, levels = ct_order, ordered = TRUE))

    # calculate totals per category
    dataset_totals <- res %>%
        group_by(plot_category, cell_type_label) %>%
        summarise(total_count = sum(n_pairs), .groups = "drop")

    # calculate proportions and prepare data for plotting

    prop_data <- res %>%
        #group_by(plot_category, cell_type_label, element_category) %>%
        mutate(category_count = n_pairs) %>%
        left_join(dataset_totals, by = c("plot_category", "cell_type_label")) %>%
        mutate(proportion = category_count / total_count,
            y_position = cumsum(proportion) - proportion / 2,
            category_use = factor(category_use, levels = names(cp), ordered = TRUE),
            plot_category = factor(plot_category, levels = cat_order, ordered = TRUE),
            cell_type_label = factor(cell_type_label, levels = ct_order, ordered = TRUE),
            category_count = format_number(category_count),
            total_count = format_number(total_count))
            #category_count = ifelse(category_count < 100e3, format_int(category_count), format_sci(category_count)),
            #total_count = ifelse(total_count < 100e3, format_int(total_count), format_sci(total_count)))

    # create stacked bar plot
    x_lab <- ""

    p <- ggplot(prop_data, aes(x = plot_category, y = proportion, fill = category_use)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = category_count), position = position_stack(vjust = 0.5), color = "#000000", size = 2.5) +
        geom_text(aes(y = -0.05, label = paste0("(", total_count, ")")), size = 3, color = "#000000", vjust = 0) +
        labs(x = x_lab, fill = "Element category", y = "Proportion of pairs with element in category") +
        facet_grid(. ~ cell_type_label, scales = "free", space = "free") +
        scale_fill_manual(values = cp) +
        theme_classic() + theme(strip.background = element_blank(), panel.grid = element_blank(),
            axis.text = element_text(size = 10, color = "#000000"), axis.text.x = element_text(angle = text_angle, hjust = 1, vjust = 1),
            axis.title = element_text(size = 12), axis.ticks = element_line(color = "#000000"), legend.position = "right",
            plot.title = element_text(size = 14, color = "#000000"), plot.subtitle = element_text(size = 12, color = "#000000"))

    return(p)
}

### --- COMPARE METRICS BETWEEN GROUPS OF DE-G PAIRS ---
# compute summary values and plot comparisons 
plot_percent_positive_enrichment_effect_size_combined <- function(crispr, group_var, direct_effect_weighted, direct_effect_threshold, all_power, out_prefix) {
  sig_vars <- c("Significant", "Regulated")

  if (direct_effect_weighted & !is.null(direct_effect_threshold)) {stop("Choose at most one of direct effect weighting or filtering!")}

  plots <- list()
  combined_pct_pos <- list()
  combined_enr <- list()
  power_enr <- list()
  combined_es <- list()
  max_prop <- list()
  max_prop_pairs <- list()
  max_es <- list()

  ds_cp <- c(K562_DC_TAP = "#006eae", WTC11_DC_TAP = "#00488d", DC_TAP = "#005a9d", Combined = "#435369",
    Gasperini2019 = "#d3a9ce", Nasser2021 = "#b778b3", Schraivogel2020 = "#a64791",
    Klann = "#d3a9ce", Morris = "#d3a9ce", Xie = "#d3a9ce") 

  for (significance_var in sig_vars) {
    z <- qnorm(0.05/2, lower.tail=FALSE)
    pos_dodge <- 0.9

    summarize_pairs_vars <- c("elementName", "measuredGeneSymbol")
    summarize_pairs_label <- "E-G pairs"

    if (group_var == "element_category") {
      category_names <- c("H3K27me3 element", "CTCF element", "High H3K27ac", "H3K27ac", "No H3K27ac")
      cp <- c("#429130", "#49bcbc", "#c5373d", "#d9694a", "#c5cad7")
      names(cp) <- category_names
      label <- "Element category"
      order_levels <- rev(names(cp))
      summarize_vars <- c("elementName")
      summarize_label <- "elements"

    } else if (group_var == "enhancerness") {
      cp <- c(`H3K27ac+ element` = "#d9694a", `Other element` = "#435369")
      label <- "Element type"
      order_levels <- rev(names(cp))
      summarize_vars <- c("elementName")
      summarize_label <- "elements" 

    } else if (group_var == "ubiq_category") {
      cp <- c("#792374", "#b778b3")
      names(cp) <- c("Ubiq. expr. gene", "Other gene")
      label <- "Promoter class"
      order_levels <- names(cp)
      summarize_vars <- c("measuredGeneSymbol")
      summarize_label <- "genes"

    } else if (group_var == "distance_category") {
      cp <- c("#002359", "#00488d", "#006eae", "#5496ce", "#9bcae9")
      names(cp) <- c("0-10 kb", "10-100 kb", "100-250 kb", "250 kb-1 Mb", "1 Mb-2 Mb")
      label <- "Distance to TSS"
      order_levels <- names(cp)
      summarize_vars <- summarize_pairs_vars
      summarize_label <- summarize_pairs_label

    } else {
      stop("Unsupported group_var")
    }

    ## prepare data
    if (direct_effect_weighted) {
        crispr_this <- crispr %>%
            mutate(positive_indicator = ifelse(!!sym(significance_var) == TRUE, direct_vs_indirect, 0))
        if (all_power) {
            crispr_powered <- crispr_this
        } else {
            crispr_powered <- crispr_this %>% filter(WellPowered == TRUE)
        }
        

    } else { # filter by threshold (not significant OR significant with high direct rate
        crispr_this <- crispr %>%
            filter(!(!!sym(significance_var) == TRUE & direct_vs_indirect < direct_effect_threshold)) %>%
            mutate(positive_indicator = ifelse(!!sym(significance_var) == TRUE, 1, 0))
        crispr_powered <- crispr_this #%>% filter(WellPowered == TRUE)
    }

    # hit rate
    pct_pos <- crispr_powered %>% 
      group_by(Dataset, pick(all_of(summarize_vars))) %>%
      #mutate(anyRegulated = any(!!sym(significance_var))) %>%
      mutate(anyRegulated = max(positive_indicator)) %>% 
      ungroup() %>%
      select(all_of(summarize_vars), Dataset, anyRegulated, category = !!sym(group_var)) %>% 
      distinct() %>% 
      group_by(Dataset, category) %>% 
      summarize(n_tested_category = n(), n_positive_category = sum(anyRegulated), .groups = "drop") %>%
      mutate(prop_positive_elements = n_positive_category / n_tested_category,
             prop_adjust = (n_positive_category + 2) / (n_tested_category + 4),
             SE_prop = sqrt(prop_adjust * (1 - prop_adjust) / (n_tested_category + 4)),
             CI_prop_low = pmax(prop_positive_elements - z * SE_prop, 0),
             CI_prop_high = prop_positive_elements + z * SE_prop,
             n_label = paste0("(", format_number(n_positive_category), ")"),
             category = factor(category, levels = order_levels, ordered = TRUE),
             Dataset = factor(Dataset, levels = names(ds_cp), ordered = TRUE),
             significance = significance_var,
             grouped_by = summarize_label) %>%
      filter(!is.na(category))
    pct_pos$y_label_pos <- max(pct_pos$CI_prop_high, na.rm = TRUE) * 1.1

    # hit rate in terms of pairs
    pct_pos_pairs <- crispr_powered %>% 
      group_by(Dataset, pick(all_of(summarize_pairs_vars))) %>%
      #mutate(anyRegulated = any(!!sym(significance_var))) %>%
      mutate(anyRegulated = max(positive_indicator)) %>% 
      ungroup() %>%
      select(all_of(summarize_pairs_vars), Dataset, anyRegulated, category = !!sym(group_var)) %>% 
      distinct() %>% 
      group_by(Dataset, category) %>% 
      summarize(n_tested_category = n(), n_positive_category = sum(anyRegulated), .groups = "drop") %>%
      mutate(prop_positive_elements = n_positive_category / n_tested_category,
             prop_adjust = (n_positive_category + 2) / (n_tested_category + 4),
             SE_prop = sqrt(prop_adjust * (1 - prop_adjust) / (n_tested_category + 4)),
             CI_prop_low = pmax(prop_positive_elements - z * SE_prop, 0),
             CI_prop_high = prop_positive_elements + z * SE_prop,
             n_label = paste0("(", format_number(n_positive_category), ")"),
             category = factor(category, levels = order_levels, ordered = TRUE),
             Dataset = factor(Dataset, levels = names(ds_cp), ordered = TRUE),
             significance = significance_var,
             grouped_by = summarize_pairs_label) %>%
      filter(!is.na(category))
    pct_pos_pairs$y_label_pos <- max(pct_pos_pairs$CI_prop_high, na.rm = TRUE) * 1.1

    # enrichment
    enr <- crispr_powered %>% 
      group_by(Dataset) %>%
      #mutate(n_tested = n(), n_positive = sum(!!sym(significance_var))) %>%
      mutate(n_tested = n(), n_positive = sum(positive_indicator)) %>% 
      group_by(Dataset, category = !!sym(group_var), n_tested, n_positive) %>%
      summarize(n_tested_category = n(),
        n_positive_category = sum(positive_indicator), # sum(!!sym(significance_var)),
        .groups = "drop") %>% 
      mutate(prop_tested = n_tested_category/n_tested,
             prop_pos = ifelse(n_positive > 0, n_positive_category/n_positive, 0),
             enrichment = ifelse(prop_tested > 0, prop_pos/prop_tested, NA),
             enrichment_label = round(enrichment, 2),
             SE_log_enr = ifelse(n_positive_category > 0 & n_tested_category > 0,
                                 sqrt(((n_positive - n_positive_category) / n_positive_category) / n_positive) +
                                 ((n_tested - n_tested_category) / n_tested_category) / n_tested,
                                 NA),
             CI_enr_low = exp(log(enrichment) - z * SE_log_enr),
             CI_enr_high = exp(log(enrichment) + z * SE_log_enr),
             p_enr = phyper(n_positive_category, n_positive, n_tested, n_positive_category + n_tested_category, log.p = FALSE, lower.tail = FALSE),
             p_adjust_enr = p.adjust(p_enr, method = "bonferroni"),
             sign_label = case_when(p_adjust_enr < 0.001 ~ "***", p_adjust_enr < 0.01 ~ "**", p_adjust_enr < 0.05 ~ "*", TRUE ~ ""),
             category = factor(category, levels = order_levels, ordered = TRUE),
             Dataset = factor(Dataset, levels = names(ds_cp), ordered = TRUE),
             significance = significance_var) %>%
      filter(!is.na(category))

    # effect size
    es_col <- ifelse(significance_var == "Regulated", "EffectSize", "absEffectSize")
    es <- crispr_powered %>% 
        filter(!!sym(significance_var) == 1) %>% 
        mutate(category = !!sym(group_var),
            absEffectSize = abs(EffectSize),
            category = factor(category, levels = order_levels, ordered = TRUE),
            Dataset = factor(Dataset, levels = names(ds_cp), ordered = TRUE),
            plotEffectSize = !!sym(es_col))

    es_smry <- es %>% 
        group_by(Dataset, category) %>% 
        summarize(n_pairs = n(),
            n_label = paste0("(", n_pairs, ")"),
            mean_EffectSize = mean(EffectSize, na.rm = TRUE), 
            mean_absEffectSize = mean(absEffectSize, na.rm = TRUE),
            sd_absEffectSize = sd(absEffectSize, na.rm = TRUE)) %>% ungroup() %>%
        mutate(se_absEffectSize = sd_absEffectSize / sqrt(n_pairs),
            CI_ES_low = mean_absEffectSize - z * se_absEffectSize,
            CI_ES_high = mean_absEffectSize + z * se_absEffectSize,
            significance = significance_var)
    es_smry$y_label_pos <- min(100, max(es$plotEffectSize, na.rm = TRUE) + 3)

    ## plot
    # params
    max_prop[[significance_var]] <- max(pct_pos$CI_prop_high, na.rm = TRUE)
    max_prop_pairs[[significance_var]] <- max(pct_pos_pairs$CI_prop_high, na.rm = TRUE)

    y_label1 <- ifelse(significance_var == "Regulated",
                       paste0("Fraction tested ", summarize_label, "\nwith 1+ downregulated hit (# hits)"),
                       paste0("Fraction tested ", summarize_label, "\nwith 1+ up/downregulated hit (# hits)"))
    
    y_label4 <- ifelse(significance_var == "Regulated",
                       paste0("Fraction tested ", summarize_pairs_label, "\nwith 1+ downregulated hit (# hits)"),
                       paste0("Fraction tested ", summarize_pairs_label, "\nwith 1+ up/downregulated hit (# hits)"))

    y_label2 <- ifelse(significance_var == "Regulated",
                       "Enrichment of downregulated hits\n(% significant / % tested)",
                       "Enrichment of up/downregulated hits\n(% significant / % tested)")
                    
    max_es[[significance_var]] <- min(100, max(es$plotEffectSize, na.rm = TRUE))
    y_label3 <- ifelse(significance_var == "Regulated",
                       "% effect size of downregulated hits (# hits)",
                       "Abs (% effect size) of up/downregulated hits (# hits)")
    # hit rate
    pct_pos_plot <- filter(pct_pos, prop_positive_elements > 0)
    p1 <- ggplot(pct_pos_plot, aes(x = Dataset, y = prop_positive_elements, fill = category)) + 
      geom_col(position = position_dodge(width = pos_dodge)) +
      geom_linerange(aes(ymin = CI_prop_low, ymax = CI_prop_high), position = position_dodge(width = pos_dodge), linewidth = 0.5, color = "black") +
      geom_text(aes(y = y_label_pos, label = n_label), size = 2, position = position_dodge(width = pos_dodge)) +
      ylim(c(0, max_prop[[significance_var]] * 1.15)) +
      scale_fill_manual(values = cp) +
      labs(x = "Dataset", y = y_label1, fill = label) +
      theme_classic() +
      theme(axis.text = element_text(size = 8, color = "#000000"), axis.title = element_text(size = 9),
        axis.ticks = element_line(color = "#000000"), axis.ticks.x = element_blank(),
        legend.position = "none")
    
    # hit rate, pairs
    pct_pos_pairs_plot <- filter(pct_pos_pairs, prop_positive_elements > 0)
    p4 <- ggplot(pct_pos_pairs_plot, aes(x = Dataset, y = prop_positive_elements, fill = category)) + 
      geom_col(position = position_dodge(width = pos_dodge)) +
      geom_linerange(aes(ymin = CI_prop_low, ymax = CI_prop_high), position = position_dodge(width = pos_dodge), linewidth = 0.5, color = "black") +
      geom_text(aes(y = y_label_pos, label = n_label), size = 2, position = position_dodge(width = pos_dodge)) +
      ylim(c(0, max_prop_pairs[[significance_var]] * 1.15)) +
      scale_fill_manual(values = cp) +
      labs(x = "Dataset", y = y_label4, fill = label) +
      theme_classic() +
      theme(axis.text = element_text(size = 8, color = "#000000"), axis.title = element_text(size = 9),
        axis.ticks = element_line(color = "#000000"), axis.ticks.x = element_blank(),
        legend.position = "none")

    # enrichment
    enr_plot <- enr %>% filter(n_positive_category > 0, CI_enr_high < 100)
    p2 <- ggplot(enr_plot, aes(x = Dataset, y = enrichment, fill = category)) + 
      geom_hline(yintercept = 1, linetype = "dashed", color = "#c5cad7") +
      geom_col(position = position_dodge(width = pos_dodge)) +
      geom_linerange(aes(ymin = CI_enr_low, ymax = CI_enr_high), position = position_dodge(width = pos_dodge), linewidth = 0.5, color = "black") +
      geom_text(aes(y = CI_enr_high + 0.2, label = sign_label), size = 2, position = position_dodge(width = pos_dodge)) +
      scale_fill_manual(values = cp) +
      labs(x = "Dataset", y = y_label2, fill = label) +
      theme_classic() +
      theme(axis.text = element_text(size = 8, color = "#000000"), axis.title = element_text(size = 9),
        axis.ticks = element_line(color = "#000000"), axis.ticks.x = element_blank(),
        legend.position = "none")
    
    # effect size
    p3 <- ggplot(es, aes(x = Dataset, y = plotEffectSize)) + 
        geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed", color = "#c5cad7") +
        geom_boxplot(aes(color = category), fill = NA, width = 0.7, outlier.shape = 16, outlier.size = 0.75, position = position_dodge(pos_dodge)) +
        geom_text(data = es_smry, aes(y = y_label_pos, group = category, label = n_label), size = 2, color = "black", position = position_dodge(pos_dodge)) +
        ylim(c(NA, max_es[[significance_var]]) + 6) +
		scale_color_manual(values = cp) + 
		labs(x = "Dataset", y = y_label3) + 
      theme_classic() +
      theme(axis.text = element_text(size = 8, color = "#000000"), axis.title = element_text(size = 9),
        axis.ticks = element_line(color = "#000000"), axis.ticks.x = element_blank(),
        legend.position = "none")

    plots[[significance_var]] <- list(p1, p2, p3, p4)
    combined_pct_pos[[significance_var]] <- rbind(pct_pos, pct_pos_pairs) %>% select(-n_label, -prop_adjust)
    combined_enr[[significance_var]] <- enr %>% select(-sign_label)
    combined_es[[significance_var]] <- es_smry %>% select(-n_label)
  }

    ## power enrichment
    enr_power <- crispr %>% 
        group_by(Dataset) %>%
        mutate(n_tested = n(), n_powered = sum(WellPowered)) %>% 
        group_by(Dataset, category = !!sym(group_var), n_tested, n_powered) %>%
        summarize(n_tested_category = n(), n_powered_category = sum(WellPowered), .groups = "drop") %>% 
        mutate(prop_tested = n_tested_category/n_tested,
            prop_pos = ifelse(n_powered > 0, n_powered_category / n_powered, 0),
            enrichment = ifelse(prop_tested > 0, prop_pos/prop_tested, NA),
            enrichment_label = round(enrichment, 2),
            SE_log_enr = ifelse(n_powered_category > 0 & n_tested_category > 0,
                                sqrt(((n_powered - n_powered_category) / n_powered_category) / n_powered) +
                                ((n_tested - n_tested_category) / n_tested_category) / n_tested,
                                NA),
            CI_enr_low = exp(log(enrichment) - z * SE_log_enr),
            CI_enr_high = exp(log(enrichment) + z * SE_log_enr),
            p_enr = phyper(n_powered_category, n_powered, n_tested, n_powered_category + n_tested_category, log.p = FALSE, lower.tail = FALSE),
            p_adjust_enr = p.adjust(p_enr, method = "bonferroni"),
            sign_label = case_when(p_adjust_enr < 0.001 ~ "***", p_adjust_enr < 0.01 ~ "**", p_adjust_enr < 0.05 ~ "*", TRUE ~ ""),
            category = factor(category, levels = order_levels, ordered = TRUE),
            Dataset = factor(Dataset, levels = names(ds_cp), ordered = TRUE),
            significance = significance_var) %>%
    filter(!is.na(category))

    # make plt for power enrichment
    p5 <- ggplot(enr_power, aes(x = Dataset, y = enrichment, fill = category)) + 
      geom_hline(yintercept = 1, linetype = "dashed", color = "#c5cad7") +
      geom_col(position = position_dodge(width = pos_dodge)) +
      geom_linerange(aes(ymin = CI_enr_low, ymax = CI_enr_high), position = position_dodge(width = pos_dodge), linewidth = 0.5, color = "black") +
      geom_text(aes(y = CI_enr_high + 0.2, label = sign_label), size = 2, position = position_dodge(width = pos_dodge)) +
      scale_fill_manual(values = cp) +
      labs(x = "Dataset", y = "Enrichment\n(% well-powered / % tested)", fill = label) +
      theme_classic() +
      theme(axis.text = element_text(size = 8, color = "#000000"), axis.title = element_text(size = 9),
        axis.ticks = element_line(color = "#000000"), axis.ticks.x = element_blank(),
        legend.position = "top")

    binary_groups <- c("ubiq_category", "enhancerness")
    n_dataset <- length(unique(crispr$Dataset))
    w <- ifelse(group_var %in% binary_groups, 12 + (n_dataset - 2) * 3, 16 + (n_dataset - 2) * 4)

    ## save plots
    ggsave(paste0(out_prefix, group_var, "_power_enrichment.pdf"), p5, height = 4, w = w/4)

    # plot combined others
  legend <- cowplot::get_plot_component(plots[["Significant"]][[1]] + theme(legend.position = "top"), 'guide-box-top')

  combined <- plot_grid(
    legend,
    plot_grid(plots[["Significant"]][[3]], plots[["Significant"]][[2]], plots[["Significant"]][[1]], plots[["Significant"]][[4]], nrow = 1, rel_widths = c(1, 1, 1, 1)),
    plot_grid(plots[["Regulated"]][[3]], plots[["Regulated"]][[2]],  plots[["Regulated"]][[1]], plots[["Regulated"]][[4]], nrow = 1, rel_widths = c(1, 1, 1, 1, 1)),
    ncol = 1, rel_heights = c(0.15, 1, 1))

  ggsave(paste0(out_prefix, group_var, "_combined_summary.pdf"), combined, height = 8, width = w)

  fwrite(enr_power, paste0(out_prefix, group_var, "_power_enrichment.tsv"), sep = "\t")
  fwrite(bind_rows(combined_pct_pos), paste0(out_prefix, group_var, "_crispr_results_combined_percent_positives.tsv"), sep = "\t")
  fwrite(bind_rows(combined_enr), paste0(out_prefix, group_var, "_crispr_results_combined_enrichment.tsv"), sep = "\t")
  fwrite(bind_rows(combined_es), paste0(out_prefix, group_var, "_crispr_results_combined_effect_sizes.tsv"), sep = "\t")
}

# compute significance values between groups
compare_metrics_across_groups <- function(crispr, group_var, out_prefix, pairs = c("Category", "Dataset"),
    direct_effect_weighted, direct_effect_threshold, all_power) {
  if (direct_effect_weighted & !is.null(direct_effect_threshold)) {stop("Choose either direct effect weighting or filtering!")}

  pairs <- match.arg(pairs)
  sig_vars <- c("Significant", "Regulated")
  summarize_pairs_vars <- c("elementName", "measuredGeneSymbol")
  summarize_pairs_label <- "pairs"

  if (group_var == "element_category") {
    category_names <- c("H3K27me3 element", "CTCF element", "High H3K27ac", "H3K27ac", "No H3K27ac")
    cp <- c("#429130", "#49bcbc", "#c5373d", "#d9694a", "#c5cad7")
    names(cp) <- category_names
    label <- "Element category"
    order_levels <- rev(names(cp))
    summarize_vars <- c("elementName")
    summarize_label <- "elements"
  } else if (group_var == "enhancerness") {
    cp <- c(`H3K27ac+ element` = "#D9694A", `Other element` = "#435369")
    label <- "Element type"
    order_levels <- rev(names(cp))
    summarize_vars <- c("elementName")
    summarize_label <- "elements"
  } else if (group_var == "ubiq_category") {
    cp <- c("#792374", "#b778b3")
    names(cp) <- c("Ubiq. expr. gene", "Other gene")
    label <- "Promoter class"
    order_levels <- names(cp)
    summarize_vars <- c("measuredGeneSymbol")
    summarize_label <- "genes"
  } else if (group_var == "distance_category") {
    cp <- c("#002359", "#00488d", "#006eae", "#5496ce", "#9bcae9")
    names(cp) <- c("0-10 kb", "10-100 kb", "100-250 kb", "250 kb-1 Mb", "1 Mb-2 Mb")
    label <- "Distance to TSS"
    order_levels <- names(cp)
    summarize_vars <- summarize_pairs_vars
    summarize_label <- summarize_pairs_label
  } else {
    stop("Unsupported group_var")
  }

  all_results <- list()

  for (significance_var in sig_vars) {

    ## prepare data
    if (direct_effect_weighted) {
        crispr_this <- crispr %>%
            mutate(positive_indicator = ifelse(!!sym(significance_var) == TRUE, direct_vs_indirect, 0))

    } else { # filter by threshold (not significant OR significant with high direct rate
        crispr_this <- crispr %>%
            filter(!(!!sym(significance_var) == TRUE & direct_vs_indirect < direct_effect_threshold)) %>%
            mutate(positive_indicator = ifelse(!!sym(significance_var) == TRUE, 1, 0))
        
    }

    if (all_power) {
        crispr_powered <- crispr_this
    } else {
        crispr_powered <- crispr_this %>% filter(WellPowered == TRUE)
    }

    crispr_sub <- crispr_powered %>%
      mutate(category = !!sym(group_var),
             absEffectSize = abs(EffectSize)) %>%
      filter(!is.na(category))

    hit_data <- crispr_sub %>%
      group_by(Dataset, pick(all_of(summarize_vars))) %>%
      #mutate(anyRegulated = any(!!sym(significance_var))) %>%
      mutate(anyRegulated = max(positive_indicator)) %>%
      ungroup() %>%
      select(all_of(summarize_vars), Dataset, anyRegulated, category) %>%
      distinct() %>%
      group_by(Dataset, category) %>%
      summarize(n_total = n(), n_hit = sum(anyRegulated), hit_rate = n_hit / n_total, .groups = "drop")
    
    hit_data_pairs <- crispr_sub %>%
      group_by(Dataset, pick(all_of(summarize_pairs_vars))) %>%
      mutate(anyRegulated = max(positive_indicator)) %>%
      ungroup() %>%
      select(all_of(summarize_vars), Dataset, anyRegulated, category) %>%
      distinct() %>%
      group_by(Dataset, category) %>%
      summarize(n_total = n(), n_hit = sum(anyRegulated), hit_rate = n_hit / n_total, .groups = "drop")

    hit_data_list <- list(smry = hit_data, pairs = hit_data_pairs)
    hit_pairs_list <- list()
    for (hd in names(hit_data_list)) {
        if (pairs == "Category") {
            hit_pairs_list[[hd]] <- hit_data_list[[hd]] %>%
                group_by(Dataset) %>%
                filter(n_distinct(category) >= 2) %>%
                group_modify(~{
                combs <- combn(unique(.x$category), 2, simplify = FALSE)
                do.call(rbind, lapply(combs, function(pair) {
                    d1 <- .x[.x$category == pair[1], ]
                    d2 <- .x[.x$category == pair[2], ]
                    if (nrow(d1) == 0 | nrow(d2) == 0) return(NULL)
                    x <- c(d1$n_hit, d2$n_hit)
                    n <- c(d1$n_total, d2$n_total)
                    p <- tryCatch(prop.test(x, n)$p.value, error = function(e) NA_real_)
                    data.frame(group1 = pair[1], group2 = pair[2],
                            group1_value = d1$hit_rate, group2_value = d2$hit_rate,
                            p_value = p, test = "prop.test")
                }))
                }, .groups = "drop") %>%
                mutate(metric = ifelse(hd == "smry", paste0("HitRate_", summarize_label),
                    paste0("HitRate_", summarize_pairs_label)))
        } else {
            hit_pairs_list[[hd]] <- hit_data_list[[hd]] %>%
                group_by(category) %>%
                filter(n_distinct(Dataset) >= 2) %>%
                group_modify(~{
                combs <- combn(unique(.x$Dataset), 2, simplify = FALSE)
                do.call(rbind, lapply(combs, function(pair) {
                    d1 <- .x[.x$Dataset == pair[1], ]
                    d2 <- .x[.x$Dataset == pair[2], ]
                    if (nrow(d1) == 0 | nrow(d2) == 0) return(NULL)
                    x <- c(d1$n_hit, d2$n_hit)
                    n <- c(d1$n_total, d2$n_total)
                    p <- tryCatch(prop.test(x, n)$p.value, error = function(e) NA_real_)
                    data.frame(group1 = pair[1], group2 = pair[2],
                            group1_value = d1$hit_rate, group2_value = d2$hit_rate,
                            p_value = p, test = "prop.test")
                }))
                }, .groups = "drop") %>%
                mutate(metric = ifelse(hd == "smry", paste0("HitRate_", summarize_label),
                    paste0("HitRate_", summarize_pairs_label)))
        }
    }

    enr_data <- crispr_sub %>%
        group_by(Dataset) %>%
        summarize(n_total = n(), n_sig = sum(positive_indicator), .groups = "drop") %>%
        right_join(crispr_sub %>%
            group_by(Dataset, category) %>%
            summarize(n_cat = n(), n_sig_cat = sum(positive_indicator), .groups = "drop"),
            by = "Dataset") %>%
        mutate(enrichment = (n_sig_cat / n_cat) / (n_sig / n_total),
            SE_log_enr = sqrt(1 / n_sig_cat - 1 / n_cat + 1 / n_sig - 1 / n_total))

    if (pairs == "Category") {
        message("enrichment pairs category")
      enr_pairs <- enr_data %>%
        group_by(Dataset) %>%
        filter(n_distinct(category) >= 2) %>%
        group_modify(~{
          combs <- combn(unique(.x$category), 2, simplify = FALSE)
          do.call(rbind, lapply(combs, function(pair) {
            d1 <- .x[.x$category == pair[1], ]
            d2 <- .x[.x$category == pair[2], ]
            d <- log(d1$enrichment / d2$enrichment)
            SE_d <- sqrt(d1$SE_log_enr^2 + d2$SE_log_enr^2)
            z <- d / SE_d
            p <- pnorm(-abs(z)) * 2
            data.frame(group1 = pair[1], group2 = pair[2],
                       group1_value = d1$enrichment, group2_value = d2$enrichment,
                       p_value = p, test = "log_enrichment_z")
          }))
        }, .groups = "drop") %>%
        mutate(metric = "Enrichment")
    } else {
        message("enrichment pairs dataset")
      enr_pairs <- enr_data %>%
        group_by(category) %>%
        filter(n_distinct(Dataset) >= 2) %>%
        group_modify(~{
          combs <- combn(unique(.x$Dataset), 2, simplify = FALSE)
          do.call(rbind, lapply(combs, function(pair) {
            d1 <- .x[.x$Dataset == pair[1], ]
            d2 <- .x[.x$Dataset == pair[2], ]
            d <- log(d1$enrichment / d2$enrichment)
            SE_d <- sqrt(d1$SE_log_enr^2 + d2$SE_log_enr^2)
            z <- d / SE_d
            p <- pnorm(-abs(z)) * 2
            data.frame(group1 = pair[1], group2 = pair[2],
                       group1_value = d1$enrichment, group2_value = d2$enrichment,
                       p_value = p, test = "log_enrichment_z")
          }))
        }, .groups = "drop") %>%
        mutate(metric = "Enrichment")
    }

    es_data <- crispr_sub %>%
      filter(!!sym(significance_var) == 1) %>%
      mutate(plotEffectSize = if (significance_var == "Regulated") EffectSize else abs(EffectSize)) %>%
      filter(!is.na(plotEffectSize))

    if (pairs == "Category") {
      es_pairs <- es_data %>%
        group_by(Dataset, category) %>%
        summarize(effect_values = list(plotEffectSize), median_effect = median(plotEffectSize), .groups = "drop") %>%
        group_by(Dataset) %>%
        filter(n_distinct(category) >= 2) %>%
        group_modify(~{
          combs <- combn(unique(.x$category), 2, simplify = FALSE)
          do.call(rbind, lapply(combs, function(pair) {
            d1 <- .x$effect_values[.x$category == pair[1]][[1]]
            d2 <- .x$effect_values[.x$category == pair[2]][[1]]
            p <- wilcox.test(d1, d2)$p.value
            data.frame(group1 = pair[1], group2 = pair[2],
                       group1_value = median(d1), group2_value = median(d2),
                       p_value = p, test = "wilcox")
          }))
        }, .groups = "drop") %>%
        mutate(metric = "EffectSize")
    } else {
      es_pairs <- es_data %>%
        group_by(Dataset, category) %>%
        summarize(effect_values = list(plotEffectSize), median_effect = median(plotEffectSize), .groups = "drop") %>%
        group_by(category) %>%
        filter(n_distinct(Dataset) >= 2) %>%
        group_modify(~{
          combs <- combn(unique(.x$Dataset), 2, simplify = FALSE)
          do.call(rbind, lapply(combs, function(pair) {
            d1 <- .x$effect_values[.x$Dataset == pair[1]][[1]]
            d2 <- .x$effect_values[.x$Dataset == pair[2]][[1]]
            p <- wilcox.test(d1, d2)$p.value
            data.frame(group1 = pair[1], group2 = pair[2],
                       group1_value = median(d1), group2_value = median(d2),
                       p_value = p, test = "wilcox")
          }))
        }, .groups = "drop") %>%
        mutate(metric = "EffectSize")
    }


    hit_pairs <- rbindlist(hit_pairs_list) %>% as.data.frame()
    combined <- bind_rows(hit_pairs, enr_pairs, es_pairs) %>%
      mutate(significance = significance_var) %>%
      group_by(metric, significance) %>%
      mutate(p_value_adj = p.adjust(p_value, method = "BH")) %>%
      ungroup()

    all_results[[significance_var]] <- combined
  }

  final <- bind_rows(all_results)
  out_file <- paste0(out_prefix, group_var, ".signficance_table.by_", pairs, ".tsv")
  fwrite(final, out_file, sep = "\t")

}

### --- PLOT ENCODE-rE2G SCORE BY CATEGORY ---
plot_score_by_category <- function(crispr, group_var, significance_var, e2g_threshold, out_prefix, scale_log10 = FALSE) {
	ds_cp <- c(Gasperini2019 = "#d3a9ce", Nasser2021 = "#b778b3", Schraivogel2020 = "#a64791",
		K562_DC_TAP = "#006eae", WTC11_DC_TAP = "#00488d", DC_TAP = "#005a9d") 

	if (group_var == "element_category") {
		category_names <- c("H3K27me3 element", "CTCF element", "High H3K27ac", "H3K27ac", "No H3K27ac")
		cp <- c("#429130", "#49bcbc", "#c5373d", "#d9694a", "#c5cad7")
		names(cp) <- category_names
		label <- "Element category"
		order_levels <- c("High H3K27ac", "H3K27ac", "H3K27me3 element", "CTCF element", "No H3K27ac")
	} else if (group_var == "enhancerness") {
		cp <- c(`H3K27ac+ element` = "#D9694A", `Other element` = "#435369")
		label <- "Element type"
		order_levels <- rev(names(cp))
	} else if (group_var == "CTCF_category") {
		cp <- c(`CTCF element` = "#49bcbc", `Other element` = "#435369")
		label <- "Element type"
		order_levels <- rev(names(cp))
    } else if (group_var == "ubiq_category") {
		cp <- c("#792374", "#b778b3")
		names(cp) <- c("Ubiq. expr. gene", "Other gene")
		label <- "Promoter class"
		order_levels <- names(cp)
	} else if (group_var == "EffectSize_5pct") {
		cp <- c("#9b241c", "#006eae")
		names(cp) <- c("Over 5%", "Under 5%")
		label <- "Effect size"
		order_levels <- names(cp)
	} else if (group_var == "pDirect_90pct") {
		cp <- c("#0e3716", "#429130", "#96a0b3")
		names(cp) <- c("Over 90%", "50-90%", "Under 50%")
		label <- "P (Direct effect)"
		order_levels <- names(cp)
	} else {
		stop("Unsupported group_var")
	}

    ## all pairs
    res_all <- crispr %>%
        mutate(category = !!sym(group_var),
            positive = !!sym(significance_var),
			category = factor(category, levels = order_levels, ordered = TRUE),
			Dataset = factor(Dataset, levels = names(ds_cp), ordered = TRUE))

    # just positives
	res <- res_all %>% 
		filter(positive == 1)

	# summary labels
	smry <- res_all %>% 
		group_by(Dataset, category) %>% 
		summarize(n_total = n(),
            n_positive = sum(positive == 1),
            n_true_positive = sum(pred_value >= e2g_threshold & positive == 1),
            n_false_positive = sum(pred_value >= e2g_threshold & positive == 0),
            n_true_negative = sum(pred_value < e2g_threshold & positive == 0),
            n_false_negative = sum(pred_value < e2g_threshold & positive == 1),
            mean_score_positives = mean(pred_value[positive == 1], na.rm = TRUE),
            .groups = "drop") %>% 
        mutate(
            precision = n_true_positive / (n_true_positive + n_false_positive),
            recall = n_true_positive / n_positive,
            recall_pct = round(100 * recall, 1),
            recall_label = paste0(recall_pct, "%"),
			n_label = paste0(n_positive))

	# Compute p-values if only 2 groups
	pval_df <- res %>%
		group_by(Dataset) %>%
		filter(n_distinct(category) == 2) %>%
		summarize(
			pval = tryCatch({
				wilcox.test(pred_value ~ category)$p.value
			}, error = function(e) NA_real_),
			.groups = "drop"
		)

	# Annotate p-values
	pval_df <- pval_df %>%
		mutate(pval_label = ifelse(!is.na(pval), paste0("p = ", signif(pval, 2)), NA),
			y = 1)

	# Plot
	pos_dodge <- 0.9
	pos_jitter <- position_jitterdodge(jitter.width = 0.2, dodge.width = pos_dodge)
    
	p <- ggplot(res, aes(x = Dataset, y = pred_value)) + 
		geom_boxplot(aes(color = category), fill = NA, width = 0.7, outlier.shape = NA, position = position_dodge(pos_dodge)) +
		geom_jitter(aes(color = category), size = 1.5, shape = 16, alpha = 0.5, position = pos_jitter) +
		geom_hline(yintercept = e2g_threshold, linewidth = 0.5, linetype = "dashed", color = "#c5cad7") +
		geom_text(data = smry, aes(y = 0.95, group = category, label = n_label), size = 4, color = "black", position = position_dodge(pos_dodge)) +
        geom_text(data = smry, aes(y = 1, group = category, label = recall_label), size = 4, color = "black", position = position_dodge(pos_dodge)) +
		#geom_text(data = pval_df, aes(x = Dataset, y = y, label = pval_label), inherit.aes = FALSE, size = 3.5, vjust = 0) +
		scale_color_manual(values = cp) + scale_size_identity() +
		labs(x = "Dataset", y = "ENCODE-rE2G score", color = label) + 
		theme_classic() +
		theme(axis.text = element_text(size = 8, color = "#000000"), axis.title = element_text(size = 9),
			axis.ticks = element_line(color = "#000000"), axis.ticks.x = element_blank(),
			legend.position = "right")

    if (scale_log10) {
        p <- p + scale_y_log10()
    }

	w <- ifelse(length(order_levels) < 4, 4.5, 5.5)
    log_addn <- ifelse(scale_log10, ".log10", "")
	out_file_bp <- paste0(out_prefix, group_var, "_by_score_", significance_var, log_addn, ".pdf")
    out_file_pval <- paste0(out_prefix, group_var, "_by_score_", significance_var, ".pvalues.tsv")
    out_file_smry <- paste0(out_prefix, group_var, "_by_score_", significance_var, ".tsv")

    fwrite(pval_df, out_file_pval, sep = "\t", quote = FALSE, col.names = TRUE)
    fwrite(smry, out_file_smry, sep = "\t", quote = FALSE, col.names = TRUE)
	ggsave(out_file_bp, p, height = 5, width = w)
}

### --- ESTIMATE UNDETECTED VERSUS DETECTED POSITIVE SIGNIFICANT PAIRS --- 
plot_positives_by_effect_size <- function(crispr, es_bins, bin_min, bin_max, significance_var, out_prefix) {
    # stacked barplot with x-axis = effect size bins, y axis = # positive DE-G pairs
    # categories for bars: "statistically significant" vs "not statistically significant"
    # second category estimated based on rate of positives and sum of power for all tested pairs for low end of effect size
    
    ## params
    power_cols <- paste0("power_at_effect_size_", bin_min); names(power_cols) <- es_bins
    category_key <- c("Statistically significant" = "#435369", "Not detected as significant" = "#c5cad7")
    dataset_key <- c(all_DC_TAP = "DC-TAP", Gasperini2019 = "Gasperini et al.")

    ## calculate total power and pairs per dataset
    bin_sumTestedPower <- crispr %>%
        mutate(!!sym(power_cols[1]) := 0) %>%
        group_by(Dataset) %>%
        mutate(nTested_total = n()) %>% 
        group_by(Dataset, nTested_total) %>%
        summarize(across(all_of(power_cols), ~sum(.x, na.rm = TRUE)), .groups = "drop") %>%
        pivot_longer(cols = all_of(names(power_cols)),
            names_to = "EffectSize_bin",
            values_to = "sum_testedPower")

    ## summarize results
    res <- crispr %>% 
        mutate(abs_EffectSize = abs(EffectSize),
            detectedPositive = !!sym(significance_var),
            EffectSize_bin = cut(abs_EffectSize, breaks = c(bin_min[1], bin_max), labels = es_bins,
                right = TRUE, include.lowest = TRUE),
            EffectSize_bin = as.character(EffectSize_bin)) %>%
        group_by(Dataset, EffectSize_bin) %>% 
        summarize(nDetectedPositive = sum(detectedPositive),
            nTested_bin = n(),
            .groups = "drop") %>%
        left_join(bin_sumTestedPower, by = c("Dataset", "EffectSize_bin")) %>% 
        mutate(positiveRate = nDetectedPositive / sum_testedPower,
            nTotalPositive = positiveRate * nTested_total,
            nUndetectedPositive = ifelse(EffectSize_bin != es_bins[1], nTotalPositive - nDetectedPositive, 0),
            EffectSize_bin = factor(EffectSize_bin, levels = es_bins, ordered = TRUE)) %>%
        ungroup() %>%
        arrange(Dataset, EffectSize_bin)

    fwrite(res, paste0(out_prefix, "positives_by_effect_size.", significance_var, ".tsv"), sep = "\t", quote = FALSE)

    ## reformat for plotting
    res_long <- res %>%
        select(Dataset, EffectSize_bin, nDetectedPositive, nUndetectedPositive) %>%
    pivot_longer(cols = c(nDetectedPositive, nUndetectedPositive),
        names_to = "DetectionCategory", values_to = "n") %>%
    mutate(DetectionCategory = ifelse(DetectionCategory == "nDetectedPositive", names(category_key)[1], names(category_key)[2]),
        DetectionCategory = factor(DetectionCategory, levels = rev(names(category_key)), ordered = TRUE),
        DatasetLabel = dataset_key[Dataset],
        DatasetLabel = factor(DatasetLabel, levels = unique(dataset_key), ordered = TRUE))

    # plot
    p1 <- ggplot(res_long, aes(x = EffectSize_bin, y = n, fill = DetectionCategory)) +
        geom_bar(stat = "identity") +
        facet_wrap(~DatasetLabel, ncol = 2, scales = "free_y") +
        scale_fill_manual(values = category_key, name = NULL) +
        labs(x = "", y = "# of positive enhancer-gene pairs", fill = "Detected?") +
        theme_classic() + theme(
            strip.background = element_blank(), panel.grid = element_blank(), strip.text = element_text(size = 12, color = "#000000"),
            axis.text = element_text(size = 10, color = "#000000"), axis.text.x = element_blank(), #axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 10), axis.ticks = element_line(color = "#000000"), 
            legend.position = "top")

    res <- res %>% 
        mutate(DatasetLabel = dataset_key[Dataset],
            DatasetLabel = factor(DatasetLabel, levels = unique(dataset_key), ordered = TRUE),
            positiveRate = ifelse(sum_testedPower == 0, 0, positiveRate))

    p2 <- ggplot(res, aes(x = EffectSize_bin, y = positiveRate)) +
        geom_bar(stat = "identity", fill = "#1c2a43") +
        facet_wrap(~DatasetLabel, ncol = 2, scales = "fixed") +
        labs(x = "Effect size bin", y = "Positive rate\n(Detected positives / sum (power of tested pairs))") +
        theme_classic() + theme(
            strip.background = element_blank(), panel.grid = element_blank(), strip.text = element_blank(),
            axis.text = element_text(size = 10, color = "#000000"), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 10), axis.ticks = element_line(color = "#000000"),
            legend.position = "none")

    gr <- plot_grid(p1, p2, nrow = 2, align = "hv", rel_heights = c(1, 1))
    ggsave2(paste0(out_prefix, "positives_by_effect_size_and_rate.", significance_var, ".pdf"), gr, height = 6, width = 7)
}

### --- DATA I/O ---
# read in crispr data processed by dc-tap pipeline and annotated with chrom-annotate pipeline
read_sceptre_crispr_data <- function(crispr_path, direct_effect_path, filter_to_random_DEG = TRUE, annot = TRUE) {
    de <- fread(direct_effect_path) %>% 
        select(element_gene_pair_identifier_hg38, cell_type,
            direct_vs_indirect_negative, direct_vs_indirect_positive, distance_category, old_element_category = element_category) %>%
        distinct()

    res <- fread(crispr_path, sep = "\t") %>%
        mutate(cell_type = ifelse(cell_type == "WTC11_for_WTC11", "WTC11", cell_type)) %>%
        left_join(de, by = c("cell_type", "element_gene_pair_identifier_hg38")) %>% 
        filter(!is.na(direct_vs_indirect_negative))
    print(nrow(res))


    if (length(unique(res$cell_type)) == 2) {
        # DC-TAP
            temp_init <- res %>% filter(Random_DistalElement_Gene == TRUE) %>% filter(significant | power_at_effect_size_15 >= 0.8)
            print(nrow(temp_init))

        if (filter_to_random_DEG) {
            res <- res %>% 
                filter(Random_DistalElement_Gene == TRUE) %>% 
                mutate(Dataset = ifelse(cell_type == "K562", "K562_DC_TAP", "WTC11_DC_TAP"))
        }
        res <- res %>% mutate(data_category = "DC_TAP_SCEPTRE") 
    } else {
        # Gasperini
        if (filter_to_random_DEG) {
            res <- res %>% 
                filter(DistalElement_Gene == TRUE) 
        }
        res <- res %>% 
            mutate(Dataset = "Gasperini2019",
                data_category = "Gasperini_SCEPTRE")
    }

    # format columns
    res <- res %>% 
        select(chr = targeting_chr_hg38, start = targeting_start_hg38, end = targeting_end_hg38, cell_type,
            measuredGeneSymbol = gene_symbol, distanceToTSS = distance_to_gencode_gene_TSS,
            direct_vs_indirect_negative, direct_vs_indirect_positive,
            distance_category, old_element_category, #any_of(c("ubiq_category")),
            pct_change_effect_size, Significant = significant, Dataset, starts_with("power_at"),
            ends_with("peak_overlap"), ends_with(".RPM"), ends_with(".RPM.expandedRegion")) %>% 
        mutate(
            elementName = paste0(cell_type, "|", chr, ":", start, "-", end),
            WellPowered = (Significant | power_at_effect_size_15 >= 0.8), 
            EffectSize = pct_change_effect_size / 100,
            Regulated = (Significant & EffectSize < 0)) %>% 
        distinct() %>% 
        select(-pct_change_effect_size)

    if (annot) {
        res <- res %>% mutate(CTCF.H3K27ac.ratio = (CTCF.RPM) / (H3K27ac.RPM + 0.001))
    }


    return(res)
}

# read in genome-wide element lists annotated with chrom-annotate pipeline
read_enh_lists <- function(file_path, remove_promoters) {
    enh <- fread(file_path, sep = "\t") 
    if (remove_promoters) {
        enh <- filter(enh, class != "promoter")
    }

    return(enh)
}