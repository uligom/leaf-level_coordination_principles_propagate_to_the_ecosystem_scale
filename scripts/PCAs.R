PCAs = function(figs = c(1:3), boots = 499) {
  #### PRINCIPAL COMPONENT ANALYSIS
  
  ### Author: Ulisse Gomarasca, Mirco Migliavacca, Talie Musavi
  ### Script options -------------------------------------------------------------
  # Measure script run time
  library(tictoc)
  tic()
  
  savedata <- T # save output
  
  
  
  ### Utilities ------------------------------------------------------------------
  ## Packages:
  packages <- c("ggplot2", "reshape2", "FactoMineR", "factoextra", "gridExtra",
                "ade4", "tidyr", "dplyr", "readr", "ggrepel", "ggridges", "cowplot", "modelr", "stringr")
  
  # Install packages not yet installed
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
  }
  
  invisible(lapply(packages, library, character.only = TRUE))
  
  ## Other:
  source("scripts/utils.R")
  source("scripts/units.R") # plot labels with variable names and units
  
  
  
  ### Import data ----------------------------------------------------------------
  df <- read.table("data/input_data.csv", header = T, sep = ",")
  
  
  
  ### Analysis options -----------------------------------------------------------
  ## Set number of runs (bootstrapping)
  tt <- boots # number of repetitions for Dray test and error bars. WARNING: time consuming
  
  
  
  ### Analysis -----------------------------------------------------------------
  ## Initialize site list
  site_list <- tibble() # initialize cumulative list of sites used in all analyses
  
  
  
  for (qq in figs) {
    ### Set Scenario for analysis ----------------------------------------------
    n_pcs <- NA # initialize number of PCs to retain
    ## 1) LEAF ECONOMICS SPECTRUM ####
    if (qq == 1){
      ### Saving options ----
      figure_name <- "Figure_1"
      
      
      ### Select variables and subset sites ------------------------------------
      ## Selecting EFPs
      EFPs_codes_4_PCA <- c("SITE_ID", "IGBP", "GPPsat", "wLL", "wNmass", "wLMA", "RECOmax")
      
      ## Filter
      df_subset <- df %>%
        dplyr::select(all_of(EFPs_codes_4_PCA)) %>%
        # dplyr::filter(IGBP != "ENF") %>% # OPTION: test bias of ENF by removing them
        tidyr::drop_na() # exclude NA (TO DO, only for most limiting variables? would need multiple imputation)
    }
    
    
    
    ## 11) LEAF ECONOMICS SPECTRUM - wNarea ####
    if (qq == 11){
      ### Saving options ----
      figure_name <- "LES_wNarea"
      
      
      ### Select variables and subset sites ------------------------------------
      ## Selecting EFPs
      EFPs_codes_4_PCA <- c("SITE_ID", "IGBP", "GPPsat", "wLL", "wNarea", "RECOmax")
      
      ## Filter
      df_subset <- df %>%
        dplyr::select(all_of(EFPs_codes_4_PCA)) %>%
        # dplyr::filter(IGBP != "ENF") %>% # OPTION: test bias of ENF by removing them
        tidyr::drop_na() # exclude NA (TO DO, only for most limiting variables? would need multiple imputation)
    }
    
    
    
    ## 12) LEAF ECONOMICS SPECTRUM - GPP/LAI ####
    if (qq == 12){
      ### Saving options ----
      figure_name <- "LES_GPP_LAI"
      
      
      ### Select variables and subset sites --------------------------------------
      ## Selecting EFPs
      EFPs_codes_4_PCA <- c("SITE_ID", "IGBP", "GPP_LAI", "wLL", "wLMA", "wNmass", "RECOmax")
      
      ## Filter
      df_subset <- df %>%
        dplyr::select(all_of(EFPs_codes_4_PCA)) %>%
        # dplyr::filter(IGBP != "ENF") %>% # OPTION: test bias of ENF by removing them
        tidyr::drop_na() # exclude NA (TO DO, only for most limiting variables? would need multiple imputation)
    }
    
    
    
    ## 2) GLOBAL SPECTRUM ####
    if (qq == 2){
      ### Saving options ----
      figure_name <- "Figure_2"
      
      
      ### Select number of PCs for plotting ----
      n_pcs <- 3 # NB: manual selection of x components
      
      
      ### Select variables and subset sites --------------------------------------
      ## Selecting EFPs
      EFPs_codes_4_PCA <- c("SITE_ID", "IGBP",
                            "GPPsat",
                            "wLMA", "wNmass", "wSSD",
                            "Hc", "LAImax")
      
      ## Filter
      df_subset <- df %>%
        dplyr::select(all_of(EFPs_codes_4_PCA)) %>%
        # dplyr::filter(IGBP != "ENF") %>% # OPTION: test bias of ENF by removing them
        tidyr::drop_na() # exclude NA (TO DO, only for most limiting variables? would need multiple imputation)
    }
    
    
    ## 3) LEAST COST HYPOTHESIS ####
    if (qq == 3){
      ### Saving options ----
      figure_name <- "Figure_3"
      
      
      ### Select variables and subset sites --------------------------------------
      ## Selecting EFPs
      EFPs_codes_4_PCA <- c("SITE_ID", "IGBP", "PNUE", "Gsmax", "WUEt", "Hc", "EF", "Ta")
      
      ## Filter
      df_subset <- df %>%
        dplyr::select(all_of(EFPs_codes_4_PCA)) %>%
        # dplyr::filter(IGBP != "ENF") %>% # OPTION: test bias of ENF by removing them
        tidyr::drop_na() # exclude NA (TO DO, only for most limiting variables? would need multiple imputation)
    }
    
    
    ## 4) LEAST COST HYPOTHESIS - PNUEarea variant ####
    if (qq == 4){
      ### Variant name (defines other options) ----
      variant <- "PNUEarea"
      
      
      ### Saving options ----
      figure_name <- glue::glue("Supplementary_Figure_4")
      
      
      ### Select variables and subset sites --------------------------------------
      ## Selecting EFPs
      EFPs_codes_4_PCA <- c("SITE_ID", "IGBP", variant, "LAImax", "Gsmax", "WUEt", "Hc", "EF", "Ta")
      
      ## Filter
      df_subset <- df %>%
        dplyr::select(all_of(EFPs_codes_4_PCA)) %>%
        # dplyr::filter(IGBP != "ENF") %>% # OPTION: test bias of ENF by removing them
        tidyr::drop_na() # exclude NA (TO DO, only for most limiting variables? would need multiple imputation)
    }
    
    
    
    ## 5) LEAST COST HYPOTHESIS - APAR_N variant ####
    if (qq == 5){
      ### Variant name (defines other options) ----
      variant <- "APAR_N"
      
      
      ### Saving options ----
      figure_name <- glue::glue("Supplementary_Figure_5")
      
      
      ### Select variables and subset sites --------------------------------------
      ## Selecting EFPs
      EFPs_codes_4_PCA <- c("SITE_ID", "IGBP", variant, "LAImax", "Gsmax", "WUEt", "Hc", "EF", "Ta")
      
      ## Filter
      df_subset <- df %>%
        dplyr::select(all_of(EFPs_codes_4_PCA)) %>%
        # dplyr::filter(IGBP != "ENF") %>% # OPTION: test bias of ENF by removing them
        tidyr::drop_na() # exclude NA (TO DO, only for most limiting variables? would need multiple imputation)
    }
    
    
    
    ## 6) LEAST COST HYPOTHESIS - WUEet variant ####
    if (qq == 6){
      ### Variant name (defines other options) ----
      variant <- "WUEet"
      
      
      ### Saving options ----
      figure_name <- glue::glue("Supplementary_Figure_6")
      
      n_pcs <- 2
      
      ### Select variables and subset sites --------------------------------------
      ## Selecting EFPs
      EFPs_codes_4_PCA <- c("SITE_ID", "IGBP", "PNUE", "Gsmax", variant, "Hc", "EF", "Ta")
      
      ## Filter
      df_subset <- df %>%
        dplyr::select(all_of(EFPs_codes_4_PCA)) %>%
        # dplyr::filter(IGBP != "ENF") %>% # OPTION: test bias of ENF by removing them
        tidyr::drop_na() # exclude NA (TO DO, only for most limiting variables? would need multiple imputation)
    }
    
    
    ## 7) LEAST COST HYPOTHESIS - uWUE variant ####
    if (qq == 7){
      ### Variant name (defines other options) ----
      variant <- "uWUE"
      
      
      ### Saving options ----
      figure_name <- glue::glue("Supplementary_Figure_7")
      
      n_pcs <- 2
      
      ### Select variables and subset sites --------------------------------------
      ## Selecting EFPs
      EFPs_codes_4_PCA <- c("SITE_ID", "IGBP", "PNUE", "Gsmax", variant, "Hc", "EF", "Ta")
      
      ## Filter
      df_subset <- df %>%
        dplyr::select(all_of(EFPs_codes_4_PCA)) %>%
        # dplyr::filter(IGBP != "ENF") %>% # OPTION: test bias of ENF by removing them
        tidyr::drop_na() # exclude NA (TO DO, only for most limiting variables? would need multiple imputation)
    }
    
    
    
    ##############################################################################
    print(glue::glue("==========================================================
                  ==========================================================
                   Running analysis for '{figure_name}'....."))
    
    
    ### Note down used sites -----------------------------------------------------
    site_list <- bind_rows(site_list, tibble(df_subset %>% select(SITE_ID), analysis = figure_name)) %>% 
      unique() %>% 
      arrange(SITE_ID)
    
    
    
    ### Bootstrap ----------------------------------------------------------------
    set.seed(8675309) # set seed for reproducibility
    
    df_boot <- df_subset %>% modelr::bootstrap(n = tt) # bootstrap
    
    
    
    ### PCA ----------------------------------------------------------------------
    ## Extract bootstrapped data
    # Initialize
    N_PCS <- tibble() # number of PCs as retained PCs by Dray method if not otherwise specified
    pca_stats <- tibble() # variable contributions and loadings
    R2 <- c() # explained variance
    
    # Loop
    for (j in 1:tt) {
      ## Extract bootstrapped data for jth bootstrap
      dat <- df_boot %>% 
        slice(j) %>% # select jth row (with bootstrapped dataframe)
        pull(strap) %>% # extract list
        as.data.frame() # extract data
      
      ## Run PCA without multiple imputation
      pca_result <- FactoMineR::PCA(dat %>% dplyr::select(-SITE_ID, -IGBP), scale.unit = T, ncp = 10, graph = F)
      
      
      
      ### Check uncertainties and significance ----------------------------------
      ## Uncertainty PCA following Dray et al., 2008 ----
      # and the methodology described in Diaz et al., 2015
      tab <- na.omit(dat %>% dplyr::select(-SITE_ID, -IGBP))
      pca1 <- ade4::dudi.pca(tab, center = TRUE, scale = TRUE, scannf = FALSE, nf = 10)
      
      ## Dray test
      pc_tested <- ade4::testdim(pca1, nrepet = 999)
      
      
      
      ### Export bootstrapped PCA outputs ----------------------------------------
      N_PCS <- N_PCS %>% 
        bind_rows(tibble(strap = j, n_pcs = pc_tested$nb.cor)) # number of PCs for jth run
      
      pca_stats <- bind_rows(pca_stats,
                             pca_result$var$contrib %>% # add contributions
                               as_tibble(rownames = "var") %>%
                               pivot_longer(cols = !var, names_to = "PC", values_to = "contrib") %>% 
                               left_join(pca_result$var$coord %>% # add loadings
                                           as_tibble(rownames = "var") %>%
                                           pivot_longer(cols = !var, names_to = "PC", values_to = "loading"),
                                         by = c("var", "PC")
                               ) %>% 
                               mutate(PC = str_sub(PC, start = 5), # PC name to digit
                                      strap = j) # bootstrap run number
      ) # variable contributions and loadings
      
      R2 <- bind_rows(R2,
                      tibble(PC = pca_result[["eig"]] %>% rownames(),
                             exp_var = pca_result[["eig"]][,2],
                             strap = j) %>% 
                        mutate(PC = str_sub(PC, start = 6)) # PC name to digit
      )
      
    }
    
    
    
    ### Statistics of bootstrapped extracted PCA outputs -------------------------
    ## Output of Dray method ----
    # Number of PCs to retain
    N_PCS <- N_PCS %>%
      group_by(n_pcs) %>% 
      summarise(n_rep = n()) %>% # count of repeated values
      mutate(retained = n_rep / tt * 100) # percentage of runs
    # print(pc_tested)
    pc_ret <- N_PCS %>% filter(retained == max(retained))
    
    print(N_PCS)
    print(paste0("Number of statistical significant components according to Dray method (Dray et al., 2008) was ",
                 pc_ret[1,1], " in ", round(pc_ret[1,3], digits = 1), "% of runs."))
    
    if (is.na(n_pcs)) {
      n_pcs <- N_PCS %>% 
        filter(retained == max(retained)) %>% 
        select(n_pcs) %>% 
        unlist() %>% unname()
    } # set number of PCs as retained PCs by Dray method if not otherwise specified
    
    
    ## Contributions and loadings ----
    pca_stats <- pca_stats %>% 
      group_by(PC, var) %>% 
      summarize(#contrib_mean = mean(contrib),
                contrib_std = sd(contrib),
                # contrib_q25 = quantile(contrib, 0.25), contrib_q75 = quantile(contrib, 0.75),
                # loading_mean = mean(loading),
                loading_std = sd(loading),
                # loading_q25 = quantile(loading, 0.25), loading_q75 = quantile(loading, 0.75),
                .groups = "drop"
      )
    
    ## Modify PC names
    pca_stats <- pca_stats %>%
      mutate(PC_name = paste0("PC", PC))
    
    
    
    ## Explained variance ----
    R2 <- R2 %>% 
      group_by(PC) %>% 
      summarize(#R2_mean = mean(exp_var),
                R2_std = sd(exp_var),
                # R2_q25 = quantile(exp_var, 0.25), R2_q75 = quantile(exp_var, 0.75),
                .groups = "drop"
      )
    
    ## Add to stats
    pca_stats <- pca_stats %>% left_join(R2, by = "PC")
    
    
    
    ### Repeat PCA on main subset ------------------------------------------------
    pca_result <- FactoMineR::PCA(df_subset %>% dplyr::select(-SITE_ID, -IGBP), scale.unit = T, ncp = 9, graph = F) # repeat PCA on actual subset (not bootstrapped)
    
    ## Include actual values for non-bootstrapped dataset
    pca_stats <- pca_stats %>% 
      dplyr::left_join( # add R2 of analysis on actual data (not the mean of bootstrapping)
        tibble(PC = pca_result[["eig"]] %>% rownames(),
               R2 = pca_result[["eig"]][,2]
        ) %>% 
          mutate(PC = str_sub(PC, start = 6)), # extract PC numbers
        by = "PC"
      ) %>% 
      dplyr::left_join( # add contribution values of analysis on actual data (not the mean of bootstrapping)
        pca_result$var$contrib %>% # add contributions
          as_tibble(rownames = "var") %>%
          pivot_longer(cols = !var, names_to = "PC", values_to = "contrib") %>% 
          mutate(PC = str_sub(PC, start = 5)), # extract PC numbers
        by = c("PC", "var")
      ) %>%
      dplyr::left_join( # add loading values of analysis on actual data (not the mean of bootstrapping)
        pca_result$var$coord %>% # add contributions
          as_tibble(rownames = "var") %>%
          pivot_longer(cols = !var, names_to = "PC", values_to = "loading") %>% 
          mutate(PC = str_sub(PC, start = 5)), # extract PC numbers
        by = c("PC", "var")
      )
    
    ## Include percentage of retention for PCs (how often PC is retained over bootstraps)
    pca_stats <- pca_stats %>% 
      dplyr::left_join(N_PCS %>% dplyr::mutate(PC = n_pcs %>% as.character) %>% dplyr::select(PC, retained),
                       by = "PC"
      ) 
    
    ## Clean output
    pca_stats <- pca_stats %>% 
      dplyr::mutate(PC = as.integer(PC)) %>% 
      dplyr::relocate(PC_name, .after = PC) %>% 
      dplyr::relocate(R2, .before = R2_std) %>% 
      dplyr::relocate(contrib, .before = contrib_std) %>% 
      dplyr::relocate(loading, .before = loading_std) %>% 
      dplyr::arrange(PC)
    
    
    
    ### Plot ---------------------------------------------------------------------
    ## Plotting settings
    # Color palette:
    CatCol <- c(
      CSH = "#586158", DBF = "#C46B39", EBF = "#4DD8C0", ENF = "#3885AB", GRA = "#9C4DC4",
      MF = "#C4AA4D", OSH = "#443396", SAV = "#CC99CC", WET = "#88C44D", WSA = "#AB3232"
    )
    
    # Text settings:
    # if (n_pcs > 3) {x_angle <- 270; x_adjust <- 0.25} else {x_angle <- 0; x_adjust <- 0} # option to change orientation of x axis text
    x_angle <- 0; x_adjust <- 0
    
    
    ## A) biplot color coded by IGBP biome type ----
    p_biplot <- fviz_pca_biplot(pca_result,
                                axes = c(1, 2),
                                fill.ind = df_subset$IGBP, #"grey50",
                                col.ind = "white",
                                geom.ind = "point",
                                palette = CatCol,#'futurama',
                                label = "var",
                                col.var = "black",
                                labelsize = 8,
                                repel = TRUE,
                                pointshape = 21,
                                pointsize = 4,
                                alpha.ind = 0.67,
                                arrowsize = 1) +
      labs(title = "",
           x = "PC1",
           y = "PC2",
           fill = "IGBP") +
      guides(fill = guide_legend(title = "")) +
      theme(title = element_blank(),
            text = element_text(size = 20),
            axis.title = element_text(size = 28, face = "bold"),
            axis.text = element_text(size = 20),
            plot.margin = unit(c(0, 0, 0, 0), "cm"),
            # legend.position = "none"
            legend.text = element_text(size = 24),
            legend.key.height = unit(15, "mm"),
            legend.key.width = unit(20, "mm")
      ) +
      NULL
    
    
    
    ## B) explained variance ----
    p_r2 <- pca_stats %>%
      dplyr::select(PC_name, PC, R2, R2_std) %>%
      unique() %>%
      dplyr::mutate(PC = as.character(PC)) %>%
      ggplot(aes(x = PC_name, y = R2, group = 1)) + # x = PC -> only numbers on axis, x = PC_name -> can give problems with PC10 being ordered before PC2;
      # group 1 to avoid warning/error for some reason
      geom_line(linewidth = 2, color = "gray70") +
      geom_errorbar(aes(ymin = R2 - R2_std, ymax = R2 + R2_std), # ymin = R2_q25, ymax = R2_q75
                    color = "gray50", linewidth = 1.25, width = 0.2) + # standard deviation from bootstrapping
      geom_point(shape = 15, size = 4, color = "gray50") +
      geom_text(aes(x = PC_name, y = R2 + R2_std + 2, label = paste0(R2 %>% round(digits = 1), "%")),
                nudge_x = 0.25, size = 7) + # values
      labs(title = "", x = "", y = "Explained variance") +
      theme_classic() +
      theme(title = element_blank(),
            # text = element_text(size = 20),
            axis.title = element_text(size = 24, face = "bold"),
            axis.title.x = element_blank(),
            axis.text = element_text(size = 20),
            # axis.line.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.margin = unit(c(1, 0, 0, 2), "cm"),
            legend.position = "none"
      ) +
      NULL
    
    
    
    ## C) loadings ----
    p_load <- pca_stats %>%
      dplyr::filter(PC <= n_pcs[1]) %>%
      ggplot(aes(x = var, y = loading)) + #color = varSignificance, fill = varSignificance)) +
      facet_grid(. ~ PC_name, scales = "free_y") +
      geom_errorbar(aes(ymin = loading - loading_std, ymax = loading + loading_std), # loading_q25, ymax = loading_q75
                    linewidth = 0.75, width = 0.75) + # standard error = std from bootstrapping
      # geom_point(color = "gray50", shape = 15, size = 3) +
      geom_bar(stat = "identity", position = position_dodge(), fill = "#b07a4f", linewidth = 0.5) + #b07a4f, #9c6a5e, #643c3c
      geom_hline(yintercept = 0) +
      coord_flip() +  # flip axis for better visualization
      scale_y_continuous(breaks = waiver(), n.breaks = 4) + # change x axis (y since it is flipped)
      labs(y = "Loadings", x = "", title = "") +
      theme_classic() +
      theme(title = element_text(size = 20, face = "bold"),
            text = element_text(size = 20),
            axis.line.y = element_blank(),
            axis.title = element_text(size = 28),
            axis.text = element_text(size = 20),
            axis.text.x = element_text(angle = x_angle, vjust = x_adjust),
            axis.ticks.y = element_blank(),
            legend.position = "none",
            legend.title = element_text(size = 28),
            legend.text = element_text(size = 24),
            legend.key.height = unit(10, "mm"),
            legend.key.width = unit(10, "mm"),
            plot.margin = unit(c(0, 0, 0, 0), "cm"),
            strip.text = element_text(face = "bold", size = 28),
            strip.background = element_blank()
      ) +
      NULL
    
    
    
    ## D) contributions ----
    p_contr <- pca_stats %>%
      dplyr::filter(PC <= n_pcs[1]) %>% # remove additional PCs
      ggplot(aes(x = var, y = contrib)) + #color = varSignificance, fill = varSignificance)) +
      facet_grid(. ~ PC_name, scales = "free_y") +
      geom_errorbar(aes(ymin = contrib - contrib_std, ymax = contrib + contrib_std), # ymin = contrib_q25, ymax = contrib_q75
                    linewidth = 0.75, width = 0.75) + # standard error = standard deviation from bootstrapping
      # geom_point(color = "gray50", shape = 15, size = 3) +
      geom_bar(stat = "identity", position = position_dodge(), fill = "#4f85b0", linewidth = 0.5) + #4f85b0, #59918e, #3c6464
      geom_hline(yintercept = 0) +
      coord_flip() + # flip axis for better visualization
      scale_y_continuous(breaks = waiver(), n.breaks = 4) + # change x axis (y since it is flipped)
      labs(y = "Contribution [%]", x = "", title = "") +
      theme_classic() +
      theme(title = element_text(size = 20, face = "bold"),
            text = element_text(size = 20),
            axis.line.y = element_blank(),
            axis.title = element_text(size = 28),
            axis.text = element_text(size = 20),
            axis.text.x = element_text(angle = x_angle, vjust = x_adjust),
            axis.ticks.y = element_blank(),
            legend.position = "none",
            legend.title = element_text(size = 28),
            legend.text = element_text(size = 24),
            legend.key.height = unit(10, "mm"),
            legend.key.width = unit(10, "mm"),
            plot.margin = unit(c(0, 0, 0, 0), "cm"),
            strip.text = element_text(face = "bold", size = 28),
            strip.background = element_blank()
      ) +
      NULL
    
    
    
    ### Figure for the article ---------------------------------------------------
    fig_top <- cowplot::plot_grid( # top row of final figure
      p_biplot, p_r2,
      align = "h",
      nrow = 1,
      ncol = 2,
      rel_widths = c(1.618, 1), # golden ratio
      labels = c("a", "b"),
      label_size = 23
    )
    fig_bottom <- cowplot::plot_grid( # bottom row of final figure
      p_load, p_contr,
      align = "h",
      axis = "l",
      nrow = 1,
      ncol = 2,
      labels = c("c", "d"),
      label_size = 24
    )
    
    
    fig_paper <- cowplot::plot_grid( # final figure
      fig_top, fig_bottom,
      align = "v",
      nrow = 2,
      ncol = 1,
      rel_heights = c(1.618, 1) # golden ratio
    )
    # print(fig_paper)
    
    R2_N_label <- paste0("Retained axis (", n_pcs[1], ") explained ", round(sum(pca_result$eig[1:n_pcs[1],2]), 1), "% of variance (N = ", nrow(df_subset), ").")
    print(R2_N_label)
    
    
    
    ### Save ---------------------------------------------------------------------
    if (savedata == T){
      ## Final figure
      ggplot2::ggsave(filename = glue::glue("results/{figure_name}.jpg"), plot = fig_paper, device = "jpeg",
                      width = 508, height = 285.75, units = "mm", dpi = 600)
      
      ## PCA statistics
      write.table(pca_stats, file = glue::glue("results/PCA_stats_{figure_name}.csv"), row.names = F, sep = ";")
      
      ## Final data input for each analysis
      write.table(df_subset, file = glue::glue("results/PCA_input_{figure_name}.csv"), row.names = F, sep = ";")
    }
  } # end loop for different scenarios (analysis sections)
  
  
  
  ### Save site list -------------------------------------------------------------
  if (savedata == T ) {
    # Save list of used sites
    write.table(site_list %>% dplyr::select(SITE_ID) %>% unique(), file = glue::glue("results/site_list.csv"), row.names = F, sep = ";")
    
    if (figs %in% c(1:3)) {
      # Save IGBP stats
      site_igbp <- site_list %>%
        dplyr::filter(analysis %in% c("Figure_1", "Figure_2", "Figure_3")) %>%
        dplyr::left_join(df %>% select(SITE_ID, IGBP)) %>%
        dplyr::group_by(analysis) %>%
        dplyr::mutate(n_sites = n()) %>%
        dplyr::group_by(IGBP, analysis) %>%
        dplyr::summarise(site_fraction = n() * 100 / n_sites) %>%
        unique() %>%
        dplyr::mutate(analysis = case_when(
          analysis == "Figure_1" ~ "Leaf economics spectrum",
          analysis == "Figure_2" ~ "Global spectrum of plant form and function",
          analysis == "Figure_3" ~ "Least cost hypothesis"
          # analysis == "Supplementary_Figure_4" ~ "Least cost hypothesis - Supplementary Figure 4",
          # analysis == "Supplementary_Figure_5" ~ "Least cost hypothesis - Supplementary Figure 5",
          # analysis == "Supplementary_Figure_6" ~ "Least cost hypothesis - Supplementary Figure 6",
          # analysis == "Supplementary_Figure_7" ~ "Least cost hypothesis - Supplementary Figure 7"
        )
        ) %>%
        dplyr::arrange(IGBP, analysis) %>%
        tidyr::pivot_wider(names_from = analysis, values_from = site_fraction)
      
      write_csv(site_igbp, file = "results/Table_S7.csv")
    }
  }
  ### End ------------------------------------------------------------------------
  toc()
  
  
}