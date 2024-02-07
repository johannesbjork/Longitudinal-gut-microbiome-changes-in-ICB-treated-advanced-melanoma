# Function for computing, from the fitted model, the marginal average for PFS>=12 (R) vs PFS<12 (NR) averaging across all levels of W1 (therapy regimen), W2 (colitis) and W3 (PPI-use), and plotting the posteriors whose "siglevel" CI does not cover 0  
# For the underlying calculations, please see the Supplementary Methods in Björk et al. (2024)
plot_RvsNR <- function(pibble_fit=pibble_fit, siglevel=siglevel) {
  
  focus_vars <- c("visit","PFS12yes","combiIOyes","colitisyes","ppiyes")

  avg_diffs <- as.data.frame.table(pibble_fit$Lambda) %>%
    pivot_wider(id_cols=c(Var1,Var3), names_from=Var2, values_from=Freq) %>%
    select(feature=Var1, b0=`(Intercept)`, contains(focus_vars)) %>% 
    mutate(feature=str_remove_all(string = feature, pattern = "f__"),
           feature=str_remove_all(string = feature, pattern = "s__"),
           feature=str_remove_all(string = feature, pattern = "t__"))

    group_by(feature) %>% 
    
    mutate(
      
      # Slopes for R
      
      slopes_NR =
        !!rlang::sym(focus_vars[1]) +
        0.5*(!!rlang::sym(paste0(focus_vars[1],":",focus_vars[3])) +
               !!rlang::sym(paste0(focus_vars[1],":",focus_vars[4])) + 
               !!rlang::sym(paste0(focus_vars[1],":",focus_vars[5]))),
      
      slopes_R =
        !!rlang::sym(focus_vars[1]) +
        !!rlang::sym(paste0(focus_vars[1],":",focus_vars[2])) +
        0.5*(!!rlang::sym(paste0(focus_vars[1],":",focus_vars[3])) +
               !!rlang::sym(paste0(focus_vars[1],":",focus_vars[4])) + 
               !!rlang::sym(paste0(focus_vars[1],":",focus_vars[5])) +
               
               !!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[3])) +
               !!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[4])) +       
               !!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[5]))),
      
      # R vs NR averaging across all W scenarios
      avgdiff_RvsNR_x_t0 =
        !!rlang::sym(focus_vars[2]) +
        
        0.5*(!!rlang::sym(paste0(focus_vars[2],":",focus_vars[3])) +
               !!rlang::sym(paste0(focus_vars[2],":",focus_vars[4])) +
               !!rlang::sym(paste0(focus_vars[2],":",focus_vars[5]))) +
        
        0*(!!rlang::sym(paste0(focus_vars[1],":",focus_vars[2]))) +
        
        0.5*(
          0*(!!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[3])) +
               !!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[4])) +
               !!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[5])))
        ),
      
      avgdiff_RvsNR_x_t1 =
        !!rlang::sym(focus_vars[2]) +
        
        0.5*(!!rlang::sym(paste0(focus_vars[2],":",focus_vars[3])) +
               !!rlang::sym(paste0(focus_vars[2],":",focus_vars[4])) +
               !!rlang::sym(paste0(focus_vars[2],":",focus_vars[5]))) +
        
        1*(!!rlang::sym(paste0(focus_vars[1],":",focus_vars[2]))) +
        
        0.5*(
          1*(!!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[3])) +
               !!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[4])) +
               !!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[5])))
        ),
      
      avgdiff_RvsNR_x_t2 =
        !!rlang::sym(focus_vars[2]) +
        
        0.5*(!!rlang::sym(paste0(focus_vars[2],":",focus_vars[3])) +
               !!rlang::sym(paste0(focus_vars[2],":",focus_vars[4])) +
               !!rlang::sym(paste0(focus_vars[2],":",focus_vars[5]))) +
        
        2*(!!rlang::sym(paste0(focus_vars[1],":",focus_vars[2]))) +
        
        0.5*(
          2*(!!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[3])) +
               !!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[4])) +
               !!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[5])))
        ),
      
      avgdiff_RvsNR_x_t3 =
        !!rlang::sym(focus_vars[2]) +
        0.5*(!!rlang::sym(paste0(focus_vars[2],":",focus_vars[3])) +
               !!rlang::sym(paste0(focus_vars[2],":",focus_vars[4])) +
               !!rlang::sym(paste0(focus_vars[2],":",focus_vars[5]))) +
        
        3*(!!rlang::sym(paste0(focus_vars[1],":",focus_vars[2]))) +
        
        0.5*(
          3*(!!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[3])) +
               !!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[4])) +
               !!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[5])))
        )
    ) %>%
    
    select(feature, slopes_R, slopes_NR, avgdiff_RvsNR_x_t0, avgdiff_RvsNR_x_t1, avgdiff_RvsNR_x_t2, avgdiff_RvsNR_x_t3)
  
  ranks_ls <- vector("list", 4)
  names(ranks_ls) <- c("avgdiff_RvsNR_x_t0","avgdiff_RvsNR_x_t1","avgdiff_RvsNR_x_t2","avgdiff_RvsNR_x_t3")
  
  for(i in names(ranks_ls)) {
    print(i)
    
    ranks_ls[[i]]$sig_increasing <-
      avg_diffs %>% 
      select(feature, !!rlang::sym(i)) %>% 
      group_by(feature) %>%
      ggdist::median_qi(!!rlang::sym(i), .width=c(0, 0.5, 0.75, 0.90, 0.95, 0.97, 1)) %>% 
      mutate(feature=reorder(factor(feature), !!rlang::sym(i))) %>% 
      pivot_wider(id_cols=feature, names_from=.width, values_from=.lower) %>%
      select(feature, p0=`0`, p50=`0.5`, p75=`0.75`, p90=`0.9`, p95=`0.95`, p97=`0.97`, p1=`1`) %>%
      filter(!! rlang::sym(siglevel) > 0) %>%
      mutate(feature=factor(feature)) %>% 
      select(feature) %>% 
      pull() %>% 
      levels() 
    
    ranks_ls[[i]]$sig_decreasing <-
      avg_diffs %>% 
      select(feature, !!rlang::sym(i)) %>% 
      group_by(feature) %>%
      ggdist::median_qi(!!rlang::sym(i), .width=c(0, 0.5, 0.75, 0.90, 0.95, 0.97, 1)) %>% 
      mutate(feature=reorder(factor(feature), !!rlang::sym(i))) %>% 
      pivot_wider(id_cols=feature, names_from=.width, values_from=.upper) %>%
      select(feature, p0=`0`, p50=`0.5`, p75=`0.75`, p90=`0.9`, p95=`0.95`, p97=`0.97`, p1=`1`) %>%
      filter(!! rlang::sym(siglevel) < 0) %>%
      mutate(feature=factor(feature)) %>% 
      select(feature) %>% 
      pull() %>% 
      levels()
  }
  
  all_taxa <- c(ranks_ls$avgdiff_RvsNR_x_t0$sig_increasing, ranks_ls$avgdiff_RvsNR_x_t0$sig_decreasing,
                ranks_ls$avgdiff_RvsNR_x_t1$sig_increasing, ranks_ls$avgdiff_RvsNR_x_t1$sig_decreasing,
                ranks_ls$avgdiff_RvsNR_x_t2$sig_increasing, ranks_ls$avgdiff_RvsNR_x_t2$sig_decreasing,
                ranks_ls$avgdiff_RvsNR_x_t3$sig_increasing, ranks_ls$avgdiff_RvsNR_x_t3$sig_decreasing)
  
  ranks_ls2 <- vector("list", 2)
  names(ranks_ls2) <- c("slopes_R","slopes_NR")
  
  for(i in names(ranks_ls2)) {
    print(i)
    
    ranks_ls2[[i]]$sig_increasing <-
      avg_diffs %>% 
      select(feature, !!rlang::sym(i)) %>% 
      group_by(feature) %>%
      ggdist::median_qi(!!rlang::sym(i), .width=c(0, 0.5, 0.75, 0.90, 0.95, 0.97, 1)) %>% 
      mutate(feature=reorder(factor(feature), !!rlang::sym(i))) %>% 
      pivot_wider(id_cols=feature, names_from=.width, values_from=.lower) %>%
      select(feature, p0=`0`, p50=`0.5`, p75=`0.75`, p90=`0.9`, p95=`0.95`, p97=`0.97`, p1=`1`) %>%
      
      #filter(!! rlang::sym(siglevel) > 0) %>% # this is to find non-zero slopes
      
      filter(feature %in% all_taxa) %>% # this is to retain differentially abundant taxa at the specified siglevel 
      
      mutate(feature=factor(feature)) %>% 
      select(feature) %>% 
      pull() %>% 
      levels() 
    
    ranks_ls2[[i]]$sig_decreasing <-
      avg_diffs %>% 
      select(feature, !!rlang::sym(i)) %>% 
      group_by(feature) %>%
      ggdist::median_qi(!!rlang::sym(i), .width=c(0, 0.5, 0.75, 0.90, 0.95, 0.97, 1)) %>% 
      mutate(feature=reorder(factor(feature), !!rlang::sym(i))) %>% 
      pivot_wider(id_cols=feature, names_from=.width, values_from=.upper) %>%
      select(feature, p0=`0`, p50=`0.5`, p75=`0.75`, p90=`0.9`, p95=`0.95`, p97=`0.97`, p1=`1`) %>%
      
      #filter(!! rlang::sym(siglevel) < 0) %>% # this is to find non-zero slopes
      
      filter(feature %in% all_taxa) %>% # this is to retain differentially abundant taxa at the specified siglevel 
      
      mutate(feature=factor(feature)) %>% 
      select(feature) %>% 
      pull() %>% 
      levels()
  }
  
  ## Color labeling features shared between different visits
  
  # Features shared between t2 and t1
  shared_t1_t2 <- Reduce(intersect, list(c(ranks_ls$avgdiff_RvsNR_x_t0$sig_increasing, ranks_ls$avgdiff_RvsNR_x_t0$sig_decreasing),
                                         c(ranks_ls$avgdiff_RvsNR_x_t1$sig_increasing, ranks_ls$avgdiff_RvsNR_x_t1$sig_decreasing)))
  
  # Features shared between t3 and t2
  shared_t2_t3 <- Reduce(intersect, list(c(ranks_ls$avgdiff_RvsNR_x_t1$sig_increasing, ranks_ls$avgdiff_RvsNR_x_t1$sig_decreasing),
                                         c(ranks_ls$avgdiff_RvsNR_x_t2$sig_increasing, ranks_ls$avgdiff_RvsNR_x_t2$sig_decreasing)))
  
  # Features shared between t4 and t3
  shared_t3_t4 <- Reduce(intersect, list(c(ranks_ls$avgdiff_RvsNR_x_t2$sig_increasing, ranks_ls$avgdiff_RvsNR_x_t2$sig_decreasing),
                                         c(ranks_ls$avgdiff_RvsNR_x_t3$sig_increasing, ranks_ls$avgdiff_RvsNR_x_t3$sig_decreasing)))
  
  # Features shared between t4, t3, t2
  shared_t2_t3_t4 <- Reduce(intersect, list(
    c(ranks_ls$avgdiff_RvsNR_x_t3$sig_increasing, ranks_ls$avgdiff_RvsNR_x_t3$sig_decreasing),
    c(ranks_ls$avgdiff_RvsNR_x_t2$sig_increasing, ranks_ls$avgdiff_RvsNR_x_t2$sig_decreasing),
    c(ranks_ls$avgdiff_RvsNR_x_t1$sig_increasing, ranks_ls$avgdiff_RvsNR_x_t1$sig_decreasing)))
  
  # Remember that t0 is visit 1 and t1 is visit 2 and so forth...
  y_labels_t2 <- setNames(rep("black",length(c(ranks_ls$avgdiff_RvsNR_x_t1$sig_increasing, ranks_ls$avgdiff_RvsNR_x_t1$sig_decreasing))),
                          rev(c(rev(ranks_ls$avgdiff_RvsNR_x_t1$sig_increasing), rev(ranks_ls$avgdiff_RvsNR_x_t1$sig_decreasing))))
  
  y_labels_t3 <- setNames(rep("black",length(c(ranks_ls$avgdiff_RvsNR_x_t2$sig_increasing, ranks_ls$avgdiff_RvsNR_x_t2$sig_decreasing))),
                          rev(c(rev(ranks_ls$avgdiff_RvsNR_x_t2$sig_increasing), rev(ranks_ls$avgdiff_RvsNR_x_t2$sig_decreasing))))
  
  y_labels_t4 <- setNames(rep("black",length(c(ranks_ls$avgdiff_RvsNR_x_t3$sig_increasing, ranks_ls$avgdiff_RvsNR_x_t3$sig_decreasing))),
                          rev(c(rev(ranks_ls$avgdiff_RvsNR_x_t3$sig_increasing), rev(ranks_ls$avgdiff_RvsNR_x_t3$sig_decreasing))))
  
  y_labels_t2[shared_t1_t2] <- "#74c476"
  y_labels_t3[shared_t2_t3] <- "#74c476"
  y_labels_t4[setdiff(shared_t3_t4, shared_t2_t3_t4)] <- "#74c476"
  y_labels_t4[shared_t2_t3_t4] <- "#005a32"
  
  plot_ls <- vector("list", 6)
  names(plot_ls) <- c("slopes_R","slopes_NR","avgdiff_RvsNR_x_t0","avgdiff_RvsNR_x_t1","avgdiff_RvsNR_x_t2","avgdiff_RvsNR_x_t3")
  
  for(p in names(plot_ls)) {
    print(p)
    
    if (p=="avgdiff_RvsNR_x_t0") {
      
      gp <- avg_diffs %>% 
        select(feature, !!rlang::sym(p)) %>% 
        group_by(feature) %>%
        ggdist::median_qi(!!rlang::sym(p), .width=c(0.5, 0.75, 0.90, 0.95, 0.97)) %>% 
        mutate(feature=reorder(factor(feature), !!rlang::sym(p))) %>% 
        filter(feature %in% c(ranks_ls[[p]]$sig_increasing, ranks_ls[[p]]$sig_decreasing)) %>%
        
        ggplot(aes(y=reorder(factor(feature), !!rlang::sym(p)), x=!!rlang::sym(p), xmin=.lower, xmax=.upper)) +
        ggdist::geom_interval(aes(alpha=.width), color="orange3") +
        scale_alpha_continuous("Credible interval", range=c(.7, .15), breaks=c(0.5, 0.75, 0.90, 0.95, 0.97)) +
        geom_point() +
        theme(
          legend.key=element_rect(fill='white'),
          legend.text=element_text(size=10, color="black"),
          strip.background=element_blank(),
          strip.text=element_blank(),
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), 
          panel.background=element_rect(fill="white"),
          panel.border=element_rect(colour="black", fill=NA, size=1),
          axis.ticks.length.y=unit(0.25,"cm"), 
          axis.text.x=element_text(size=12, color="black"),
          axis.text.y=element_text(size=7, color="black")) +
        labs(x="Log-Ratio Value", y=NULL, title=paste0("T0: R vs NR")) + 
        geom_vline(xintercept=0, linetype="dashed", color="darkgray")
      
      plot_ls[[p]] <- gp
      
    } else if (p=="avgdiff_RvsNR_x_t1") {
      
      gp <- avg_diffs %>% 
        select(feature, !!rlang::sym(p)) %>% 
        group_by(feature) %>%
        ggdist::median_qi(!!rlang::sym(p), .width=c(0.5, 0.75, 0.90, 0.95, 0.97)) %>% 
        mutate(feature=reorder(factor(feature), !!rlang::sym(p))) %>% 
        filter(feature %in% c(ranks_ls[[p]]$sig_increasing, ranks_ls[[p]]$sig_decreasing)) %>%
        
        ggplot(aes(y=reorder(factor(feature), !!rlang::sym(p)), x=!!rlang::sym(p), xmin=.lower, xmax=.upper)) +
        ggdist::geom_interval(aes(alpha=.width), color="orange3") +
        scale_alpha_continuous("Credible interval", range=c(.7, .15), breaks=c(0.5, 0.75, 0.90, 0.95, 0.97)) +
        geom_point() +
        theme(
          legend.key=element_rect(fill='white'),
          legend.text=element_text(size=10, color="black"),
          strip.background=element_blank(),
          strip.text=element_blank(),
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), 
          panel.background=element_rect(fill="white"),
          panel.border=element_rect(colour="black", fill=NA, size=1),
          axis.ticks.length.y=unit(0.25,"cm"), 
          axis.text.x=element_text(size=12, color="black"),
          axis.text.y=ggtext::element_markdown(color=unname(y_labels_t2), size=7)) +
        labs(x="Log-Ratio Value", y=NULL, title=paste0("T1: R vs NR")) + 
        geom_vline(xintercept=0, linetype="dashed", color="darkgray")
      
      plot_ls[[p]] <- gp
      
    } else if (p=="avgdiff_RvsNR_x_t2") {
      
      gp <- avg_diffs %>% 
        select(feature, !!rlang::sym(p)) %>% 
        group_by(feature) %>%
        ggdist::median_qi(!!rlang::sym(p), .width=c(0.5, 0.75, 0.90, 0.95, 0.97)) %>% 
        mutate(feature=reorder(factor(feature), !!rlang::sym(p))) %>% 
        filter(feature %in% c(ranks_ls[[p]]$sig_increasing, ranks_ls[[p]]$sig_decreasing)) %>%
        
        ggplot(aes(y=reorder(factor(feature), !!rlang::sym(p)), x=!!rlang::sym(p), xmin=.lower, xmax=.upper)) +
        ggdist::geom_interval(aes(alpha=.width), color="orange3") +
        scale_alpha_continuous("Credible interval", range=c(.7, .15), breaks=c(0.5, 0.75, 0.90, 0.95, 0.97)) +
        geom_point() +
        theme(
          legend.key=element_rect(fill='white'),
          legend.text=element_text(size=10, color="black"),
          strip.background=element_blank(),
          strip.text=element_blank(),
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), 
          panel.background=element_rect(fill="white"),
          panel.border=element_rect(colour="black", fill=NA, size=1),
          axis.ticks.length.y=unit(0.25,"cm"), 
          axis.text.x=element_text(size=12, color="black"),
          axis.text.y=ggtext::element_markdown(color=unname(y_labels_t3), size=7)) +
        labs(x="Log-Ratio Value", y=NULL, title=paste0("T2: R vs NR")) + 
        geom_vline(xintercept=0, linetype="dashed", color="darkgray")
      
      plot_ls[[p]] <- gp
      
    } else if (p=="avgdiff_RvsNR_x_t3") { 
      
      gp <- avg_diffs %>% 
        select(feature, !!rlang::sym(p)) %>% 
        group_by(feature) %>%
        ggdist::median_qi(!!rlang::sym(p), .width=c(0.5, 0.75, 0.90, 0.95, 0.97)) %>% 
        mutate(feature=reorder(factor(feature), !!rlang::sym(p))) %>% 
        filter(feature %in% c(ranks_ls[[p]]$sig_increasing, ranks_ls[[p]]$sig_decreasing)) %>%
        
        ggplot(aes(y=reorder(factor(feature), !!rlang::sym(p)), x=!!rlang::sym(p), xmin=.lower, xmax=.upper)) +
        ggdist::geom_interval(aes(alpha=.width), color="orange3") +
        scale_alpha_continuous("Credible interval", range=c(.7, .15), breaks=c(0.5, 0.75, 0.90, 0.95, 0.97)) +
        geom_point() +
        theme(
          legend.key=element_rect(fill='white'),
          legend.text=element_text(size=10, color="black"),
          strip.background=element_blank(),
          strip.text=element_blank(),
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), 
          panel.background=element_rect(fill="white"),
          panel.border=element_rect(colour="black", fill=NA, size=1),
          axis.ticks.length.y=unit(0.25,"cm"), 
          axis.text.x=element_text(size=12, color="black"),
          axis.text.y=ggtext::element_markdown(color=unname(y_labels_t4), size=7)) +
        labs(x="Log-Ratio Value", y=NULL, title=paste0("T3: R vs NR")) + 
        geom_vline(xintercept=0, linetype="dashed", color="darkgray")
      
      plot_ls[[p]] <- gp
      
    } else if (p=="slopes_R") { 
      
      gp <- avg_diffs %>% 
        select(feature, !!rlang::sym(p)) %>% 
        group_by(feature) %>%
        ggdist::median_qi(!!rlang::sym(p), .width=c(0.5, 0.75, 0.90, 0.95, 0.97)) %>% 
        mutate(feature=reorder(factor(feature), !!rlang::sym(p))) %>% 
        filter(feature %in% c(ranks_ls2[[p]]$sig_increasing, ranks_ls2[[p]]$sig_decreasing)) %>%
        
        ggplot(aes(y=reorder(factor(feature), !!rlang::sym(p)), x=!!rlang::sym(p), xmin=.lower, xmax=.upper)) +
        ggdist::geom_interval(aes(alpha=.width), color="orange3") +
        scale_alpha_continuous("Credible interval", range=c(.7, .15), breaks=c(0.5, 0.75, 0.90, 0.95, 0.97)) +
        geom_point() +
        theme(
          legend.key=element_rect(fill='white'),
          legend.text=element_text(size=10, color="black"),
          strip.background=element_blank(),
          strip.text=element_blank(),
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), 
          panel.background=element_rect(fill="white"),
          panel.border=element_rect(colour="black", fill=NA, size=1),
          axis.ticks.length.y=unit(0.25,"cm"), 
          axis.text.x=element_text(size=12, color="black"),
          axis.text.y=ggtext::element_markdown(color=unname(y_labels_t4), size=7)) +
        labs(x="Log-Ratio Value", y=NULL, title=paste0("Slope R")) + 
        geom_vline(xintercept=0, linetype="dashed", color="darkgray")
      
      plot_ls[[p]] <- gp
      
    } else if (p=="slopes_NR") { 
      
      gp <- avg_diffs %>% 
        select(feature, !!rlang::sym(p)) %>% 
        group_by(feature) %>%
        ggdist::median_qi(!!rlang::sym(p), .width=c(0.5, 0.75, 0.90, 0.95, 0.97)) %>% 
        mutate(feature=reorder(factor(feature), !!rlang::sym(p))) %>% 
        filter(feature %in% c(ranks_ls2[[p]]$sig_increasing, ranks_ls2[[p]]$sig_decreasing)) %>%
        
        ggplot(aes(y=reorder(factor(feature), !!rlang::sym(p)), x=!!rlang::sym(p), xmin=.lower, xmax=.upper)) +
        ggdist::geom_interval(aes(alpha=.width), color="orange3") +
        scale_alpha_continuous("Credible interval", range=c(.7, .15), breaks=c(0.5, 0.75, 0.90, 0.95, 0.97)) +
        geom_point() +
        theme(
          legend.key=element_rect(fill='white'),
          legend.text=element_text(size=10, color="black"),
          strip.background=element_blank(),
          strip.text=element_blank(),
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), 
          panel.background=element_rect(fill="white"),
          panel.border=element_rect(colour="black", fill=NA, size=1),
          axis.ticks.length.y=unit(0.25,"cm"), 
          axis.text.x=element_text(size=12, color="black"),
          axis.text.y=ggtext::element_markdown(color=unname(y_labels_t4), size=7)) +
        labs(x="Log-Ratio Value", y=NULL, title=paste0("Slope NR")) + 
        geom_vline(xintercept=0, linetype="dashed", color="darkgray")
      
      plot_ls[[p]] <- gp
    }
  }
  return(plot_ls)
}

# Then we run this function for the "siglevel" that we want. 
avg_RvsNR_p90 <- plot_RvsNR(pibble_fit=fit_species_clr, siglevel="p90")
