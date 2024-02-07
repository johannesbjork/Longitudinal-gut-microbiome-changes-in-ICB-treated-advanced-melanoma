# Function for computing and plotting the marginal average for PFS12>=12 (R) vs PFS<12 (NR) on combination therapy (W1=1), no colitis (W2=0) and no PPIs (W3=0)
plot_RvsNR_combitherapy <- function(pibble_fit=pibble_fit, focus_var=focus_var, siglevel=siglevel) {
  
  focus_vars <- str_split(focus_var,":")[[1]]
  focus_vars <- c(focus_vars[1],paste0(focus_vars[2:3], "yes"))
  print(focus_vars) 
  
  intercepts_and_slopes <- 
    as.data.frame.table(pibble_fit$Lambda) %>%
    pivot_wider(c(Var1,Var3), names_from=Var2, values_from=Freq) %>%
    select(feature=Var1, b0=`(Intercept)`, contains(focus_vars)) %>% 
    group_by(feature) %>% 
    mutate(
      ## Intercepts
      # Z=0; W=1 (PFS12=0; e.g. W=1)
      intercept2a = b0 + !!rlang::sym(focus_vars[3]), 
      # Z=1; W=1 (PFS12=1; e.g. W=1)
      intercept2b = b0 + !!rlang::sym(focus_vars[2]) + !!rlang::sym(focus_vars[3]) + !!rlang::sym(paste0(focus_vars[2],":",focus_vars[3])), # 
      
      ## Slopes without intercepts
      # Z=0; W=1 (PFS12=0; e.g. W1=1)
      slope2a = !!rlang::sym(focus_vars[1]) + !!rlang::sym(paste0(focus_vars[1],":",focus_vars[3])), 
      # Z=1; W=1 (PFS12=1; e.g. W1=1)
      slope2b = !!rlang::sym(focus_vars[1]) + !!rlang::sym(paste0(focus_vars[1],":",focus_vars[2])) + !!rlang::sym(paste0(focus_vars[1],":",focus_vars[3])) + !!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[3])),
      
      # Z=1 vs Z=0 for e.g. W=1 (e.g. W1=1) 
      avgdiff2 = !!rlang::sym(focus_vars[2]) + 
        !!rlang::sym(paste0(focus_vars[2],":",focus_vars[3])) +
        !!rlang::sym(paste0(focus_vars[1],":",focus_vars[2])) +
        !!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[3])),
      
      avgdiff2_x_t0 = 
        !!rlang::sym(focus_vars[2]) + # b2
        !!rlang::sym(paste0(focus_vars[2],":",focus_vars[3])) + 
        !!rlang::sym(paste0(focus_vars[1],":",focus_vars[2])) * 0 +
        !!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[3])) * 0, # at t1 (visit=0)
      
      avgdiff2_x_t1 = 
        !!rlang::sym(focus_vars[2])  + # b2
        !!rlang::sym(paste0(focus_vars[2],":",focus_vars[3])) + 
        !!rlang::sym(paste0(focus_vars[1],":",focus_vars[2])) * 1 +
        !!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[3])) * 1, # at t1 (visit=1)
      
      avgdiff2_x_t2 = 
        !!rlang::sym(focus_vars[2])  + # b2
        !!rlang::sym(paste0(focus_vars[2],":",focus_vars[3])) + 
        !!rlang::sym(paste0(focus_vars[1],":",focus_vars[2])) * 2 +
        !!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[3])) * 2, # at t2 (visit=2)
      
      avgdiff2_x_t3 = 
        !!rlang::sym(focus_vars[2])  + 
        !!rlang::sym(paste0(focus_vars[2],":",focus_vars[3])) + 
        !!rlang::sym(paste0(focus_vars[1],":",focus_vars[2])) * 3 + 
        !!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[3])) * 3 # at t3 (visit=3)
      ) %>%  
        
    select(feature,
           slope2a, 
           slope2b, 
           avgdiff2_x_t0,
           avgdiff2_x_t1, 
           avgdiff2_x_t2,
           avgdiff2_x_t3)
  
  ranks_ls <- vector("list", 4)
  names(ranks_ls) <- c("avgdiff2_x_t0","avgdiff2_x_t1","avgdiff2_x_t2","avgdiff2_x_t3")
  
  for(i in names(ranks_ls)) {
    print(i)
    
    ranks_ls[[i]]$sig_increasing <-
      intercepts_and_slopes %>% 
      select(feature, !!rlang::sym(i)) %>% 
      group_by(feature) %>%
      ggdist::median_qi(!!rlang::sym(i), .width=c(0, 0.5, 0.75, 0.90, 0.95, 0.97)) %>% 
      mutate(feature=reorder(factor(feature), !!rlang::sym(i))) %>% 
      pivot_wider(feature, names_from=.width, values_from=.lower) %>%
      select(feature, p0=`0`, p50=`0.5`, p75=`0.75`, p90=`0.9`, p95=`0.95`, p97=`0.97`) %>%
      filter(!! rlang::sym(siglevel) > 0) %>%
      #filter(str_detect(feature, paste(literature_species, collapse="|"))) %>% 
      mutate(feature=factor(feature)) %>% 
      select(feature) %>% 
      pull() %>% 
      levels() 
    
    ranks_ls[[i]]$sig_decreasing <-
      intercepts_and_slopes %>% 
      select(feature, !!rlang::sym(i)) %>% 
      group_by(feature) %>%
      ggdist::median_qi(!!rlang::sym(i), .width=c(0, 0.5, 0.75, 0.90, 0.95, 0.97)) %>% 
      mutate(feature=reorder(factor(feature), !!rlang::sym(i))) %>% 
      pivot_wider(feature, names_from=.width, values_from=.upper) %>%
      select(feature, p0=`0`, p50=`0.5`, p75=`0.75`, p90=`0.9`, p95=`0.95`, p97=`0.97`) %>%
      filter(!! rlang::sym(siglevel) < 0) %>%
      #filter(str_detect(feature, paste(literature_species, collapse="|"))) %>% 
      mutate(feature=factor(feature)) %>% 
      select(feature) %>% 
      pull() %>% 
      levels()
  }
  
  ranks_ls2 <- vector("list", 2)
  
  names(ranks_ls2) <- c("slope2a","slope2b")
  
  included_taxa <- c(ranks_ls$avgdiff2_x_t0$sig_increasing, ranks_ls$avgdiff2_x_t0$sig_decreasing,
                     ranks_ls$avgdiff2_x_t1$sig_increasing, ranks_ls$avgdiff2_x_t1$sig_decreasing,
                     ranks_ls$avgdiff2_x_t2$sig_increasing, ranks_ls$avgdiff2_x_t2$sig_decreasing,
                     ranks_ls$avgdiff2_x_t3$sig_increasing, ranks_ls$avgdiff2_x_t3$sig_decreasing)
  
  for(i in names(ranks_ls2)) {
    print(i)
    
    ranks_ls2[[i]]$sig_increasing <-
      intercepts_and_slopes %>% 
      select(feature, !!rlang::sym(i)) %>% 
      group_by(feature) %>%
      ggdist::median_qi(!!rlang::sym(i), .width=c(0, 0.5, 0.75, 0.90, 0.95, 0.97)) %>% 
      mutate(feature=reorder(factor(feature), !!rlang::sym(i))) %>% 
      pivot_wider(feature, names_from=.width, values_from=.lower) %>%
      select(feature, p0=`0`, p50=`0.5`, p75=`0.75`, p90=`0.9`, p95=`0.95`, p97=`0.97`) %>%
      #filter(!! rlang::sym(siglevel) > 0) %>% # Either focus on features whoes posterior was differentially abundant (included_taxa) or those that have a posterior slope that does not cover 0 given the "siglevel"
      filter(feature %in% included_taxa) %>% 
      mutate(feature=factor(feature)) %>% 
      select(feature) %>% 
      pull() %>% 
      levels() 
    
    ranks_ls2[[i]]$sig_decreasing <-
      intercepts_and_slopes %>% 
      select(feature, !!rlang::sym(i)) %>% 
      group_by(feature) %>%
      ggdist::median_qi(!!rlang::sym(i), .width=c(0, 0.5, 0.75, 0.90, 0.95, 0.97)) %>% 
      mutate(feature=reorder(factor(feature), !!rlang::sym(i))) %>% 
      pivot_wider(feature, names_from=.width, values_from=.upper) %>%
      select(feature, p0=`0`, p50=`0.5`, p75=`0.75`, p90=`0.9`, p95=`0.95`, p97=`0.97`) %>%
      #filter(!! rlang::sym(siglevel) < 0) %>%
      filter(feature %in% included_taxa) %>% 
      mutate(feature=factor(feature)) %>% 
      select(feature) %>% 
      pull() %>% 
      levels()
  }
  
  ## Color labeling features shared between different visits
  
  # Z=0 vs Z=1 for W=1
  
  # Features shared between t1 and t0
  shared_t1_t0_a <- Reduce(intersect, list(c(ranks_ls$avgdiff2_x_t0$sig_increasing, ranks_ls$avgdiff2_x_t0$sig_decreasing),
                                           c(ranks_ls$avgdiff2_x_t1$sig_increasing, ranks_ls$avgdiff2_x_t1$sig_decreasing)))
  
  # Features shared between t2 and t1
  shared_t1_t2_a <- Reduce(intersect, list(c(ranks_ls$avgdiff2_x_t1$sig_increasing, ranks_ls$avgdiff2_x_t1$sig_decreasing),
                                           c(ranks_ls$avgdiff2_x_t2$sig_increasing, ranks_ls$avgdiff2_x_t2$sig_decreasing)))
  
  # Features shared between t3 and t2
  shared_t2_t3_a <- Reduce(intersect, list(c(ranks_ls$avgdiff2_x_t2$sig_increasing, ranks_ls$avgdiff2_x_t2$sig_decreasing),
                                           c(ranks_ls$avgdiff2_x_t3$sig_increasing, ranks_ls$avgdiff2_x_t3$sig_decreasing)))
  
  # Features shared between t3, t2, t1
  shared_t1_t2_t3_a <- Reduce(intersect, list(
    c(ranks_ls$avgdiff2_x_t3$sig_increasing, ranks_ls$avgdiff2_x_t3$sig_decreasing),
    c(ranks_ls$avgdiff2_x_t2$sig_increasing, ranks_ls$avgdiff2_x_t2$sig_decreasing),
    c(ranks_ls$avgdiff2_x_t1$sig_increasing, ranks_ls$avgdiff2_x_t1$sig_decreasing)))
  
  y_labels_t1_a <- setNames(rep("black",length(c(ranks_ls$avgdiff2_x_t1$sig_increasing, ranks_ls$avgdiff2_x_t1$sig_decreasing))),
                            rev(c(rev(ranks_ls$avgdiff2_x_t1$sig_increasing), rev(ranks_ls$avgdiff2_x_t1$sig_decreasing))))
  
  y_labels_t2_a <- setNames(rep("black",length(c(ranks_ls$avgdiff2_x_t2$sig_increasing, ranks_ls$avgdiff2_x_t2$sig_decreasing))),
                            rev(c(rev(ranks_ls$avgdiff2_x_t2$sig_increasing), rev(ranks_ls$avgdiff2_x_t2$sig_decreasing))))
  
  y_labels_t3_a <- setNames(rep("black",length(c(ranks_ls$avgdiff2_x_t3$sig_increasing, ranks_ls$avgdiff2_x_t3$sig_decreasing))),
                            rev(c(rev(ranks_ls$avgdiff2_x_t3$sig_increasing), rev(ranks_ls$avgdiff2_x_t3$sig_decreasing))))
  
  y_labels_t1_a[shared_t0_t1_a] <- "#74c476"
  y_labels_t2_a[shared_t1_t2_a] <- "#74c476"
  y_labels_t3_a[setdiff(shared_t2_t3_a, shared_t1_t2_t3_a)] <- "#74c476"
  y_labels_t3_a[shared_t1_t2_t3_a] <- "#005a32"
  
  plot_ls <- vector("list", 6)
  names(plot_ls) <- c("slope2a","slope2b","avgdiff2_x_t0","avgdiff2_x_t1","avgdiff2_x_t2","avgdiff2_x_t3")
  
  for(p in names(plot_ls)) {
    print(p)
    
    if (p=="avgdiff2_x_t1") {
      
      gp <- intercepts_and_slopes %>% 
        select(feature, !!rlang::sym(p)) %>% 
        group_by(feature) %>%
        ggdist::median_qi(!!rlang::sym(p), .width=c(0.5, 0.75, 0.90, 0.95, 0.97)) %>% 
        mutate(feature=reorder(factor(feature), !!rlang::sym(p))) %>% 
        filter(feature %in% c(ranks_ls[[p]]$sig_increasing, ranks_ls[[p]]$sig_decreasing)) %>%
        
        mutate(feature=str_remove_all(string = feature, pattern = "f__"),
               feature=str_remove_all(string = feature, pattern = "s__"),
               feature=str_remove_all(string = feature, pattern = "t__")) %>% 
        
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
          axis.text.y=ggtext::element_markdown(color=unname(y_labels_t2_a), size=7)) +
        labs(x="Log-Ratio Value", y=NULL, title=paste0("T1: ", str_split(focus_var,":")[[1]][2]," 1 vs 0 ", paste0(str_split(focus_var,":")[[1]][3],"=0"))) + 
        geom_vline(xintercept=0, linetype="dashed", color="darkgray")
      
      plot_ls[[p]] <- gp
      
    } else if (p=="avgdiff2_x_t2") {
      
      gp <- intercepts_and_slopes %>% 
        select(feature, !!rlang::sym(p)) %>% 
        group_by(feature) %>%
        ggdist::median_qi(!!rlang::sym(p), .width=c(0.5, 0.75, 0.90, 0.95, 0.97)) %>% 
        mutate(feature=reorder(factor(feature), !!rlang::sym(p))) %>% 
        filter(feature %in% c(ranks_ls[[p]]$sig_increasing, ranks_ls[[p]]$sig_decreasing)) %>%
        
        mutate(feature=str_remove_all(string = feature, pattern = "f__"),
               feature=str_remove_all(string = feature, pattern = "s__"),
               feature=str_remove_all(string = feature, pattern = "t__")) %>% 
        
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
          axis.text.y=ggtext::element_markdown(color=unname(y_labels_t3_a), size=7)) +
        labs(x="Log-Ratio Value", y=NULL, title=paste0("T2: ", str_split(focus_var,":")[[1]][2]," 1 vs 0 ", paste0(str_split(focus_var,":")[[1]][3],"=0"))) + 
        geom_vline(xintercept=0, linetype="dashed", color="darkgray")
      
      plot_ls[[p]] <- gp
      
    } else if (p=="avgdiff2_x_t3") { 
      
      gp <- intercepts_and_slopes %>% 
        select(feature, !!rlang::sym(p)) %>% 
        group_by(feature) %>%
        ggdist::median_qi(!!rlang::sym(p), .width=c(0.5, 0.75, 0.90, 0.95, 0.97)) %>% 
        mutate(feature=reorder(factor(feature), !!rlang::sym(p))) %>% 
        filter(feature %in% c(ranks_ls[[p]]$sig_increasing, ranks_ls[[p]]$sig_decreasing)) %>%
        
        mutate(feature=str_remove_all(string = feature, pattern = "f__"),
               feature=str_remove_all(string = feature, pattern = "s__"),
               feature=str_remove_all(string = feature, pattern = "t__")) %>% 
        
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
          axis.text.y=ggtext::element_markdown(color=unname(y_labels_t4_a), size=7)) +
        labs(x="Log-Ratio Value", y=NULL, title=paste0("T3: ", str_split(focus_var,":")[[1]][2]," 1 vs 0 ", paste0(str_split(focus_var,":")[[1]][3],"=0"))) + 
        geom_vline(xintercept=0, linetype="dashed", color="darkgray")
      
      plot_ls[[p]] <- gp
      
    } else if (p=="slope2a") {
      
      gp <- intercepts_and_slopes %>% 
        select(feature, !!rlang::sym(p)) %>% 
        group_by(feature) %>%
        ggdist::median_qi(!!rlang::sym(p), .width=c(0.5, 0.75, 0.90, 0.95, 0.97)) %>% 
        mutate(feature=reorder(factor(feature), !!rlang::sym(p))) %>% 
        filter(feature %in% c(ranks_ls2[[p]]$sig_increasing, ranks_ls2[[p]]$sig_decreasing)) %>%
        
        mutate(feature=str_remove_all(string = feature, pattern = "f__"),
               feature=str_remove_all(string = feature, pattern = "s__"),
               feature=str_remove_all(string = feature, pattern = "t__")) %>% 
        
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
          axis.text.y=ggtext::element_markdown(color=unname(y_labels_t2_a), size=7)) +
        labs(x="Log-Ratio Value", y=NULL, title=paste0("T1: ", str_split(focus_var,":")[[1]][2]," 1 vs 0 ", paste0(str_split(focus_var,":")[[1]][3],"=0"))) + 
        geom_vline(xintercept=0, linetype="dashed", color="darkgray")
      
      plot_ls[[p]] <- gp
      
    } else if (p=="slope2b") {
      
      gp <- intercepts_and_slopes %>% 
        select(feature, !!rlang::sym(p)) %>% 
        group_by(feature) %>%
        ggdist::median_qi(!!rlang::sym(p), .width=c(0.5, 0.75, 0.90, 0.95, 0.97)) %>% 
        mutate(feature=reorder(factor(feature), !!rlang::sym(p))) %>% 
        filter(feature %in% c(ranks_ls2[[p]]$sig_increasing, ranks_ls2[[p]]$sig_decreasing)) %>%
        
        mutate(feature=str_remove_all(string = feature, pattern = "f__"),
               feature=str_remove_all(string = feature, pattern = "s__"),
               feature=str_remove_all(string = feature, pattern = "t__")) %>% 
        
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
          axis.text.y=ggtext::element_markdown(color=unname(y_labels_t3_a), size=7)) +
        labs(x="Log-Ratio Value", y=NULL, title=paste0("T2: ", str_split(focus_var,":")[[1]][2]," 1 vs 0 ", paste0(str_split(focus_var,":")[[1]][3],"=0"))) + 
        geom_vline(xintercept=0, linetype="dashed", color="darkgray")
      
      plot_ls[[p]] <- gp
    }
  }
  return(plot_ls)
}

avg_RvsNR_combitherapy_p90 <- plot_RvsNR_combitherapy(pibble_fit=fit_species_clr, focus_var="visit:PFS12:combiIO", siglevel="p90")
