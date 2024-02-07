# Function for computing and plotting the marginal average for colitis yes vs no, averaging over Z (PFS12), (therapy regimen) and W3 (PPI-use)

plot_colitisyes_vs_colitisno <- function(pibble_fit=pibble_fit, siglevel=siglevel) {
  
  focus_vars <- c("visit","PFS12yes","combiIOyes","colitisyes","ppiyes")
  avg_diffs <- as.data.frame.table(pibble_fit$Lambda) %>%
    pivot_wider(c(Var1,Var3), names_from=Var2, values_from=Freq) %>%
    select(feature=Var1, b0=`(Intercept)`, contains(focus_vars)) %>% 
    mutate(feature=str_remove_all(string = feature, pattern = "f__"),
           feature=str_remove_all(string = feature, pattern = "s__"),
           feature=str_remove_all(string = feature, pattern = "t__")) %>%
    
    group_by(feature) %>% 
    
    mutate(
      
      # Coltis=yes vs Colitis=no (so W2=1 vs W2=0) averaging across the levels of Z (PFS12), W1 (therapy regimen) and W3 (PPI-use)
      
      slopes_colitisno = 
        !!rlang::sym(focus_vars[1]) +
        0.5*(!!rlang::sym(paste0(focus_vars[1],":",focus_vars[2])) +
               !!rlang::sym(paste0(focus_vars[1],":",focus_vars[3])) +
               !!rlang::sym(paste0(focus_vars[1],":",focus_vars[5]))) +
        0.25*(!!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[3])) + 
                !!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[5]))),
      
      slopes_colitisyes = 
        !!rlang::sym(focus_vars[1]) +
        !!rlang::sym(paste0(focus_vars[1],":",focus_vars[4])) + 
        0.5*(!!rlang::sym(paste0(focus_vars[1],":",focus_vars[2])) +
               !!rlang::sym(paste0(focus_vars[1],":",focus_vars[3])) +
               !!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[4])) +
               !!rlang::sym(paste0(focus_vars[1],":",focus_vars[5]))) +
        0.25*(!!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[3])) + 
                !!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[5]))),
      
      avgdiff_coltisyesVScolitisno_x_t0 =
        !!rlang::sym(focus_vars[4]) +
        0.5*!!rlang::sym(paste0(focus_vars[2],":",focus_vars[4])) +
        
        0*(!!rlang::sym(paste0(focus_vars[1],":",focus_vars[4]))) +
        0.5*(
          0*(!!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[4])))
        ),
      
      avgdiff_coltisyesVScolitisno_x_t1 =
        !!rlang::sym(focus_vars[4]) +
        0.5*!!rlang::sym(paste0(focus_vars[2],":",focus_vars[4])) +
        
        1*(!!rlang::sym(paste0(focus_vars[1],":",focus_vars[4]))) +
        0.5*(
          1*(!!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[4])))
        ),
      
      avgdiff_coltisyesVScolitisno_x_t2 =
        !!rlang::sym(focus_vars[4]) +
        0.5*!!rlang::sym(paste0(focus_vars[2],":",focus_vars[4])) +
        
        2*(!!rlang::sym(paste0(focus_vars[1],":",focus_vars[4]))) +
        0.5*(
          2*(!!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[4])))
        ),
      
      avgdiff_coltisyesVScolitisno_x_t3 =
        !!rlang::sym(focus_vars[4]) +
        0.5*!!rlang::sym(paste0(focus_vars[2],":",focus_vars[4])) +
        
        3*(!!rlang::sym(paste0(focus_vars[1],":",focus_vars[4]))) +
        0.5*(
          3*(!!rlang::sym(paste0(focus_vars[1],":",focus_vars[2],":",focus_vars[4])))
        )
    ) %>%
    select(feature, 
           slopes_colitisno, slopes_colitisyes, 
           avgdiff_coltisyesVScolitisno_x_t0, avgdiff_coltisyesVScolitisno_x_t1, avgdiff_coltisyesVScolitisno_x_t2, avgdiff_coltisyesVScolitisno_x_t3)
  
  ranks_ls <- vector("list", 4)
  names(ranks_ls) <- c("avgdiff_coltisyesVScolitisno_x_t0","avgdiff_coltisyesVScolitisno_x_t1","avgdiff_coltisyesVScolitisno_x_t2","avgdiff_coltisyesVScolitisno_x_t3")
  
  for(i in names(ranks_ls)) {
    print(i)
    
    ranks_ls[[i]]$sig_increasing <-
      avg_diffs %>% 
      select(feature, !!rlang::sym(i)) %>% 
      group_by(feature) %>%
      ggdist::median_qi(!!rlang::sym(i), .width=c(0, 0.5, 0.75, 0.90, 0.95, 0.97)) %>% 
      mutate(feature=reorder(factor(feature), !!rlang::sym(i))) %>% 
      pivot_wider(feature, names_from=.width, values_from=.lower) %>%
      select(feature, p0=`0`, p50=`0.5`, p75=`0.75`, p90=`0.9`, p95=`0.95`, p97=`0.97`) %>%
      filter(!! rlang::sym(siglevel) > 0) %>%
      mutate(feature=factor(feature)) %>% 
      select(feature) %>% 
      pull() %>% 
      levels() 
    
    ranks_ls[[i]]$sig_decreasing <-
      avg_diffs %>% 
      select(feature, !!rlang::sym(i)) %>% 
      group_by(feature) %>%
      ggdist::median_qi(!!rlang::sym(i), .width=c(0, 0.5, 0.75, 0.90, 0.95, 0.97)) %>% 
      mutate(feature=reorder(factor(feature), !!rlang::sym(i))) %>% 
      pivot_wider(feature, names_from=.width, values_from=.upper) %>%
      select(feature, p0=`0`, p50=`0.5`, p75=`0.75`, p90=`0.9`, p95=`0.95`, p97=`0.97`) %>%
      filter(!! rlang::sym(siglevel) < 0) %>%
      mutate(feature=factor(feature)) %>% 
      select(feature) %>% 
      pull() %>% 
      levels()
  }
  
  included_taxa <- c(
    ranks_ls$avgdiff_coltisyesVScolitisno_x_t0$sig_increasing, ranks_ls$avgdiff_coltisyesVScolitisno_x_t0$sig_decreasing,
    ranks_ls$avgdiff_coltisyesVScolitisno_x_t1$sig_increasing, ranks_ls$avgdiff_coltisyesVScolitisno_x_t1$sig_decreasing,
    ranks_ls$avgdiff_coltisyesVScolitisno_x_t2$sig_increasing, ranks_ls$avgdiff_coltisyesVScolitisno_x_t2$sig_decreasing,
    ranks_ls$avgdiff_coltisyesVScolitisno_x_t3$sig_increasing, ranks_ls$avgdiff_coltisyesVScolitisno_x_t3$sig_decreasing)
  
  ranks_ls2 <- vector("list", 2)
  names(ranks_ls2) <- c("slopes_colitisyes","slopes_colitisno")
  
  for(i in names(ranks_ls2)) {
    print(i)
    
    ranks_ls2[[i]]$sig_increasing <-
      avg_diffs %>% 
      select(feature, !!rlang::sym(i)) %>% 
      group_by(feature) %>%
      ggdist::median_qi(!!rlang::sym(i), .width=c(0, 0.5, 0.75, 0.90, 0.95, 0.97)) %>% 
      mutate(feature=reorder(factor(feature), !!rlang::sym(i))) %>% 
      pivot_wider(feature, names_from=.width, values_from=.lower) %>%
      select(feature, p0=`0`, p50=`0.5`, p75=`0.75`, p90=`0.9`, p95=`0.95`, p97=`0.97`) %>%
      #filter(!! rlang::sym(siglevel) > 0) %>%
      filter(feature %in% included_taxa) %>% 
      mutate(feature=factor(feature)) %>% 
      select(feature) %>% 
      pull() %>% 
      levels() 
    
    ranks_ls2[[i]]$sig_decreasing <-
      avg_diffs %>% 
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
  
  plot_ls <- vector("list", 6)
  names(plot_ls) <- c("slopes_colitisyes","slopes_colitisno",
                      "avgdiff_coltisyesVScolitisno_x_t0","avgdiff_coltisyesVScolitisno_x_t1","avgdiff_coltisyesVScolitisno_x_t2","avgdiff_coltisyesVScolitisno_x_t3")
  
  for(p in names(plot_ls)) {
    print(p)
    
    if (p=="avgdiff_coltisyesVScolitisno_x_t0") {
      
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
        labs(x="Log-Ratio Value", y=NULL, title=paste0("T0: colitisyes vs colitisno")) + 
        geom_vline(xintercept=0, linetype="dashed", color="darkgray")
      
      plot_ls[[p]] <- gp
      
    } else if (p=="avgdiff_coltisyesVScolitisno_x_t1") {
      
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
        labs(x="Log-Ratio Value", y=NULL, title=paste0("T1: colitisyes vs colitisno")) + 
        geom_vline(xintercept=0, linetype="dashed", color="darkgray")
      
      plot_ls[[p]] <- gp
      
    } else if (p=="avgdiff_coltisyesVScolitisno_x_t2") {
      
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
        labs(x="Log-Ratio Value", y=NULL, title=paste0("T2: colitisyes vs colitisno")) + 
        geom_vline(xintercept=0, linetype="dashed", color="darkgray")
      
      plot_ls[[p]] <- gp
      
    } else if (p=="avgdiff_coltisyesVScolitisno_x_t3") { 
      
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
        labs(x="Log-Ratio Value", y=NULL, title=paste0("T3: colitisyes vs colitisno")) + 
        geom_vline(xintercept=0, linetype="dashed", color="darkgray")
      
      plot_ls[[p]] <- gp
      
    } else if (p=="slopes_colitisyes") { 
      
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
          axis.text.y=element_text(size=7, color="black")) +
        labs(x="Log-Ratio Value", y=NULL, title=paste0("Slope colitisyes")) + 
        geom_vline(xintercept=0, linetype="dashed", color="darkgray")
      
      plot_ls[[p]] <- gp
      
    } else if (p=="slopes_colitisno") { 
      
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
          axis.text.y=element_text(size=7, color="black")) +
        labs(x="Log-Ratio Value", y=NULL, title=paste0("Slope colitisno")) + 
        geom_vline(xintercept=0, linetype="dashed", color="darkgray")
      
      plot_ls[[p]] <- gp
    }
  }
  return(plot_ls)
}

colitisyes_vs_colitisno_p90 <- plot_colitisyes_vs_colitisno(pibble_fit=fit_species_clr, siglevel="p90")
