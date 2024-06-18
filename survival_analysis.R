# Code for survival analysis
ftbl <- read.csv("sgbs_baseline_survival.csv", check.names = F, row.names=1)
mdat <- read.csv("mdat_baseline_survival.csv", check.names = F)

ftbl <- as.matrix(zCompositions::cmultRepl(ftbl, label=0, method="CZM", z.delete=F)) # Bayesian 0 imputation

# Select SGBs for the "Longitudinal" balance
numerator_taxa <- c("f__Ruminococcaceae | s__Agathobaculum_butyriciproducens | t__SGB14993_group",
                    "f__Peptostreptococcaceae | s__Intestinibacter_bartlettii | t__SGB6140",
                    "f__Lachnospiraceae | s__Dorea_sp_AF24_7LB | t__SGB4571",
                    "f__Lactobacillaceae | s__Lactobacillus_gasseri | t__SGB7038_group",
                    "f__Lachnospiraceae | s__Lacrimispora_celerecrescens | t__SGB4868")

numerator_taxa_index <- match(numerator_taxa, colnames(ftbl))

denominator_taxa <- c("f__Ruminococcaceae | s__Ruthenibacterium_lactatiformans | t__SGB15271",
                      "f__Ruminococcaceae | s__Ruminococcaceae_unclassified_SGB15265 | t__SGB15265_group",
                      "f__Prevotellaceae | s__Prevotella_copri_clade_A | t__SGB1626",
                      "f__FGB602 | s__GGB1420_SGB1957 | t__SGB1957")

denominator_taxa_index <- match(denominator_taxa, colnames(ftbl))

balance_df <- mdat %>% mutate(balance_value=NA) %>% mutate(to_rowNames=sampleid) %>% column_to_rownames("to_rowNames")

# Compute balance
for(sample_i in rownames(ftbl)){
  balance_df[sample_i,]$balance_value <- log(exp(mean(log(ftbl[sample_i,numerator_taxa_index])))) - log(exp(mean(log(ftbl[sample_i,denominator_taxa_index]))))  
}

# Categorize balance based on median 
balance_df$balance_cat <- ifelse(balance_df$balance_value>quantile(balance_df$balance_value)[3], "high", "low")

# Plot KM curves
p <- 
  survfit2(Surv(os_months, status) ~ balance_cat, data = balance_df) %>% 
  ggsurvfit(linewidth=2) +
  add_censor_mark(color="black", size=1, shape=3) +
  add_risktable(risktable_stats = "n.risk") +
  labs(y="Overall survival probability", x="Survival time (months)") +
  scale_y_continuous(expand = c(0.025, 0), limits = c(0, 1), label = scales::label_percent()) +
  scale_x_continuous(expand = c(0.035, 0), breaks=seq(from=0, to=93, by=10)) +
  scale_color_manual(values=c("#bfb31d","#0cc1ab")) + 
  theme(
    legend.position = "none",
    legend.text=element_text(size=15, color="black"),
    strip.background=element_blank(),
    strip.text=element_text(size=15, color="black"),
    panel.grid.major=element_blank(), 
    panel.grid.minor=element_blank(), 
    panel.background=element_rect(fill="white"),
    panel.border=element_rect(colour="black", fill=NA, size=1),
    plot.title=element_text(size=18, color="black"),
    axis.title=element_text(size=18, color="black"),
    axis.text.y=element_text(size=15,color="black"),
    axis.text.x=element_text(size=15,color="black")) 

# Fit multivariable Cox regression
fit <- coxph(Surv(os_months, status) ~ balance_cat + sex + age + bmi + ppi + antibiotics + previous_therapy, data = balance_df)


