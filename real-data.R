library(arrow)
library(dplyr)
library(tidyr)
library(hyperinf)
library(ggrepel)
library(ggpubr)
library(parallel)

# see details
# https://huggingface.co/datasets/ayates/amr_portal/blob/main/README.md?utm_source=chatgpt.com
# important fields here: antibiotic_name, BioSample_ID, resistance_phenotype

system("wget https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/2025-11/phenotype.parquet")
df = read_parquet("phenotype.parquet")

# subset our bug of interest
bug.df = as.data.frame(df[df$species=="Klebsiella pneumoniae",])

# reframe resistance scores
df_bin <- bug.df %>%
  filter(resistance_phenotype %in% c("resistant", "susceptible")) %>%
  mutate(value = ifelse(resistance_phenotype == "resistant", 1, 0)) %>%
  group_by(BioSample_ID, species, antibiotic_name) %>%
  summarise(value = max(value), .groups = "drop")

# count how many drug entries we have
tmp = df_bin
counts = tmp %>% count(species, antibiotic_name) 
drugs = c()
ggplot(counts, aes(x=species, y=antibiotic_name, fill=n)) + geom_tile()

# count appearances for each drug and produce an ordered list
appears = data.frame()
for(drug in unique(df_bin$antibiotic_name)) {
  appears = rbind(appears, data.frame(drug=drug,
                                      species=length(which(counts$antibiotic_name==drug)),
                                      min = min(counts$n[counts$antibiotic_name==drug]),
                                      sum = sum(counts$n[counts$antibiotic_name==drug]))
  )
}
appears = appears[order(-appears$min),]

drugs = appears$drug[1:5]
wide_all_unc <- tmp[tmp$antibiotic_name %in% drugs,] %>%
  pivot_wider(names_from = antibiotic_name,
              values_from = value) 

wide_all_unc <- wide_all_unc %>%
  dplyr::select(all_of(c("BioSample_ID")), all_of(drugs))

all.df <- wide_all_unc %>%
  drop_na()

set.seed(1)

results = list()

for(expt in 1:6) {
  # get a random subset and look at its properties
  all.small = all.df[sample(1:nrow(all.df), 20),]
  # enforce start and end states so we can compare methods
  all.small[20,2:ncol(all.small)] = 1
  all.small[19,2:ncol(all.small)] = 0
  colSums(all.small[,2:ncol(all.small)])
  
  # run reversible and irreversible HyperMk and HyperHMM fits
  job1 = mcparallel({hyperinf(all.small, reversible=TRUE)})
  job2 = mcparallel({hyperinf(all.small, method="hypermk", reversible=FALSE)})
  job3 = mcparallel({hyperinf(all.small)})
  job4 = mcparallel({hyperinf(all.small, method="hypertraps")})
  results[[expt]] <- mccollect(list(job1, job2, job3, job4))
}

all.plots = ggarrange(
  plot_hyperinf_data(all.small),
  plot_hyperinf_comparative(results[[6]][c(1,3)], bend = 1, style="full", expt.names = c("R", "IR")) + theme(legend.position="none"), 
  plot_hyperinf_bubbles(results[[6]][c(1,3)], sqrt.trans = TRUE, p.scale = 0.2, expt.names = c("R", "IR")) ,
  nrow=1, widths=c(1,2,3), labels=c("A", "B", "C")
)

all.plots

sf = 3
png("rev-irrev-kp-compare.png", width=600*sf, height=240*sf, res=72*sf)
print(all.plots)
dev.off()

all.plots.full = ggarrange(
  plot_hyperinf_comparative(results[[1]][c(1,3)], bend = 1, style="full", expt.names = c("R", "IR")) + theme(legend.position="none"), 
  plot_hyperinf_bubbles(results[[1]][c(1,3)], sqrt.trans = TRUE, p.scale = 0.2, expt.names = c("R", "IR")) ,
  plot_hyperinf_comparative(results[[2]][c(1,3)], bend = 1, style="full", expt.names = c("R", "IR")) + theme(legend.position="none"), 
  plot_hyperinf_bubbles(results[[2]][c(1,3)], sqrt.trans = TRUE, p.scale = 0.2, expt.names = c("R", "IR")) ,
  plot_hyperinf_comparative(results[[3]][c(1,3)], bend = 1, style="full", expt.names = c("R", "IR")) + theme(legend.position="none"), 
  plot_hyperinf_bubbles(results[[3]][c(1,3)], sqrt.trans = TRUE, p.scale = 0.2, expt.names = c("R", "IR")) ,
  
  nrow=3, ncol=2, widths=c(2,3)
)

sf = 3
png("rev-irrev-kp-compare-full.png", width=500*sf, height=3*200*sf, res=72*sf)
print(all.plots.full)
dev.off()

results[[1]][[1]]$fitted_mk$AIC
results[[1]][[2]]$fitted_mk$AIC
results[[2]][[1]]$fitted_mk$AIC
results[[2]][[2]]$fitted_mk$AIC
results[[3]][[1]]$fitted_mk$AIC
results[[3]][[2]]$fitted_mk$AIC
