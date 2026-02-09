# Code for Extended P-scores paper

# Libraries ---------------------------

l<-c("netmeta","readxl","mvtnorm", "dplyr", "tidyr", "ggplot2", "xtable")
lapply(l, require, character.only = TRUE)

#load the functions required
source("Pscores_function.R")
source("league_table.R")
source("league_table_var.R")
source("intersect2.R")
source("Prepare_Multi.R")
source("Prepare_Single.R")
source("Prepare_function.R")
source("pscrs.R")
source("pscore_graph.R")
source("augustine-funcs.R")

# Load example dataset
data <- read_excel("Schizophrenia dataset.xls")

# One outcome: Efficacy --------------------------------

data_eff <-pairwise(treat=list(TreatA, TreatB),mean=list(PA,PB),sd=list(sdA,sdB),n=list(Na,Nb),
                          data = data, studlab = data$`study ID`, sm = "SMD")
net_eff<-netmeta(TE,seTE,treat1,treat2,studlab,sm="SMD",data=data_eff,
                   tol.multiarm = 0.1,ref="PBO")

eff <- list(net_eff)

# P-scores vs. CIV

# pscore_graph_t(x=eff,beta_=seq(0,1,0.01),type="B") # I get an out of bounds error when I run this - might want to double check this code

pscores_1outcome <- pscore_civs(eff, CIVs = list(efficacy = seq(0,0.5,length.out = 50)),
                                type = 'H')

tab <- subset(pscores_1outcome$all, efficacy %in% c(0, 0.5)) %>% group_by(efficacy) %>%
  mutate(mp = mean(Pscore), resid = Pscore-mp) %>% arrange(desc(Pscore)) %>% ungroup() %>% select(-c(poth, mp))%>%
  pivot_wider(values_from = c(Pscore, ranking, resid), names_from = efficacy, names_vary = "slowest")
print(xtable(tab, digits = 3), include.rownames = F)

# Plot P-scores vs CIV

plot_pscores(pscores_1outcome, room = c(0,0.15), title = "One Outcome: Efficacy")

# Plot of POTH vs CIV
plot_pothciv(pscores_1outcome, title = "One Outcome: Efficacy")

# POTH vs CIV over a wider range of CIVs
pscores_wide <- pscore_civs(eff, CIVs = list(efficacy = seq(-1,1,length.out = 50)),
                            type = 'H')
plot_pothciv(pscores_wide, title = "One Outcome: Efficacy")

# Summarise performance using AUPC
summ_eff <- aupc(pscores_1outcome)
plot_civsummary(summ_eff, title = "One Outcome: Efficacy")

# One outcome: Weight ---------------------------------------------

data_weight<-pairwise (treat=list(TreatA, TreatB),mean=list(meanweightA,meanweightB),sd=list(SDweightA,SDweightB),n=list(nweightA,nweightB),
                       data = data, studlab = data$`study ID`, sm = "SMD")
net_weight<-netmeta(TE,seTE,treat1,treat2,studlab,sm="SMD",data=data_weight,tol.multiarm = 0.5,tol.multiarm.se=0.5,ref="PBO")

wt <- list(net_weight)

pscores_wt <- pscore_civs(wt, CIVs = list(weight = seq(-0.5, 0, length.out = 50)),
                                type = 'H')

aupc(pscores_wt)

tab <- subset(pscores_wt$all, weight %in% c(0, 0.5)) %>% group_by(weight) %>%
  mutate(mp = mean(Pscore), resid = Pscore-mp) %>%
  arrange(desc(Pscore)) %>% ungroup() %>% select(-c(poth, mp))%>%
  pivot_wider(values_from = c(Pscore, ranking, resid), names_from = weight, names_vary = "slowest")
print(xtable(tab, digits = 3), include.rownames = F)
tab
# Two outcomes: Efficacy and weight gain -------------------------------------



correlation<-matrix(c(1,-0.5,-0.5,1),2,2)

pscores_2outcomes <- pscore_civs(list(net_eff, net_weight),
                                 CIVs = list(efficacy = seq(0, 0.5,length.out = 50),
                                             weight = seq(0, -0.5, length.out = 50)),
                                 correlation = correlation,
                                 type = c("H", "H"))

# Heat plot
pscores_heatplot(pscores_2outcomes, title = "Two Outcomes: Efficacy and Weight Gain")

# POTH bubble plot
ok <- pscores_pothplot(pscores_2outcomes, newgridsize = 11, title = "Two Outcomes: Efficacy and Weight Gain")

subset(ok$grid, POTH == max(POTH))

# Summary over 2 ranges of CIVs
summ_2outcomes <- aupc(pscores_2outcomes)

# Try comparing two hierarchies - compare a single CIV versus a small range
o1 <- pscore_civs(eff, CIVs = list(efficacy = seq(0.2, 0.4,length.out = 50)),
                  type = 'H')

h1 <- aupc(o1)
o2 <- pscore_civs(eff, CIVs = list(efficacy = 0.3),
                  type = 'H')

h2 <- o2$all

pscores_compare(x1 = h1, x2 = h2, name1 = "AUPC Rank,\nCIV (efficacy) 0.20-0.40",
                name2 = "Extended P-score Rank,\nCIV (efficacy) = 0.30")


# See the new Hierarchy
data.frame(Treatment = summ_2outcomes$Treatment,
           AUPC = round(summ_2outcomes$AUPC, digits = 5)) %>% arrange(-AUPC)

plot_civsummary(summ_2outcomes, title = "Two Outcomes: Efficacy and Weight Gain")
# can keep order the same to make comparison easier between different CIV ranges/number of outcomes
plot_civsummary(summ_2outcomes, ordering = "asis", title = "Two Outcomes: Efficacy and Weight Gain")

# Three outcomes: Efficacy, acceptability, weight gain ---------------------------
#Drop out all cause (Acceptability)
pair_drp<-pairwise (treat=list(TreatA, TreatB),event=list(DOtotalAP1,DOtotalAP2),n=list(NAP1r,NAP2r),
                    data = data, studlab = data$`study ID`, sm = "OR")
#transform OR to SMD
pair_drp_sel<-data.frame((pair_drp$TE*sqrt(3))/pi,((pair_drp$TE-1.96*pair_drp$seTE)*sqrt(3))/pi,((pair_drp$TE+1.96*pair_drp$seTE)*sqrt(3))/pi,pair_drp$studlab,pair_drp$treat1,pair_drp$treat2)
names(pair_drp_sel)<-c("TE","CIL","CIU","studlab","treat1","treat2")
pair_drp_sel$seTE<-(pair_drp_sel$CIU-pair_drp_sel$CIL)/3.92
net_drp<-netmeta(TE,seTE,treat1,treat2,studlab,sm="SMD",data=pair_drp_sel,tol.multiarm = 0.5,tol.multiarm.se=0.5,ref="PBO")

cn<-matrix(c(1,-0.5,0.2,-0.5,1,0,0.2,0,1),3,3)

pscores_3outcomes <- pscore_civs(x=list(net_eff, net_weight, net_drp),
                                 CIVs = list(efficacy = 0,
                                             weight = 0,
                                             dropout = 0),
                                 correlation =cn,
                                 type=c("H","H","H"))

# Try a range of CIVs for some outcomes

pscores_3outcomes2 <- pscore_civs(x=list(net_eff, net_weight, net_drp),
                                 CIVs = list(efficacy = seq(0,-0.3, length.out = 10),
                                             weight = 0,
                                             dropout = 0),
                                 correlation =cn,
                                 type=c("H","H","H"))

# Try again...

p_scores(x = list(net_eff, net_weight, net_drp), c(-0.267,0,0), correlation =cn,type=c("H","H","H")) # OLA is NaN?

subset(pscores_3outcomes2$all, is.nan(Pscore)) # get NaN P-scores with some treatments...
summ_3outcomes <- aupc(pscores_3outcomes2)
plot_civsummary(summ_3outcomes, title = "Three Outcomes: Efficacy (CIV [0, 0.5]), Weight (CIV = 0) and Dropout (CIV = 0)")

# Look at POTH values
pscores_3outcomes2$all %>% summarise(POTH = unique(poth))

# POTH versus CIV for the Only varying outcome
plot_pothciv(pscores_3outcomes2, title = "Three Outcomes: Efficacy, Weight (CIV = 0) and Dropout (CIV = 0)")


