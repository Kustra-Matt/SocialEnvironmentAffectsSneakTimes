#Code for:
#Social environment influences the temporal dynamics of sneak-spawning in a fish with alternative reproductive tactics.

#for questions contact Matt Kustra at mkustra@ucsc.edu

###Most analyses are stored as rds files and read in for faster visualizations. You can uncomment the model running code to run the actual models but they may take a while to run!
# 1. Load up required libraries -------------------------------------------
library(tidyverse)
library(factoextra)
library(FactoMineR)
library(patchwork)
library(data.table)
library(brms)
library(bayesplot)
library(tidybayes)
library(grid)
library(scales)
library(viridis)
library(bayestestR)
#writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

# 1. Theme setting for plots ---------------------------------------------
mytheme <- theme_classic()+ theme(
  legend.position = "bottom", # this puts legend on the bottom
 # this makes the axis titles in bold,
  axis.line = element_line(color = "black", size = 1.5), # Makes the axis line black and  thicker
  text = element_text(size = 12,family="sans"),
  axis.text = element_text(color="black")
) # makes all the text larger and bold

theme_set(mytheme)
# 2. Load up data ---------------------------------------------------------
Data_all <- read.csv("Data/Data.csv")
# glimpse at data
glimpse(Data_all)
#sample size by year for SI Table 5
Data_all%>%
  group_by(Year)%>%
  summarise(length(unique(File.Name)))

#sample size by year and male type for SI Table 5
Data_all%>%
  group_by(Year,Male)%>%
  count()

# 3. Separating Nest Data for PCA ----------------------------------------------------
Data_PCA <- Data_all %>%
  group_by(File.Name) %>%
  select(TOTAL.FEMALE.VISITS, TOTAL.FEMALES.SPAWNING, TOTAL.SPAWNS, SN_Sneak_N, SAT_Sneak_N, Sneak_N, Nest.male.to.SAT_AG, SAT.to.SN_AG, Sat.to.NM_SUB, Nest.male.to.SN_AG, Average.Sn) %>%
  #all values are the same for a nest, so median is just a shortcut to get one observation per nest
  summarise_all(median)

# 4. PCA analysis------------------------------------------------------------------
# Perform PCA analysis on selected variables
all.pca <- prcomp(Data_PCA[, c("TOTAL.FEMALE.VISITS", "TOTAL.FEMALES.SPAWNING", "TOTAL.SPAWNS", "SN_Sneak_N", "SAT_Sneak_N", "Sneak_N", "Nest.male.to.SAT_AG", "SAT.to.SN_AG", "Sat.to.NM_SUB", "Nest.male.to.SN_AG", "Average.Sn")], scale = TRUE, center = TRUE)
# quick visualization
fviz_eig(all.pca)
# get eigenvalues and variance explained for each PC
eig.val <- get_eigenvalue(all.pca)
eig.val
# get how much variation of PC are explained by variables
all.var <- get_pca_var(all.pca)
all.var$coord # Coordinates
all.var$contrib # Contributions to the PCs
all.var$cos2 # Quality of representation

all.ind <- get_pca_ind(all.pca)
# quick visualization of PCA
fviz_pca_var(all.pca,
  col.var = "contrib", # Color by contributions to the PC
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  repel = TRUE, # Avoid text overlapping
  axes = c(1, 2)
)

# Merge PCA's to main data frame
Data_PCA <- merge(Data_PCA, all.ind$coord, by = "row.names")
Data_all_PCA <- left_join(Data_all, Data_PCA)
# FLIP PCA SIGN FOR INTEPRETABILITY
Data_all_PCA$Dim.1 <- Data_all_PCA$Dim.1 * -1
# 5. Resampling PCA ------------------------------------------------------
## Now we need to test if PCA makes sense using methods from the following paper:
# https://onlinelibrary.wiley.com/doi/10.1111/evo.13835#evo13835-bib-0010
# psi test statistic, higher values indicate there is structure to the data so PCA is informative
psi <- function(vectoreig) {
  sum((vectoreig - 1)^2)
}
# phi test statistic, higher values indicate there is structure to the data so PCA is informative
phi <- function(vectoreig) {
  p <- length(vectoreig)
  ((sum(vectoreig^2) - p) / (p * (p - 1)))^0.5
}
## resampling function to test for significance in structure
# Data is the dataframe
# vars is a vector of variable names that will be in PCA
# reps is how many Reps you are running
# returns a list of (1) randomized overall pca scores
# (2) a list of the actual simulated draws
resamp_PCA <- function(Data, vars, reps = 10000) {
  # Allocating memory
  # psi
  ResultsPP <- data.frame(Psi = numeric(reps), Phi = numeric(reps), Sim = seq(1, reps))
  # EV results (eigen value)
  ResultsEV <- data.frame(Dim = numeric(reps * length(vars)), eigenvalue = numeric(reps * length(vars)), variance.percent = numeric(reps * length(vars)), cumulative.variance.percent = numeric(reps * length(vars)), Sim = 0)
  # correlations with variables
  ResultsC <- data.frame(matrix(ncol = (2 + length(vars)), nrow = (reps * length(vars))))
  names(ResultsC) <- c("Var", paste("Dim", seq(1:length(vars)), sep = "."), "Sim")
  #for loop to run the simulation
  for (i in 1:reps) {
    # Temp data frame by randomizing
    temp <- data.frame(lapply(Data[, vars], sample))
    # PCA on temp using the variables
    temp.pca <- prcomp(temp[, vars], scale = TRUE, center = TRUE)
    # calculate eigen values from variables
    eig.temp <- get_eigenvalue(temp.pca)
    # calculate Index loading vector
    temp.var <- data.frame((temp.pca$rotation^2) %*% diag(eig.temp$eigenvalue^2))
    # making rownames into an actual column
    setDT(eig.temp, keep.rownames = "Dim")
    # adding the sim reps
    eig.temp$Sim <- i
    # adding results of eigenvalues
    ResultsEV[((i - 1) * length(vars) + 1):(i * length(vars)), ] <- eig.temp
    # calculating overall PCA scores
    psiT <- psi(eig.temp$eigenvalue)
    phiT <- phi(eig.temp$eigenvalue)
    ResultsPP[i, ] <- list(psiT, phiT, i)
    # Store/name results into existing data frames
    names(temp.var) <- paste("F", names(data.frame(temp.pca$rotation)), "F", sep = "_")
    setDT(temp.var, keep.rownames = "Var")
    temp.var$Sim <- i
    ResultsC[((i - 1) * length(vars) + 1):(i * length(vars)), ] <- temp.var
  }
  # Summary things of actual data
  real.pca <- prcomp(Data[, vars], scale = TRUE, center = TRUE)
  eig.real <- get_eigenvalue(real.pca)
  # loadings
  real.var <- data.frame((real.pca$rotation^2) %*% diag(eig.real$eigenvalue^2))
  # PP
  psiR <- psi(eig.real$eigenvalue)
  phiR <- phi(eig.real$eigenvalue)
  SumPP <- ResultsPP %>%
    summarize(F_Psi_M = median(Psi), F_Psi_L = quantile(Psi, 0.025), F_Psi_U = quantile(Psi, 0.975), F_Phi_M = median(Phi), F_Phi_L = quantile(Phi, 0.025), F_Phi_U = quantile(Phi, 0.975)) %>%
    cbind(data.frame(real_Psi_V = psiR, real_Phi_V = phiR))
  # EV
  names(eig.real) <- paste("R", names(eig.real), "R", sep = "_")
  setDT(eig.real, keep.rownames = "Dim")
  SumEV <- ResultsEV %>%
    group_by(Dim) %>%
    summarize(F_eigenvalue_M = median(eigenvalue), F_eigenvalue_L = quantile(eigenvalue, 0.025), F_eigenvalue_U = quantile(eigenvalue, 0.975), F_variance.percent_M = median(variance.percent), F_variance.percent_L = quantile(variance.percent, 0.025), F_variance.percent_U = quantile(variance.percent, 0.975)) %>%
    left_join(eig.real)
  # Contributions
  # real.var
  names(real.var) <- paste(paste("Dim", seq(1:length(vars)), sep = "."), "R", "R", sep = "_")
  setDT(real.var, keep.rownames = "Var")
  SumC <- ResultsC %>%
    group_by(Var) %>%
    select(-Sim) %>%
    summarise_at(vars(-group_cols()), list(Fake_Median = ~ quantile(., probs = 0.5), Fake_L = ~ quantile(., probs = 0.025), Fake_U = ~ quantile(., probs = 0.975))) %>%
    left_join(real.var)

  return(list("Sums" = list("PP" = SumPP, "EV" = SumEV, "Con" = SumC), "Sims" = list("PPSims" = ResultsPP, "EVSims" = ResultsEV, "ConSims" = ResultsC)))
}

# 6. Resampling Plots ----------------------------------------------------
# Make vector of variables for PCA
vars <- c("TOTAL.FEMALE.VISITS", "TOTAL.FEMALES.SPAWNING", "TOTAL.SPAWNS", "SN_Sneak_N", "SAT_Sneak_N", "Sneak_N", "Nest.male.to.SAT_AG", "SAT.to.SN_AG", "Sat.to.NM_SUB", "Nest.male.to.SN_AG", "Average.Sn")
# subset data to just those columns
Data <- Data_PCA[, vars]
# set seed for analysis and replication
set.seed(10)
# Uncomment the line below to actually run the resampling function which takes a few minutes
#Resamps <- resamp_PCA(Data, vars)
#alternatively just read in the simulation results for faster checking
Resamps <-readRDS("ModelResults/PCA_Randomization_Results.rds")
# extract the summary from psi/phi (aka is the data structured?)
PP <- Resamps$Sums$PP
# extract the summary from eigenvalue and variance explained
EV <- Resamps$Sums$EV
# extract the summary from how variables contribute to PC's
Con <- Resamps$Sums$Con
# 7. Psi and Phi (statistics High values indicate the data is structured) --------------------------------------------------------------
# Making summary into long format for ploting
PPP <- PP %>%
  pivot_longer(cols = everything(), names_to = c("Cat", "Statistic", "Measurement"), names_pattern = "^(.*)_(.*)_(.*)", values_to = "Value") %>%
  pivot_wider(names_from = Measurement, values_from = Value) %>%
  mutate(Vals = ifelse(Cat == "F", M, V))
# Fake median psi, fake median phi, real psi, real phi
PPP$Vals
# Upper 97.5% quantile of psi, Upper 97.5% quantile of phi
PPP$U

# quick visualization of Psi and Phi (not shown in figure)
ggplot(data = PPP, aes(x = Cat, y = Vals, color = Cat)) +
  geom_point() +
  facet_wrap(~Statistic, scales = "free") +
  geom_errorbar(aes(ymin = L, ymax = U)) +
  ylab("Value")

# 8. EV helps us select PCs to use SI Figure 2 --------------------------------------------------------------
# Formating results into long format for plotting
EVP <- EV %>%
  pivot_longer(cols = c(-Dim), names_to = c("Cat", "Statistic", "Measurement"), names_pattern = "^(.*)_(.*)_(.*)", values_to = "Value") %>%
  pivot_wider(names_from = Measurement, values_from = Value) %>%
  mutate(Vals = ifelse(Cat == "F", M, R), Dim = as.numeric(str_extract(Dim, "[^.]*$"))) %>%
  filter(Statistic != "cumulative.variance.percent")

# Making SI figure 2
(pcaplot1 <- ggplot(data = EVP[EVP$Statistic == "variance.percent", ], aes(x = Dim, y = Vals, color = Cat)) +
  geom_point() +
  geom_line() +
  ylab("Percent variance explained") +
  geom_errorbar(aes(ymin = L, ymax = U)) +
  scale_x_continuous(breaks = unique(EVP$Dim)) +
  scale_color_manual(values = c("#FF0015", "#ACD9C3"), labels = c("Randomized", "Actual"), name = "") +
  xlab("Principal Component") +
  theme(legend.position = c(0.8,0.87)))

#ggsave("SIFigure2.pdf",height = 160,width=183,units="mm")
# 9. Contributions of variables to each PC SI Figure 4 --------------------------------------------------------------
# Formating data frame so contributions can be plotted
ConP <- Con %>%
  pivot_longer(cols = c(-Var), names_to = c("Statistic", "Cat", "Measurement"), names_pattern = "^(.*)_(.*)_(.*)", values_to = "Value") %>%
  pivot_wider(names_from = Measurement, values_from = Value) %>%
  mutate(Vals = ifelse(Cat == "Fake", Median, R), Statistic = as.numeric(str_extract(Statistic, "[^.]*$")))

# setting up labels for plots
varLabels <- c("TOTAL.SPAWNS" = "# Spawns", "TOTAL.FEMALES.SPAWNING" = "# Spawning females", "TOTAL.FEMALE.VISITS" = "# Female visits", "Sneak_N" = "# Sneaks", "SN_Sneak_N" = "# Sneaker sneaks", "SAT_Sneak_N" = "# Satellite sneaks", "Average.Sn" = "Average # sneakers", "Nest.male.to.SAT_AG" = "# NM to SAT agression", "SAT.to.SN_AG" = "# SAT to SN agression", "Nest.male.to.SN_AG" = "# NM to SN agression", "Sat.to.NM_SUB" = "# SAT to NM submission")

# SI Figure 4.
ggplot(data = ConP[ConP$Statistic %in% c(3), ], aes(x = factor(Var, levels = rev(names(varLabels))), y = Vals, color = Cat)) +
  geom_point() +
  geom_errorbar(aes(ymin = L, ymax = U), show.legend = F) +
  ylab("Value") +
  facet_wrap(~Statistic, scales = "free_x", labeller = labeller(Statistic = c("1" = "PC 1", "2" = "PC 2", "3" = "PC 3"))) +
  coord_flip() +
  scale_x_discrete(labels = varLabels) +
  scale_color_manual(values = c("#FF0015", "#ACD9C3"), labels = c("Randomized", "Actual"), name = "") +
  xlab("") +
  ylab("Index loading")
# ggsave("SIFigure4.pdf",height = 160,width=183,units="mm")

# 10. Making SI Figure 3 (please note colors were changed manually after saving) -----------------------------------------------------
# PC contributions for biological interpretations
(SIFig3A <- ggplot(data = ConP[ConP$Statistic %in% c(1, 2), ], aes(x = factor(Var, levels = rev(names(varLabels))), y = Vals, color = Cat)) +
  geom_point() +
  geom_errorbar(aes(ymin = L, ymax = U), show.legend = F) +
  ylab("Value") +
  facet_wrap(~Statistic, scales = "free_x", labeller = labeller(Statistic = c("1" = "PC 1", "2" = "PC 2"))) +
  coord_flip() +
  scale_x_discrete(labels = varLabels) +
  scale_color_manual(values = c("#FF0015", "#ACD9C3"), labels = c("Randomized", "Actual"), name = "") +
  xlab("") +
  ylab("Index loading") +
  theme(legend.position = "top"))
# variable labels
varLabels <- c("TOTAL.SPAWNS" = "# Spawns", "TOTAL.FEMALES.SPAWNING" = "# Spawning females", "TOTAL.FEMALE.VISITS" = "# Female visits", "Sneak_N" = "# Sneaks", "SN_Sneak_N" = "# Sneaker sneaks", "SAT_Sneak_N" = "# Satellite sneaks", "Average.Sn" = "Average # sneakers", "Nest.male.to.SAT_AG" = "# NM to SAT agression", "SAT.to.SN_AG" = "# SAT to SN agression", "Nest.male.to.SN_AG" = "# NM to SN agression", "Sat.to.NM_SUB" = "# SAT to NM submission")

# changing row names for PC plot
row.names(all.pca$rotation) <- c("# FM visits", "# Spawning FM", "# Spawns", "# SN sneaks", "# SAT sneaks", "# Sneaks", "# NM to SAT AG", "# SAT to SN AG", "# SAT to NM SB", "# NM to SN AG", "Average # SN")
# FLIP SIGN
all.pca$rotation[, 1] <- all.pca$rotation[, 1] * -1
# PCA corelation circle
(SIFig3B <- fviz_pca_var(all.pca,
  labelsize = 3,
  repel = TRUE, # Avoid text overlapping
  col.var = "black",
  axes = c(1, 2), select.var = list(name = c("# Spawns", "# Spawning FM", "# FM visits", "# Sneaks", "# SN sneaks", "# SAT sneaks", "Average # SN", "# SAT to SN AG", "# SAT to NM SB"))
) + ggtitle("") + mytheme + xlab("Nest activity (PC 1)") + ylab("Male interactions (PC 2)"))

# Put it together
SIFig3A + SIFig3B + plot_annotation(tag_levels = c("A"), tag_suffix = ")")

# save file
# ggsave("SIFigure3.pdf",height = 160,width=183,units="mm")
#Please note that label and arrow color was edited in illustrator.

# 11.a Analysis 1 (Table 1) Sat vs SN sneaking delays  -------------------------------------------------------
#Filter out only sneaker spawns for many category
#this is because this analysis focuses on satellite male and sneaker differences.
Data_all2<-Data_all%>%
  filter(!(SN_Sneak_S>1&SAT_Sneak_S<1))
#get sample size by male type
Data_all2%>%
  group_by(Male)%>%
  count()
#Get sample sizes for each category
Data_all2%>%
  group_by(AloneS,Male)%>%
  count()

# defining "uninformative" prior
prior <- c(prior(normal(0, 1000), class = Intercept), prior(normal(0, 1000), class = b), prior(inv_gamma(0.001, 0.001), class = sd), prior(gamma(0.01, 0.01), class = shape))
# reordering male ID so everything is compared to Sneaker
Data_all2$Male <- factor(Data_all2$Male, levels = c("SN", "SAT"))
# To run the model uncomment the next line.

#Lagtime_modelAM <- brm(Lag.Time ~ Male * AloneS+(1 | Year/File.Name / SpawnNumber), data = Data_all2, family = Gamma(link = "log"), cores = 4, chains = 4, seed = 12, iter = 5000, control = list(adapt_delta = 0.999), prior = prior, warmup = 1000)

#Since these models can take a while you can just load up the results of the model
Lagtime_modelAM<-readRDS("ModelResults/Model_1.rds")
# look at model outputs
sum <- summary(Lagtime_modelAM)
print(sum, 4)

# Calculate evidence ratio for parameters
print(hypothesis(Lagtime_modelAM, "Intercept<0"), 4)
print(hypothesis(Lagtime_modelAM, "MaleSAT<0"), 4)
print(hypothesis(Lagtime_modelAM, "AloneSMany<0"), 4)
print(hypothesis(Lagtime_modelAM, "AloneSPaired>0"), 4)
print(hypothesis(Lagtime_modelAM, "MaleSAT:AloneSMany<0"), 4)
print(hypothesis(Lagtime_modelAM, "MaleSAT:AloneSPaired<0"), 4)

#calculate Rsquared
bayes_R2(Lagtime_modelAM)

# Model fit diagnostics
# Look at chains
# RHAT
rhat_vals <- rhat(Lagtime_modelAM)
mcmc_rhat(rhat_vals)

# Look at how predicted draws match original distribution
pp_check(Lagtime_modelAM, nsamples = 500)

# use shinystan to do further analyses on model fit
# launch_shinystan(Lagtime_modelAM)


# 11.b.1 Alone ------------------------------------------------------------
#calculating independent contrasts for figure 2.
# Difference between sat and sneaker when spawning alone.
HovA <- hypothesis(Lagtime_modelAM, "exp(MaleSAT+Intercept) - exp(Intercept )<0")
HovA$hypothesis
# make data frame
DatahovA <- data.frame(HovA$samples)
# plot
ggplot(DatahovA, aes(x = H1)) +
  geom_density(fill = "#796FBD", alpha = 0.5, size = 1) +
  xlab("Difference in lagtime between \n satelite and sneakers") +
  ylab("Density") +
  geom_vline(aes(xintercept = HovA$hypothesis$Estimate, color = "estimate"), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = HovA$hypothesis$CI.Lower, color = "95%CI"), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = HovA$hypothesis$CI.Upper, color = "95%CI"), linetype = "dashed", size = 1) +
  scale_color_manual(name = "", values = c(estimate = "black", "95%CI" = "red"))

# sneaker lag time alone
HovASN <- hypothesis(Lagtime_modelAM, "exp(Intercept)>0")
HovASN$hypothesis
DatahovASN <- data.frame(HovASN$samples)
ggplot(DatahovASN, aes(x = H1)) +
  geom_density(fill = "#796FBD", alpha = 0.5, size = 1) +
  xlab("Sneaker Lag Time") +
  ylab("Density") +
  geom_vline(aes(xintercept = HovASN$hypothesis$Estimate, color = "estimate"), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = HovASN$hypothesis$CI.Lower, color = "95%CI"), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = HovASN$hypothesis$CI.Upper, color = "95%CI"), linetype = "dashed", size = 1) +
  scale_color_manual(name = "", values = c(estimate = "black", "95%CI" = "red"))

# Satelitte lag time alone
HovAST <- hypothesis(Lagtime_modelAM, "exp(MaleSAT+Intercept)>0")
HovAST$hypothesis
DatahovAST <- data.frame(HovAST$samples)
ggplot(DatahovAST, aes(x = H1)) +
  geom_density(fill = "#796FBD", alpha = 0.5, size = 1) +
  xlab("Satellite Lag Time") +
  ylab("Density") +
  geom_vline(aes(xintercept = HovAST$hypothesis$Estimate, color = "estimate"), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = HovAST$hypothesis$CI.Lower, color = "95%CI"), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = HovAST$hypothesis$CI.Upper, color = "95%CI"), linetype = "dashed", size = 1) +
  scale_color_manual(name = "", values = c(estimate = "black", "95%CI" = "red"))
# 11.b.2 Many ------------------------------------------------------------
# Difference between sat and sneaker when spawning with many sneakers
HovM <- hypothesis(Lagtime_modelAM, "exp(MaleSAT+Intercept+ AloneSMany+MaleSAT:AloneSMany ) - exp(Intercept + AloneSMany )<0")
HovM$hypothesis
DatahovM <- data.frame(HovM$samples)
ggplot(DatahovM, aes(x = H1)) +
  geom_density(fill = "#796FBD", alpha = 0.5, size = 1) +
  xlab("Difference in lagtime between \n satelite and sneakers") +
  ylab("Density") +
  geom_vline(aes(xintercept = HovM$hypothesis$Estimate, color = "estimate"), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = HovM$hypothesis$CI.Lower, color = "95%CI"), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = HovM$hypothesis$CI.Upper, color = "95%CI"), linetype = "dashed", size = 1) +
  scale_color_manual(name = "", values = c(estimate = "black", "95%CI" = "red"))
# sneaker time
HovMSN <- hypothesis(Lagtime_modelAM, "exp(Intercept + AloneSMany )>0")
HovMSN$hypothesis
DatahovMSN <- data.frame(HovMSN$samples)
ggplot(DatahovMSN, aes(x = H1)) +
  geom_density(fill = "#796FBD", alpha = 0.5, size = 1) +
  xlab("Difference in lagtime between \n satelite and sneakers") +
  ylab("Density") +
  geom_vline(aes(xintercept = HovMSN$hypothesis$Estimate, color = "estimate"), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = HovMSN$hypothesis$CI.Lower, color = "95%CI"), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = HovMSN$hypothesis$CI.Upper, color = "95%CI"), linetype = "dashed", size = 1) +
  scale_color_manual(name = "", values = c(estimate = "black", "95%CI" = "red"))

# satellite time
HovMST <- hypothesis(Lagtime_modelAM, "exp(MaleSAT+Intercept+ AloneSMany+MaleSAT:AloneSMany )>0")
HovMST$hypothesis
DatahovMST <- data.frame(HovMST$samples)
ggplot(DatahovMST, aes(x = H1)) +
  geom_density(fill = "#796FBD", alpha = 0.5, size = 1) +
  xlab("Difference in lagtime between \n satelite and sneakers") +
  ylab("Density") +
  geom_vline(aes(xintercept = HovMST$hypothesis$Estimate, color = "estimate"), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = HovMST$hypothesis$CI.Lower, color = "95%CI"), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = HovMST$hypothesis$CI.Upper, color = "95%CI"), linetype = "dashed", size = 1) +
  scale_color_manual(name = "", values = c(estimate = "black", "95%CI" = "red"))
# 11.b.3 Paired ------------------------------------------------------------
# difference between paired sat and sn
HovP <- hypothesis(Lagtime_modelAM, "exp(MaleSAT+Intercept+ AloneSPaired+MaleSAT:AloneSPaired) - exp(Intercept + AloneSPaired )<0")
HovP$hypothesis
DatahovP <- data.frame(HovP$samples)
ggplot(DatahovP, aes(x = H1)) +
  geom_density(fill = "#796FBD", alpha = 0.5, size = 1) +
  xlab("Difference in lagtime between \n satelite and sneakers") +
  ylab("Density") +
  geom_vline(aes(xintercept = HovP$hypothesis$Estimate, color = "estimate"), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = HovP$hypothesis$CI.Lower, color = "95%CI"), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = HovP$hypothesis$CI.Upper, color = "95%CI"), linetype = "dashed", size = 1) +
  scale_color_manual(name = "", values = c(estimate = "black", "95%CI" = "red"))
# SN time when spawning paired with sat
HovPSN <- hypothesis(Lagtime_modelAM, " exp(Intercept + AloneSPaired )>0")
HovPSN$hypothesis
DatahovPSN <- data.frame(HovPSN$samples)
ggplot(DatahovPSN, aes(x = H1)) +
  geom_density(fill = "#796FBD", alpha = 0.5, size = 1) +
  xlab("Sneaker lag-time") +
  ylab("Density") +
  geom_vline(aes(xintercept = HovPSN$hypothesis$Estimate, color = "estimate"), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = HovPSN$hypothesis$CI.Lower, color = "95%CI"), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = HovPSN$hypothesis$CI.Upper, color = "95%CI"), linetype = "dashed", size = 1) +
  scale_color_manual(name = "", values = c(estimate = "black", "95%CI" = "red"))
# Sat time when spawning with one sneaker
HovPST <- hypothesis(Lagtime_modelAM, "exp(MaleSAT+Intercept+ AloneSPaired+MaleSAT:AloneSPaired) >0")
HovPST$hypothesis
DatahovPST <- data.frame(HovPST$samples)
ggplot(DatahovPST, aes(x = H1)) +
  geom_density(fill = "#796FBD", alpha = 0.5, size = 1) +
  xlab("Satellite lag time") +
  ylab("Density") +
  geom_vline(aes(xintercept = HovPST$hypothesis$Estimate, color = "estimate"), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = HovPST$hypothesis$CI.Lower, color = "95%CI"), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = HovPST$hypothesis$CI.Upper, color = "95%CI"), linetype = "dashed", size = 1) +
  scale_color_manual(name = "", values = c(estimate = "black", "95%CI" = "red"))


# comparing paried and alone
HovPSN$hypothesis$Estimate
HovASN$hypothesis$Estimate
HovPST$hypothesis$Estimate
HovAST$hypothesis$Estimate
HovMSN$hypothesis$Estimate
HovMST$hypothesis$Estimate
#comparing sneaking delays
#pairwise snekaer
#Paired - alone
print(hypothesis(Lagtime_modelAM, "(exp(Intercept + AloneSPaired ))-(exp(Intercept )) >0"),4)
#paired - many
print(hypothesis(Lagtime_modelAM, "(exp(Intercept + AloneSPaired ))-(exp(Intercept + AloneSMany )) >0"),4)
#alone -many 
print(hypothesis(Lagtime_modelAM, "(exp(Intercept))-(exp(Intercept + AloneSMany )) >0"),4)
#pairwise Sat
#Paired - alone
print(hypothesis(Lagtime_modelAM, "(exp(MaleSAT+Intercept+ AloneSPaired+MaleSAT:AloneSPaired))-(exp(MaleSAT+Intercept)) <0"),4)
#paired - many
print(hypothesis(Lagtime_modelAM, "(exp(MaleSAT+Intercept+ AloneSPaired+MaleSAT:AloneSPaired))-(exp(MaleSAT+Intercept+ AloneSMany+MaleSAT:AloneSMany)) >0"),4)
#alone -many 
print(hypothesis(Lagtime_modelAM, "(exp(MaleSAT+Intercept))-(exp(MaleSAT+Intercept+ AloneSMany+MaleSAT:AloneSMany)) >0"),4)
#pairwise Diff
#Paired - alone
print(hypothesis(Lagtime_modelAM, "(exp(MaleSAT+Intercept+ AloneSPaired+MaleSAT:AloneSPaired) - exp(Intercept + AloneSPaired ))-(exp(MaleSAT+Intercept) - exp(Intercept ))<0"),4)
#paired - many
print(hypothesis(Lagtime_modelAM, "(exp(MaleSAT+Intercept+ AloneSPaired+MaleSAT:AloneSPaired) - exp(Intercept + AloneSPaired ))-(exp(MaleSAT+Intercept+ AloneSMany+MaleSAT:AloneSMany ) - exp(Intercept + AloneSMany ))<0"),4)
# 11.c.1 Plotting All differences (Fig 2B) ------------------------------------------------------------
DatahovA$Type <- "Single sneaker \n or satelitte"
DatahovM$Type <- "Multiple sneakers \n and satelitte"
DatahovP$Type <- "Single sneaker \n and satelitte"
Data_allP <- rbind(DatahovA, DatahovM, DatahovP) %>%
  select(H1, Type) %>%
  mutate(Type = factor(Type, levels = c("Single sneaker \n or satelitte","Single sneaker \n and satelitte", "Multiple sneakers \n and satelitte")))
Data_sum <- data.frame(Type = c("Single sneaker \n or satelitte", "Multiple sneakers \n and satelitte","Single sneaker \n and satelitte"), Low = c(HovA$hypothesis$CI.Lower, HovM$hypothesis$CI.Lower, HovP$hypothesis$CI.Lower), Median = c(HovA$hypothesis$Estimate, HovM$hypothesis$Estimate, HovP$hypothesis$Estimate), High = c(HovA$hypothesis$CI.Upper, HovM$hypothesis$CI.Upper, HovP$hypothesis$CI.Upper)) %>% mutate(Type = factor(Type, levels = c("Single sneaker \n or satelitte", "Single sneaker \n and satelitte", "Multiple sneakers \n and satelitte")))
Data_sum
#change alpah to 0.5 and color =white
(Dif <- ggplot(Data_allP, aes(x = H1)) +
    geom_density(fill = "#796FBD", alpha = 0.5, size = 1,color="black") +
    xlab("Difference in sneaking delay (satellite - sneaker)") +
    ylab("Probability density") +
    facet_grid(Type ~ .) +
    geom_vline(data = Data_sum, aes(xintercept = Median, color = "Estimate"), linetype = "dashed", size = 1) +
    geom_vline(data = Data_sum, aes(xintercept = Low, color = "95%CI"), linetype = "dashed", size = 1) +
    geom_vline(data = Data_sum, aes(xintercept = High, color = "95%CI"), linetype = "dashed", size = 1) +
    scale_color_manual(name = "", values = c(Estimate = "black", "95%CI" = "red")) +
    geom_vline(aes(xintercept = 0), size = 1, linetype = 3, color = "grey") +
    theme(legend.position = "top") +
    scale_y_continuous(n.breaks = 4))

# 11.c.2 Plotting All sn and sat separate(Fig 2A) ------------------------------------------------------------
# getting the estimates in difference scenarios
HovPSN$hypothesis$Estimate
HovASN$hypothesis$Estimate
HovPST$hypothesis$Estimate
HovAST$hypothesis$Estimate
HovMSN$hypothesis$Estimate
HovMST$hypothesis$Estimate
# adding column to categorize the spawning situation
DatahovASN$Type <- "Single sneaker \n or satelitte"
DatahovMSN$Type <- "Multiple sneakers \n and satelitte"
DatahovPSN$Type <- "Single sneaker \n and satelitte"
DatahovAST$Type <- "Single sneaker \n or satelitte"
DatahovMST$Type <- "Multiple sneakers \n and satelitte"
DatahovPST$Type <- "Single sneaker \n and satelitte"
# adding column to categorize the male type
DatahovASN$Male <- "Sneaker"
DatahovMSN$Male <- "Sneaker"
DatahovPSN$Male <- "Sneaker"
DatahovAST$Male <- "Satellite"
DatahovMST$Male <- "Satellite"
DatahovPST$Male <- "Satellite"
# Merging all the datasets together
Data_allPSNST <- rbind(DatahovASN, DatahovMSN, DatahovPSN, DatahovAST, DatahovMST, DatahovPST) %>%
  select(H1, Type, Male) %>%
  mutate(Type = factor(Type, levels = c("Single sneaker \n or satelitte", "Single sneaker \n and satelitte", "Multiple sneakers \n and satelitte")))
# creating summary data frame for credible intervals
Data_sumSNST <- data.frame(Type = c("Single sneaker \n or satelitte", "Multiple sneakers \n and satelitte", "Single sneaker \n and satelitte", "Single sneaker \n or satelitte", "Multiple sneakers \n and satelitte", "Single sneaker \n and satelitte"), Male = c("Sneaker", "Sneaker", "Sneaker", "Satellite", "Satellite", "Satellite"), Low = c(HovASN$hypothesis$CI.Lower, HovMSN$hypothesis$CI.Lower, HovPSN$hypothesis$CI.Lower, HovAST$hypothesis$CI.Lower, HovMST$hypothesis$CI.Lower, HovPST$hypothesis$CI.Lower), Median = c(HovASN$hypothesis$Estimate, HovMSN$hypothesis$Estimate, HovPSN$hypothesis$Estimate, HovAST$hypothesis$Estimate, HovMST$hypothesis$Estimate, HovPST$hypothesis$Estimate), High = c(HovASN$hypothesis$CI.Upper, HovMSN$hypothesis$CI.Upper, HovPSN$hypothesis$CI.Upper, HovAST$hypothesis$CI.Upper, HovMST$hypothesis$CI.Upper, HovPST$hypothesis$CI.Upper)) %>% mutate(Type = factor(Type, levels = c("Single sneaker \n or satelitte", "Single sneaker \n and satelitte", "Multiple sneakers \n and satelitte")))

# make lagtime  plots (Fig 2A)
(lt <- ggplot(Data_allPSNST, aes(x = H1, fill = Male)) +
    geom_density(alpha = 0.5, size = 1,color="black") +
    xlab("Sneaking delay (seconds)") +
    ylab("Probability density") +
    facet_grid(Type ~ .) +
    geom_vline(data = Data_sumSNST, aes(xintercept = Median, color = "Estimate"), linetype = "dashed", size = 1) +
    geom_vline(data = Data_sumSNST, aes(xintercept = Low, color = Male), linetype = "dashed", size = 1) +
    geom_vline(data = Data_sumSNST, aes(xintercept = High, color = Male), linetype = "dashed", size = 1) +
    scale_color_manual(name = "", c("Estimate", "95%CI", "95%CI"), values = c(Estimate = "black", Sneaker = "#D4BE62", Satellite = "#13D1E2")) +
    theme(legend.position = "top") +
    scale_fill_manual(c("Satellite", "Sneaker"), values = c("#13D1E2", "#D4BE62"), name = "")+theme(text=element_text(family="Arial")))
lt
# Put Figure 2 together
lt / Dif + plot_annotation(tag_levels = c("A"),tag_suffix =")") & theme(legend.position = "top")

# Save Figure 2
#ggsave("Plots/Figure2_unedited.pdf",height = 160*1.7,width=140*1.5,units="mm")
#please note that fish images and order of legend were later edited in illustrator. 

# 11.d Analysis 1 Supplement (SI Table 4 sneaker spawns only) -------------------------------------------------------
#Supplemental analysis just looking at sneakers 

#first make graph of sample sizes for number of spawns at different numbers of males.

# Number of events 
(SIFigure1a<-Data_all%>%
  distinct(Year,File.Name,SpawnNumber,.keep_all = T)%>%
  group_by(SN_Sneak_S,SAT_Sneak_S,TotalSneaks_S)%>%
  summarise(Count=n())%>%
  mutate(`Sat. present?`=ifelse(SAT_Sneak_S==0,"No","Yes"))%>%
  filter(SN_Sneak_S >1|(SN_Sneak_S ==1&SAT_Sneak_S>0))%>%
  ggplot(aes(x=TotalSneaks_S,y=Count,fill=`Sat. present?`))+
  geom_bar(stat="identity",position = position_dodge())+
  xlab("Number of parasitic males")+
  ylab("Number of spawns")+
  scale_fill_manual(c("Sneaker only", "Sneakers and satellite"), values = c("#D4BE62", "#13D1E2"), name = "")+
  scale_x_continuous(breaks=c(2,3,4,5,6,7)))
(SIFigure1b<-Data_all%>%
    distinct(Year,File.Name,SpawnNumber,.keep_all = T)%>%
    group_by(SN_Sneak_S,SAT_Sneak_S,TotalSneaks_S)%>%
    summarise(Count=n())%>%
    filter(SN_Sneak_S >1|(SN_Sneak_S ==1&SAT_Sneak_S>0))%>%
    mutate(`Sat. present?`=ifelse(SAT_Sneak_S==0,"No","Yes"))%>%
    mutate(Category=ifelse(TotalSneaks_S==2,"2","3+"))%>%
    group_by(`Sat. present?`,Category)%>%
    summarize(Count=sum(Count))%>%
    ggplot(aes(x=Category,y=Count,fill=`Sat. present?`))+
    geom_bar(stat="identity",position = position_dodge())+
    xlab("Number of parasitic male category")+
    ylab("Number of spawns")+
    scale_fill_manual(c("Sneaker only", "Sneakers and satellite"), values = c("#D4BE62", "#13D1E2"), name = ""))
# Put together figure.
((SIFigure1a+theme(legend.position = "none"))|(SIFigure1b+theme(legend.position = c(0.8,1))))+
  plot_annotation(tag_levels = c("A"),tag_suffix =")")
#Save Figure
#ggsave("SIFigure1.pdf",height = 120*1.5,width=160*1.5,units="mm")

#Filter out satellite observations and redefine categories
Data_all$TotalSneaks_S
Data_allSNO <- Data_all %>%
  filter(Male != "SAT") %>%
  mutate(Cat = ifelse(
    SAT_Sneak_S == 0 & SN_Sneak_S == 1,"Sn only",
    ifelse(SAT_Sneak_S == 1 & SN_Sneak_S == 1,"Sn + Sat",
    ifelse(SAT_Sneak_S == 1 & SN_Sneak_S > 1,"Many Sn + Sat",
    ifelse(SAT_Sneak_S == 0 & SN_Sneak_S == 2, "Two Sn", "Many Sn")
      )))) %>%
  mutate(Cat = factor(Cat,
    levels = c("Sn only", "Sn + Sat", "Two Sn", "Many Sn + Sat", "Many Sn"))) %>%
  filter(Cat %in% c("Sn only", "Two Sn", "Many Sn"))


# defining "uninformative" prior
prior <- c(prior(normal(0, 1000), class = Intercept), prior(normal(0, 1000), class = b), prior(inv_gamma(0.001, 0.001), class = sd), prior(gamma(0.01, 0.01), class = shape))

# Uncomment the line below to run the model

#Lagtime_modelAMNum<- brm(Lag.Time ~ Cat+(1 | Year/File.Name / SpawnNumber), data = Data_allSNO, family = Gamma(link = "log"), cores = 4, chains = 4, seed = 12, iter = 5000, control = list(adapt_delta = 0.999), prior = prior, warmup = 1000)
#load up model for faster visualization/checking

Lagtime_modelAMNum<-readRDS("ModelResults/Model_SN_Only.rds")
# look at model outputs
sum <- summary(Lagtime_modelAMNum)
print(sum, 4)

# Calculate evidence ratio for parameters
print(hypothesis(Lagtime_modelAMNum, "Intercept<0"), 4)
print(hypothesis(Lagtime_modelAMNum, "CatTwoSn>0"), 4)
print(hypothesis(Lagtime_modelAMNum, "CatManySn<0"), 4)

#calculate Rsquared
bayes_R2(Lagtime_modelAMNum)

# Model fit diagnostics
# Look at chains
# RHAT
rhat_vals <- rhat(Lagtime_modelAMNum)
mcmc_rhat(rhat_vals)

# Look at how predicted draws match original distribution
pp_check(Lagtime_modelAMNum, nsamples = 500)

# use shinystan to do further analyses on model fit
# launch_shinystan(Lagtime_modelAMNum)

# 11.e.1 Alone ------------------------------------------------------------
#Getting independent contrasts for SI Figure 5 based on the sneaker only model.

# sneaker lag time alone
HovASN2 <- hypothesis(Lagtime_modelAMNum, "exp(Intercept)>0")
HovASN2$hypothesis
DatahovASN2 <- data.frame(HovASN2$samples)
ggplot(DatahovASN2, aes(x = H1)) +
  geom_density(fill = "#796FBD", alpha = 0.5, size = 1) +
  xlab("Sneaker Lag Time") +
  ylab("Density") +
  geom_vline(aes(xintercept = HovASN2$hypothesis$Estimate, color = "estimate"), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = HovASN2$hypothesis$CI.Lower, color = "95%CI"), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = HovASN2$hypothesis$CI.Upper, color = "95%CI"), linetype = "dashed", size = 1) +
  scale_color_manual(name = "", values = c(estimate = "black", "95%CI" = "red"))

# 11.e.2 Many ------------------------------------------------------------
# Difference between sat and sneaker when spawning with many sneakers
# sneaker time
HovMSN2 <- hypothesis(Lagtime_modelAMNum, "exp(Intercept + CatManySn )>0")
DatahovMSN2 <- data.frame(HovMSN2$samples)
ggplot(DatahovMSN2, aes(x = H1)) +
  geom_density(fill = "#796FBD", alpha = 0.5, size = 1) +
  xlab("Difference in lagtime between \n satelite and sneakers") +
  ylab("Density") +
  geom_vline(aes(xintercept = HovMSN2$hypothesis$Estimate, color = "estimate"), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = HovMSN2$hypothesis$CI.Lower, color = "95%CI"), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = HovMSN2$hypothesis$CI.Upper, color = "95%CI"), linetype = "dashed", size = 1) +
  scale_color_manual(name = "", values = c(estimate = "black", "95%CI" = "red"))
# 11.e.3 Paired ------------------------------------------------------------
# difference between paired sat and sn
# SN time when spawning with another sneaker
HovPSN2 <- hypothesis(Lagtime_modelAMNum, " exp(Intercept + CatTwoSn)>0")
HovPSN2$hypothesis
DatahovPSN2 <- data.frame(HovPSN2$samples)
ggplot(DatahovPSN2, aes(x = H1)) +
  geom_density(fill = "#796FBD", alpha = 0.5, size = 1) +
  xlab("Sneaker lag-time") +
  ylab("Density") +
  geom_vline(aes(xintercept = HovPSN2$hypothesis$Estimate, color = "estimate"), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = HovPSN2$hypothesis$CI.Lower, color = "95%CI"), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = HovPSN2$hypothesis$CI.Upper, color = "95%CI"), linetype = "dashed", size = 1) +
  scale_color_manual(name = "", values = c(estimate = "black", "95%CI" = "red"))
# 11.f Plotting All sn categories (SI Fig 5) ------------------------------------------------------------
# adding column to categorize the spawning situation
DatahovASN2$Type <- "Single sneaker"
DatahovMSN2$Type <- "More than two sneakers"
DatahovPSN2$Type <- "Two sneakers"
# adding column to categorize the male type
DatahovASN2$Male <- "Sneaker"
DatahovMSN2$Male <- "Sneaker"
DatahovPSN2$Male <- "Sneaker"
# Merging all the datasets together
Data_allPSNST2 <- rbind(DatahovASN2, DatahovMSN2, DatahovPSN2) %>%
  select(H1, Type, Male) %>%
  mutate(Type = factor(Type, levels = c("Single sneaker", "Two sneakers", "More than two sneakers")))
# creating summary data frame for credible intervals
Data_sumSNST2 <- data.frame(Type = c("Single sneaker", "More than two sneakers", "Two sneakers"), Male ="Sneaker", Low = c(HovASN2$hypothesis$CI.Lower, HovMSN2$hypothesis$CI.Lower, HovPSN2$hypothesis$CI.Lower), Median = c(HovASN2$hypothesis$Estimate, HovMSN2$hypothesis$Estimate, HovPSN2$hypothesis$Estimate), High = c(HovASN2$hypothesis$CI.Upper, HovMSN2$hypothesis$CI.Upper, HovPSN2$hypothesis$CI.Upper)) %>% mutate(Type = factor(Type, levels = c("Single sneaker", "Two sneakers", "More than two sneakers")))

# make lagtime  plots (Fig 2A)
(lt <- ggplot(Data_allPSNST2, aes(x = H1, fill = Male)) +
    geom_density(alpha = 0.5, size = 1,color="black",fill="#D4BE62") +
    xlab("Sneaking delay (seconds)") +
    ylab("Probability density") +
    facet_grid(Type ~ .) +
    geom_vline(data = Data_sumSNST2, aes(xintercept = Median, color = "Estimate"), linetype = "dashed", size = 1) +
    geom_vline(data = Data_sumSNST2, aes(xintercept = Low, color = Male), linetype = "dashed", size = 1) +
    geom_vline(data = Data_sumSNST2, aes(xintercept = High, color = Male), linetype = "dashed", size = 1) +
    scale_color_manual(name = "", c("Estimate", "95%CI", "95%CI"), values = c(Estimate = "black", Sneaker = "#D4BE62")) +
    theme(legend.position = "top") )

#ggsave("SIFigure5.pdf",height = 120*1.5,width=160*1.5,units="mm")

#pairwise sneaker comparisons
#Paired - alone
print(hypothesis(Lagtime_modelAMNum, "(exp(Intercept + CatTwoSn ))-(exp(Intercept )) >0"),4)
#paired - many
print(hypothesis(Lagtime_modelAMNum, "(exp(Intercept + CatTwoSn ))-(exp(Intercept + CatManySn )) >0"),4)
#alone -many 
print(hypothesis(Lagtime_modelAMNum, "(exp(Intercept))-(exp(Intercept + CatManySn )) >0"),4)

# 11.g Regularized prior sensitivity analysis (SI Table 1) ----------------
#Repeating analysis 1 (11A) but using a regularizing prior as a sensitivity analysis
# Set regularized prior
priorR <- c(prior(normal(0, 3), class = Intercept), prior(normal(0, 3), class = b), prior(inv_gamma(0.001, 0.001), class = sd), prior(gamma(0.01, 0.01), class = shape))

# Uncomment the line below to run the model. 
#Lagtime_modelAMD <- brm(Lag.Time ~ Male * AloneS + (1 | File.Name / SpawnNumber), data = Data_all2, family = Gamma(link = "log"), cores = 4, chains = 4, seed = 12, iter = 5000, control = list(adapt_delta = 0.999), warmup = 1000, prior = priorR)

#Load up the model instead of running since it can take a while to run.
Lagtime_modelAMD <-readRDS("ModelResults/Model_1_Regularized_Prior.rds")
# Look at model outps
sumD <- summary(Lagtime_modelAMD)
print(sumD, 4)
sumD$fixed
sumD$spec_pars
# Calculate evidence ratio for model parameters
print(hypothesis(Lagtime_modelAMD, "Intercept<0"), 4)
print(hypothesis(Lagtime_modelAMD, "MaleSAT<0"), 4)
print(hypothesis(Lagtime_modelAMD, "AloneSMany<0"), 4)
print(hypothesis(Lagtime_modelAMD, "AloneSPaired>0"), 4)
print(hypothesis(Lagtime_modelAMD, "MaleSAT:AloneSMany<0"), 4)
print(hypothesis(Lagtime_modelAMD, "MaleSAT:AloneSPaired<0"), 4)

# calculate bayesian R squared
bayes_R2(Lagtime_modelAMD)
# Model diagnostics
# rhat
rhat_valsD <- rhat(Lagtime_modelAMD)
mcmc_rhat(rhat_valsD)

#Check how it compares to distribution 
pp_check(Lagtime_modelAMD, nsamples = 500)
# 12.a Analysis 2: variation (Table 2; Figure 3) -----------------------------------------------
#second main analysis to look at how nest principal components influence the time-delays.
# get sample sizes
Data_all_PCA%>%
  group_by(Male)%>%
  count()
# setting priors
prior <- c(prior(normal(0, 1000), class = Intercept), prior(normal(0, 1000), class = b), prior(inv_gamma(0.001, 0.001), class = sd), prior(gamma(0.01, 0.01), class = shape))
# reorder factors so sneaker is baseline
Data_all_PCA$Male <- factor(Data_all_PCA$Male, levels = c("SN", "SAT"))

#dim 1 = nest activity, dim2=male interactions.
# Uncomment the line below to run the model.
#Lagtime_model2 <- brm(Lag.Time ~ (Male * Dim.1 ) +(Male *Dim.2 )+(1 | Year/File.Name / SpawnNumber), data = Data_all_PCA, family = Gamma(link = "log"), cores = 4, chains = 4, seed = 12, iter = 5000, warmup = 1000, control = list(adapt_delta = 0.99), prior = prior)

#read in the stored model results instead for faster visualization.
Lagtime_model2 <- readRDS("ModelResults/Model_2.rds")
# Summary of model
sum2 <- summary(Lagtime_model2)
print(sum2, 4)

  
# Calculate model fit Bayes Rsquared
bayes_R2(Lagtime_model2)

# Calculate evidence ratios of model parameters
print(hypothesis(Lagtime_model2, "Intercept<0"), 4)
print(hypothesis(Lagtime_model2, "MaleSAT<0"), 4)
print(hypothesis(Lagtime_model2, "Dim.1<0"), 4)
print(hypothesis(Lagtime_model2, "Dim.2<0"), 4)
print(hypothesis(Lagtime_model2, "MaleSAT:Dim.1>0"), 4)
print(hypothesis(Lagtime_model2, "MaleSAT:Dim.2<0"), 4)

# Look at model diagnostics
# RHAT
rhat_vals2 <- rhat(Lagtime_model2)
mcmc_rhat(rhat_vals2)

# How well predicted results match distribution of real data.
pp_check(Lagtime_model2, nsamples = 500)

# Create Figure 3
#marginal effect of nest activity + raw data
DF<-as.data.frame(conditional_effects(Lagtime_model2)$`Dim.1:Male`)
#change name for plotting
Data_all_PCA$effect2__<-Data_all_PCA$Male
#actual plot
(a<-ggplot(DF,aes(x=Dim.1,y=estimate__,group=effect2__))+
    geom_point(data=Data_all_PCA,mapping=aes(x=Dim.1,y=Lag.Time,color=effect2__,alpha=effect2__,shape=effect2__),size=2.5,alpha=0.8)+
  geom_ribbon(aes(ymin=lower__, ymax=upper__, fill = effect2__),alpha=0.5)+
  geom_line(size=1, aes(color=effect2__, linetype=effect2__,alpha=effect2__))+
  xlab("Nest activity (PC1)") + 
  ylab("Sneaking delay (seconds)") + 
    scale_shape_manual(name="",c("Sneaker", "Satellite"),values=c(16,17))+
    scale_linetype_manual(name="",c("Sneaker", "Satellite"),values=c(1,3))+
  scale_color_manual(name = "", c("Sneaker", "Satellite"), values = c("#D4BE62", "#13D1E2")) +
  scale_fill_manual(c("Sneaker", "Satellite"), values = c("#D4BE62", "#13D1E2"), name = "")+ 
  scale_alpha_manual(name="",c("Sneaker", "Satellite"), values = c(0.7,1.0))+
  theme(legend.position = "none",text=element_text(size=16)))
  


#marginal effect of male interactions + raw data
DF2<-as.data.frame(conditional_effects(Lagtime_model2)$`Dim.2:Male`)

(b<-ggplot(DF2,aes(x=Dim.2,y=estimate__,group=effect2__))+
    geom_point(data=Data_all_PCA,mapping=aes(x=Dim.2,y=Lag.Time,color=effect2__,alpha=effect2__,shape=effect2__),size=2.5,alpha=0.8)+
    geom_ribbon(aes(ymin=lower__, ymax=upper__, fill = effect2__),alpha=0.5)+
    geom_line(size=1, aes(color=effect2__, linetype=effect2__,alpha=effect2__))+
    xlab("Male interactions (PC2)") + 
    ylab("Sneaking delay (seconds)") + 
    scale_shape_manual(name="",c("Sneaker", "Satellite"),values=c(16,17))+
    scale_linetype_manual(name="",c("Sneaker", "Satellite"),values=c(1,3))+
    scale_color_manual(name = "", c("Sneaker", "Satellite"), values = c("#D4BE62", "#13D1E2")) +
    scale_fill_manual(c("Sneaker", "Satellite"), values = c("#D4BE62", "#13D1E2"), name = "")+ 
    scale_alpha_manual(name="",c("Sneaker", "Satellite"), values = c(0.7,1.0))+
    theme(legend.position = "none",text=element_text(size=16)))

#put together the plots using patchwork
(((a+theme(legend.position = c(0.7,0.9))) | b) &theme(text=element_text(size=16),axis.text = element_text(color="black"),aspect.ratio = 1)) + plot_annotation(tag_levels = c("A"), tag_suffix = ")")

#ggsave("Figure3.pdf", height = 120*1.1, width = 200*1.1, units = "mm")

# 12.b Regularized model 2 (SI Table 2)----------------------
#sensitivity analysis of using a regularizing prior on model reuslts 
#define the regualrizing priors
priorR <- c(prior(normal(0, 3), class = Intercept), prior(normal(0,  3), class = b), prior(inv_gamma(0.001, 0.001), class = sd), prior(gamma(0.01, 0.01), class = shape))

#need larger iterations for better effective sample size
#dim 1 = nest activity, dim2=male interactions.
#Uncomment the line below to run the model
#Lagtime_model2D <- brm(Lag.Time ~ Male * Dim.1 + (Male * Dim.2) + (1 | Year/File.Name / SpawnNumber), data = Data_all_PCA, family = Gamma(link = "log"), cores = 4, chains = 4, seed = 12, iter = 7000, warmup = 1000, control = list(adapt_delta = 0.99), prior = priorR)

#read in the stored model results instead for faster visualization.
Lagtime_model2D<-readRDS("ModelResults/Model_2_Regularized_Prior.rds")
#get summary of model
sum2D <- summary(Lagtime_model2D)
print(sum2D, 4)

#get rsquared value
bayes_R2(Lagtime_model2D)

#get evidence ratios for the different model effects
print(hypothesis(Lagtime_model2D, "Intercept<0"), 4)
print(hypothesis(Lagtime_model2D, "MaleSAT<0"), 4)
print(hypothesis(Lagtime_model2D, "Dim.1<0"), 4)
print(hypothesis(Lagtime_model2D, "Dim.2<0"), 4)
print(hypothesis(Lagtime_model2D, "MaleSAT:Dim.1>0"), 4)
print(hypothesis(Lagtime_model2D, "MaleSAT:Dim.2<0"), 4)

## Checking model fit 
# Look at chains
# RHAT
rhat_vals2D <- rhat(Lagtime_model2D)
mcmc_rhat(rhat_vals2D)

# Make sure predicted values is similar to real data
pp_check(Lagtime_model2D, nsamples = 500)


# 13. SI analysis all (SI Table 3)-----------------------------------------

#Third analysis with all model effects.
#reorder factors for analysis
Data_all_PCA$Male <- factor(Data_all_PCA$Male, levels = c("SN", "SAT"))
#dim 1 = nest activity, dim2=male interactions.
#need to filter out when there are no satillites in the many comparison
# Uncomment the line below to run the model
#Lagtime_modelAll <- brm(Lag.Time ~ Male * AloneS+Male*Dim.1+Male*Dim.2 + (1 | Year/File.Name / SpawnNumber), data = Data_all_PCA%>%filter(!(SN_Sneak_S>1&SAT_Sneak_S<1)), family = Gamma(link = "log"), cores = 4, chains = 4, seed = 12, iter = 5000, control = list(adapt_delta = 0.999), prior = prior, warmup = 1000)

#read in the stored model results instead for faster visualization.
Lagtime_modelAll<-readRDS("ModelResults/Model_3.rds")

# look at model outputs
sumAll <- summary(Lagtime_modelAll)
print(sumAll, 4)

#get rsquared value
bayes_R2(Lagtime_modelAll)

# Calculate evidence ratio for parameters
print(hypothesis(Lagtime_modelAll, "Intercept<0"), 4)
print(hypothesis(Lagtime_modelAll, "MaleSAT<0"), 4)
print(hypothesis(Lagtime_modelAll, "AloneSMany<0"), 4)
print(hypothesis(Lagtime_modelAll, "AloneSPaired>0"), 4)
print(hypothesis(Lagtime_modelAll, "MaleSAT:AloneSMany<0"), 4)
print(hypothesis(Lagtime_modelAll, "MaleSAT:AloneSPaired<0"), 4)
print(hypothesis(Lagtime_modelAll, "Dim.1<0"), 4)
print(hypothesis(Lagtime_modelAll, "Dim.2<0"), 4)
print(hypothesis(Lagtime_modelAll, "MaleSAT:Dim.1>0"), 4)
print(hypothesis(Lagtime_modelAll, "MaleSAT:Dim.2<0"), 4)

###Model diagnostics
# Look at chains
# RHAT
rhat_valsAll <- rhat(Lagtime_modelAll)
mcmc_rhat(rhat_valsAll)

# Make sure predicted values is similar to real data
pp_check(Lagtime_modelAll, nsamples = 500)

