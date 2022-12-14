library(CNORode)
library(parallel)
library(ggrepel)
library(Metrics)
library(ComplexHeatmap)
library(tidyverse)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(circlize)

rm(list = ls())

#IMPORT CNOLIST AND PKNMODEL AND LOOKUP LIST
cnolist = CNOlist("99quant_per_Pep_cleaned_ERK.csv")

#Import PKN model
pknmodel = readSIF("sig.SIF")

#VIEW THE MODEL IN A NETWORK FASHION AND DO THE THREE STEP PRE-PROCESSING
model = preprocessing(cnolist, pknmodel, compression=FALSE)
plotModel(pknmodel, cnolist)

#FROM HERE CNORode WILL START BY SETTING DEFAULT PARAMETERS
ode_parameters=createLBodeContPars(model, LB_n = 1, LB_k = 0.001,
                                   LB_tau = 0.001, UB_n = 5, UB_k = 3.0, UB_tau = 10, default_n = 3,
                                   default_k = 0.5, default_tau = 1, opt_n = FALSE, opt_k = TRUE,
                                   opt_tau = TRUE, random = FALSE)


#VISUALIZE INITIAL SOLUTION
modelSim=plotLBodeModelSim(cnolist, model, ode_parameters,
                           timeSignals=c(0, 5, 15, 30, 60))

## Parameter Optimization
# essm
requireNamespace("MEIGOR")

paramsSSm=defaultParametersSSm()
paramsSSm$local_solver = "DHC"
paramsSSm$maxtime = Inf;#36000
paramsSSm$maxeval = 250000;
paramsSSm$atol=1e-6;
paramsSSm$reltol=1e-6;
paramsSSm$nan_fac=1000;
paramsSSm$dim_refset=30;
paramsSSm$n_diverse=1000;
paramsSSm$maxStepSize=Inf;
paramsSSm$maxNumSteps=10000;
transferFun=5;
paramsSSm$transfer_function = transferFun;

paramsSSm$lambda_tau=-1e-5
paramsSSm$lambda_k=5e-5
paramsSSm$bootstrap=F
paramsSSm$SSpenalty_fac=10#10
paramsSSm$SScontrolPenalty_fac=1000#1000

#RUN THE OPTIMIZATION PROCESS. TAKES A LONG TIME IF PARAMETERS ABOVE ARE NOT TUNED
opt_pars=parEstimationLBode(cnolist,model, method="essm", ode_parameters=ode_parameters, paramsSSm=paramsSSm)

#PLOT DATA WITH OPTIMIZED PARAMETERS
dir.create("final")
pdf("final/Optimization_results_final.pdf", width = 40, height = 20)
#pdf("ERK_feedback/Optimization_results_final_corrected_size.pdf", width = 200, height = 40)
simulatedData=plotLBodeFitness(cnolist, model, transfer_function=transferFun, ode_parameters=opt_pars, reltol = 1e-6, atol = 1e-6, maxStepSize = 1)
dev.off()

#simulatedData contains values of predicted curve

#Get optimization data
estim_pars = opt_pars$parValues[opt_pars$index_opt_pars]
estim_pars_names = opt_pars$parNames[opt_pars$index_opt_pars]

#Get df
dfCNORode <- data.frame(name = opt_pars$parNames, value = opt_pars$parValues)
dfCNORode <- subset(dfCNORode, value != 3.00000000)

#Clean df
dfCNORode <- dfCNORode %>% rowwise() %>% mutate(Kinase_1 = ifelse(str_detect(name, "_k_"), str_split(name, "_k_", simplify = TRUE)[1], str_split(name, "tau_", simplify = TRUE)[2]),
                                                Kinase_2 = ifelse(str_detect(name, "_k_"), str_split(name, "_k_", simplify = TRUE)[2], str_split(name, "tau_", simplify = TRUE)[2]),
                                                type = ifelse(str_detect(name, "_k_"), "k", "tau"))

#Also get TAU
returnTau <- function(x, y){
  temp <- data.frame(value = x, type = y)
  temp <- subset(temp, type == "tau")
  
  return(temp$value[1])
}

dfCNORode <- dfCNORode %>% group_by(Kinase_1) %>% mutate(Tau_1 = returnTau(x = value, y = type))
dfCNORode <- dfCNORode %>% group_by(Kinase_2) %>% mutate(Tau_2 = returnTau(x = value, y = type))

dfCNORodeK <- dfCNORode %>% filter(type == "k")
dfCNORodeK <- dfCNORodeK %>% filter(Tau_2 < 9)

#Add pathway information
RAS <- c("MP2K1-S222;MP2K2-S226",
         "MAPK1-Y187/T190",
         "CDK2-T160",
         "CDK7-S164",
         "MAPKAKP-T185",
         "CDK11B-T595",
         "CDK11A-S577;CDK11B-S589",
         "CDK1-T161")

PI3K <- c("PDPK1-S241",
          "PRKACA-T200;PRKACB-T202;PRKACG-T202",
          "PRKACA-T196;PRKACB-T198;PRKACG-T198",
          "PAK4-S474",
          "PKN2-T820",
          "PKN2-T816",
          "MARK1-T215;MARK2-T208;MARK3-T211",
          "MARK1-T219;MARK2-T212;MARK3-T215",
          "GSK3A-Y279")

PLCy <- c("PRKCD-T507",
          "PRKCG-T518",
          "PRKD1-S742;PRKD3-S731",
          "PRKD2-S710")

FGFR <- c("FGFR")

dfCNORodeK <- dfCNORodeK %>% rowwise() %>% mutate(pathway = ifelse(Kinase_2 %in% RAS, "RAS/RAF",
                                                                   ifelse(Kinase_2 %in% PI3K, "PI3K", 
                                                                          ifelse(Kinase_2 %in% PLCy, "PLCy", 
                                                                                 ifelse(Kinase_2 %in% FGFR, "FGFR", "RAS/RAF")))))
T
dfCNORodeK$name <- gsub("_k_", ">", dfCNORodeK$name)

####Scatter plot for tau & k comparison####
p <- ggplot(dfCNORodeK, aes(x = Tau_2, y = value, color = pathway)) + geom_point(show.legend = FALSE, size = 3) +
  stat_chull(aes(color = pathway, fill = pathway), alpha = 0.1, geom = "polygon", show.legend = FALSE)+ geom_text_repel(aes(label=name), show.legend = FALSE, color = "black", size = 3)
p + theme_bw() + theme(text = element_text(size = 20), panel.border = element_blank()) + facet_wrap(facets = vars(pathway)) +
  xlab("Node life-time (tau)") + ylab("Edge strength (k)")

ggsave("final/k_and_tau_transtest_final.pdf", width = 11, height = 6)

####Calculate RMSE####
#This creates a dataframe for the RMSE measurements
RMSEdf <- data.frame(FGF = character(),
                     Protein = character(),
                     Time = integer(),
                     simulatedValues = double())

timeVector <- c(0, 5, 15, 30, 60)
FGFs <- c("FGF2", "FGF3", "FGF4", "FGF10", "FGF19")
#loop per time
for (i in 1:length(timeVector)){
  #loop over the columns
  for (j in 1:length(colnames(simulatedData[[i]]))){
    print(paste("iteration ", i, " - colnumber ", j, sep =""))
    tempDF <- data.frame(FGF = FGFs,
                         Protein = rep(colnames(simulatedData[[i]])[j], 5),
                         Time = rep(timeVector[i], 5),
                         simulatedValues = simulatedData[[i]][,j])
    
    RMSEdf <- rbind(RMSEdf, tempDF)}}

#Add also the measurement values
#Signals are stored: cnolist@signals
signalVector <- c(cbind(as.vector(cnolist@signals$'0'),
                       as.vector(cnolist@signals$'5'),
                       as.vector(cnolist@signals$'15'),
                       as.vector(cnolist@signals$'30'),
                       as.vector(cnolist@signals$'60')))

RMSEdf$measuredValues <- signalVector

####calculate RMSE
#Sort by protein and FGF
#Measure RMSE by metrics

RMSECNORode <- function(t, meas, sim){
  measdf <- data.frame(tim = as.vector(t), measured = as.vector(meas), simulated = as.vector(sim))
  val <- rmse(measdf$meas, measdf$sim)
  return(val)
}

RMSEdf <- RMSEdf %>% group_by(Protein, FGF) %>% mutate(RMSE = RMSECNORode(t = Time, meas = measuredValues, sim = simulatedValues))

####Start with the dataframe of RMSE
RMSEdfMtrx <- RMSEdf %>% filter(Time == 0) %>% select(FGF, Protein, RMSE) %>% pivot_wider(values_from = RMSE, names_from = Protein)

mtrx <- data.matrix(RMSEdfMtrx[,2:ncol(RMSEdfMtrx)])
rownames(mtrx) <- RMSEdfMtrx$FGF

col <- colorRampPalette(brewer.pal(10, "RdYlBu"), bias = 1)(256)
col_fun = colorRamp2(c(seq(1, 0.2, length.out = 128), seq(0.25, 0, length.out = 128)), col)

map <- Heatmap(mtrx, col = col_fun, cluster_columns = TRUE, cluster_rows = FALSE, column_names_rot = 45)
pdf("final/RMSE_plots_CNORode_final.pdf", width = 15, height = 9)
map
dev.off()

save.image("final/Rdata.RData")

####Scatter plot predicted vs measured####
scatPredExp <- RMSEdf

cor.test(scatPredExp$measuredValues, scatPredExp$simulatedValues, method = "pearson", conf.level = 0.95)
cor(scatPredExp$measuredValues, scatPredExp$simulatedValues) ^ 2

p <- ggplot(scatPredExp, aes(x = measuredValues, y = simulatedValues)) + geom_polygon(data = shapeDF, mapping = aes(x  = x, y = y ), alpha = 0.05, color= NA) + geom_point() + geom_smooth(method="lm", se = FALSE, color = "grey", size = 3)
p <- ggplot(scatPredExp, aes(x = measuredValues, y = simulatedValues)) + geom_point() + geom_smooth(method="lm", se = FALSE, color = "grey", size = 3)
p + theme_bw() + theme(text = element_text(size = 20), panel.grid = element_blank()) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey", size = 1.2) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey", size = 1.2) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 1)) +
  xlab("Normalized kinase activity") +
  ylab("Simulated kinase activity") +
  ggtitle(print(paste("R^2 = ", summary(simpleLm)$r.squared, ". r = ", cor.test(scatPredExp$measuredValues, scatPredExp$simulatedValues, method = "pearson", conf.level = 0.95)$estimate, sep="")))

ggsave("Linear_V1/Scatterplot_measured_predicted.pdf", width = 14, height = 9)


stop()
####Calculate lm + interpolate graph
simpleLm <- lm(simulatedValues ~ measuredValues, data = scatPredExp)
simpleLm
summary(simpleLm)$r.squared
summary(simpleLm)$adj.r.squared
summary(simpleLm)

#The formula is
#measuredValues = simulatedValues * 0.6858 + 0.1398 (0.095)
0.6858 * 1 + 0.1398 +0.095
#x = c(0, 0, 1, 1)
#y = c(0.057, 0.257, 0.7365, 0.9365)
shapeDF <- data.frame(x  = c(0, 1, 1, 0), y =c(0.0448, 0.7306, 0.9206, 0.2348))
print(paste("R^2 = ", summary(simpleLm)$r.squared, ". r = ", cor.test(scatPredExp$measuredValues, scatPredExp$simulatedValues, method = "pearson", conf.level = 0.95)$estimate, sep=""))

#Calculate average RMSE
mean(mtrx) #0.1681661
shapeDFR <- data.frame(x  = c(0, 1, 1, 0), y =c(0, 0.6447339, 0.9810661, 0.3138661))
