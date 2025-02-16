# used libraries
library(ape)
library(geiger)
library(nlme)
library(phytools)
library(DHARMa)
library(ade4)
library(multcomp)
library(FSA)
library(rcompanion)
library(fields)
library(reshape2)
library(ggplot2)
library(patchwork)

#=================================================================================================
# import table with germination data and table with trait data
data <- read.table("data1.csv", header=T, sep=",")
traity <- read.table("traity.csv", header=T, sep=",")

# add data about dormancy, oxygen requirements and alternating temperature requirement to the table with germination data
dorm <- NULL
for (radek in (1:nrow(data))){
    dorm <- rbind(dorm,unlist(traity[traity$Species==data[radek,1],][2:4]))
}
data <- cbind(data,dorm)

# replace a dot in the species name with a dash to match the Daphne format
druhy <- data[,1]
for (d in (1:length(druhy))){
    a <- unlist(strsplit(druhy[d], split="\\."))
    aa <- paste(a[1], a[2], sep="_")
    druhy[d] <- aa
}
data[,1] <- druhy

# make a new structure of data table so that each water treatment in one line

water <- colnames(data)[-1]
data_new <- NULL

for (d in (1:nrow(data))){
    for (w in 1:4){

        radek <- c(data[d,1],water[w],data[d,1+w], data[d,6], data[d,7], data[d,8])
        data_new <- rbind(data_new,radek)
    }
}

kliceni <- as.numeric(data_new[,3])
data_new <- data.frame(data_new[,-3],kliceni)
rownames(data_new) <- NULL
colnames(data_new) <- c("species","water","dorm","temp","ox","germination")

# merge PD, PY and MPD type of dormancy into one category
data_new$dorm[data_new$dorm=="PD"] <- "yes"
data_new$dorm[data_new$dorm=="PY"] <- "yes"
data_new$dorm[data_new$dorm=="MPD"] <- "yes"

# small correction of names for two species
data_new$species[data_new$species=="Alisma_plantago"] <- "Alisma_plantago-aquatica"
data_new$species[data_new$species=="Polygonum_persicaria"] <- "Persicaria_maculosa"

#=================================================================================================
# linking germination data with data from Daphne
tr=read.tree("DaPhnE_01.tre")
taxa=tr$tip.label

# find matches in original data (to find species which are written incorectly of are missing in Daphne)
spc = as.character(data_new[,1])
firstmatch = match(spc,taxa)
chybejici <- spc[is.na(firstmatch)]

# make a tree 
todel = which(is.na(match(taxa,spc))) 
subtree = drop.tip(tr,todel)
taxa_subtree=subtree$tip.label

# plot a tree
plot(subtree,cex=0.5)
axisPhylo()

# vezme z tabulky s daty pouze ty, ktere jsou v Daphne (9 species excluded - mostly those which we were not able to distinguis into the species level and therefore were written in the original data as "Carex sp.", "Chara sp.",...)
dta <- data_new[!data_new$species %in% chybejici,]

#=================================================================================================
# models
# checking phylogenetic structure according to: https://theoreticalecology.github.io/AdvancedRegressionModels/4C-CorrelationStructures.html#phylogenetic-structures-pgls

# for dormancy
data_interD <- with(dta, interaction(dta$water, dta$dorm, sep = "x"))
modelD <- glm(dta$germination ~ data_interD)
res <- simulateResiduals(modelD)
res2 <- recalculateResiduals(res, group=subtree$tip.label)
testSpatialAutocorrelation(res2, distMat = cophenetic(subtree))

summary(modelD)
confint(modelD)
l2 <- glht(modelD, linfct = mcp(data_interD = "Tukey"))
pismenaD <- cld(l2)

# for alternating temperature requirement
data_interT <- with(dta, interaction(dta$water, dta$temp, sep = "x"))
modelT <- glm(dta$germination ~ data_interT)
res <- simulateResiduals(modelT)
res2 <- recalculateResiduals(res, group=subtree$tip.label)
testSpatialAutocorrelation(res2, distMat = cophenetic(subtree))

summary(modelT)
confint(modelT)
l2 <- glht(modelT, linfct = mcp(data_interT = "Tukey"))
pismenaT <- cld(l2)

# for oxygen requirement
data_interO <- with(dta, interaction(dta$water, dta$ox, sep = "x"))
modelO <- glm(dta$germination ~ data_interO)
res <- simulateResiduals(modelO)
res2 <- recalculateResiduals(res, group=subtree$tip.label)
testSpatialAutocorrelation(res2, distMat = cophenetic(subtree))

summary(modelO)
confint(modelO)
l2 <- glht(modelO, linfct = mcp(data_interO = "Tukey"))
pismenaO <- cld(l2)

#=================================================================================================
# Drawing the figure
# Dormancy
dta$water <- factor(dta$water, levels = c("Dry_count","Wet_count","F10_count","F40_count"))
p <- ggplot(na.omit(dta[,-(c(4,5))]), aes(x=dorm, y = log(germination), fill=water)) + 
geom_boxplot() + scale_fill_brewer(palette="Blues") + 
theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) +
theme(legend.position = "none") +
theme(panel.background = element_rect(fill = "white")) +
theme(panel.border = element_rect(fill = "transparent", color = "black", size = 1))

p <- p + labs(x="Dormancy", y="Germination [%]")

labels <- c("No (N=28)","Yes (N=29)")

p <- p + scale_x_discrete(label = labels)

pism <- c(pismenaD[[10]]$Letters[1],pismenaD[[10]]$Letters[4],pismenaD[[10]]$Letters[2],pismenaD[[10]]$Letters[3],pismenaD[[10]]$Letters[5],pismenaD[[10]]$Letters[8],pismenaD[[10]]$Letters[6],pismenaD[[10]]$Letters[7])

pozice <- c(0.7,0.9,1.1,1.3, 1.7,1.9,2.1,2.3)

for (i in 1:8){
    p <- p + geom_text(x = pozice[i], y = 1.02, label = pism[i], size=4)
}

pD <- p

# Temperature
dta$water <- factor(dta$water, levels = c("Dry","Wet","F10","F40"))
p <- ggplot(na.omit(dta[,-(c(3,5))]), aes(x=temp, y = germination, fill=water)) + 
geom_boxplot() + scale_fill_brewer(palette="Blues") + 
theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) +
theme(legend.position = "none") + 
theme(panel.background = element_rect(fill = "white")) +
theme(panel.border = element_rect(fill = "transparent", color = "black", size = 1))

p <- p + labs(x="Temperature requirements", y="")

labels <- c("No (N=22)","Yes (N=26)")

p <- p + scale_x_discrete(label = labels)

pism <- c(pismenaT[[10]]$Letters[1],pismenaT[[10]]$Letters[4],pismenaT[[10]]$Letters[2],pismenaT[[10]]$Letters[3],pismenaT[[10]]$Letters[5],pismenaT[[10]]$Letters[8],pismenaT[[10]]$Letters[6],pismenaT[[10]]$Letters[7])

pozice <- c(0.7,0.9,1.1,1.3, 1.7,1.9,2.1,2.3)

for (i in 1:8){
    p <- p + geom_text(x = pozice[i], y = 1.02, label = pism[i], size=4)
}

pT <- p

# Oxygen
dta$water <- factor(dta$water, levels = c("Dry","Wet","F10","F40"))
dta$ox <- factor(dta$ox, levels = c("HY","OX","both"))


p <- ggplot(na.omit(dta[,-(c(3,4))]), aes(x=ox, y = germination, fill=water)) + 
geom_boxplot() + scale_fill_brewer(palette="Blues", labels = c("Dry to \n moist", "Wet", "Flooded \n 10cm", "Flooded \n 40cm")) + 
theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + 
guides(fill=guide_legend(title="Hydrological \n treatment")) +
theme(panel.background = element_rect(fill = "white")) +
theme(panel.border = element_rect(fill = "transparent", color = "black", size = 1))

p <- p + labs(x="Oxygen requirements", y="Germination [%]")

labels <- c("Hypoxic (N=6)", "Oxic (N=15)","Both (N=6)")

p <- p + scale_x_discrete(label = labels)

pism <- c(pismenaO[[10]]$Letters[5],pismenaO[[10]]$Letters[8],pismenaO[[10]]$Letters[6],pismenaO[[10]]$Letters[7],pismenaO[[10]]$Letters[9],pismenaO[[10]]$Letters[12],pismenaO[[10]]$Letters[10],pismenaO[[10]]$Letters[11],pismenaO[[10]]$Letters[1],pismenaO[[10]]$Letters[4],pismenaO[[10]]$Letters[2],pismenaO[[10]]$Letters[3])

pozice <- c(0.7, 0.9, 1.1, 1.3, 1.7, 1.9, 2.1, 2.3, 2.7, 2.9, 3.1, 3.3)

for (i in 1:12){
    p <- p + geom_text(x = pozice[i], y = 1.02, label = pism[i], size=4)
}

pO <- p

png(filename="Traits_boxplots.png",width=21,height=14,
	units="cm",pointsize=10,res=120)

(pD | pT) / pO

dev.off()
