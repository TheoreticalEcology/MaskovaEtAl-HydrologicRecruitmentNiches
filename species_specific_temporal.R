# used libraries
library(FSA)
library(rcompanion)
library(fields)
library(reshape2)
library(ggplot2)

# import the data
data <- read.table("data.csv", header=T, sep=",")
head(data)

# 3 treatments (Flooded40, DrainageJune and DrainageAugust, others are excluded)
# only species with more than 10 seedling in total
# sum up seedlings in each treatment -> result is in pocty_final2
# sum up seedling in each treatment and divide them by total number of seedlings -> result is in pocty_final2_st

d2 <- data[data$hydrological_variant=="Flooded40" | data$hydrological_variant=="DrainedJune" | data$hydrological_variant=="DrainedAugust",]

regime <- d2$hydrological_variant
site <- d2$area

regime <- as.factor(regime)
site = as.factor(site)

d2 <- d2[,-c(1,2)]
d2 <- d2[,-which(colSums(d2)<10)]

pocty_final2 <- NULL
for(i in 1:ncol(d2)){
    druh <- d2[,i]

    pocty <- aggregate(druh, list(regime), sum)

    pocty_st <- c(pocty[3,2],pocty[2,2], pocty[1,2])

    pocty_final2 <- rbind(pocty_final2,pocty_st)
    rownames(pocty_final2)[i] <-names(d2)[i]
}


pocty_final2_st <- NULL
for(i in 1:ncol(d2)){
    druh <- d2[,i]

    pocty <- aggregate(druh, list(regime), sum)

    pocty_st <- c(pocty[3,2],pocty[2,2], pocty[1,2])/sum(pocty[,2])

    pocty_final2_st <- rbind(pocty_final2_st,pocty_st)
    rownames(pocty_final2_st)[i] <-names(d2)[i]
}

#save new dataset
P <- data.frame(rownames(pocty_final2_st), pocty_final2_st)
rownames(P) <- NULL
colnames(P) <- c("Species","Flooded40","DrainageJ", "DrainageA")

write.table(P,"data2.csv", sep=",", row.names=F)


# models for individual species
pismena <- NULL

for(j in 1:ncol(d2)){
    druh <- d2[,j]
    site_p <- aggregate(druh, list(site), sum)
    d = data.frame(druh=druh,regime=regime, site=site, siteSum = site_p$x[as.integer(site)]
)
#   kruskal.test(druh ~ regime, data = d)

    kruskal_ph <- dunnTest(druh ~ regime, data = d, method="bonferroni")
    vysl <- kruskal_ph$res
    pp <- cldList(P.adj ~ Comparison, data = vysl, threshold = 0.05)$Letter    
    pism <- c(pp[3],pp[2],pp[1])
 
    pismena <- rbind(pismena, pism)
}

# reordering species according to results for drawing the figure
pism_stejna <- NULL
pocty_stejna <- NULL

pism_1 <- NULL
pocty_1 <- NULL

pism_ruzna <- NULL
pocty_ruzna <- NULL

for(i in 1:nrow(pismena)){

    if(length(unique(pismena[i,])) == 1){
        pism_stejna <- rbind(pism_stejna,pismena[i,])
        pocty_stejna <- rbind(pocty_stejna,as.data.frame(pocty_final2)[i,])
        
        }else{
            if(sum(table(pismena[i,]) == 1) == 1){
                pism_1 <- rbind(pism_1,pismena[i,])
                pocty_1 <- rbind(pocty_1,as.data.frame(pocty_final2)[i,])
            }else{
                pism_ruzna <- rbind(pism_ruzna,pismena[i,])
                pocty_ruzna <- rbind(pocty_ruzna,as.data.frame(pocty_final2)[i,]) 
            }
        }

}

nuly <- pocty_1[pism_1[,2]=="a",1]==0
nenuly <- pocty_1[pism_1[,2]=="a",1]!=0

pism_1final <- rbind(pism_1[pism_1[,2]=="b" & pism_1[,1]=="b",], pism_1[pism_1[,2]=="a",][nuly,], pism_1[pism_1[,2]=="a",][nenuly,], pism_1[pism_1[,1]=="a",])

pocty_1final <- rbind(pocty_1[pism_1[,2]=="b" & pism_1[,1]=="b",], pocty_1[pism_1[,2]=="a",][nuly,], pocty_1[pism_1[,2]=="a",][nenuly,], pocty_1[pism_1[,1]=="a",])

pocty_final <- rbind(pocty_1final[8,],pocty_1final[c(12:15,17:24),],pocty_1final[c(3:6,7,9,10),],pocty_1final[2,],pocty_stejna[-1,],pocty_ruzna[c(1,2,4,6:8,5),])

pism_final <- rbind(pism_1final[8,],pism_1final[c(12:15,17:24),],pism_1final[c(3:6,7,9,10),],pism_1final[2,],pism_stejna[-1,],pism_ruzna[c(1,2,4,6:8,5),])


prazdny_radek <- c(0,0,0)

pr_radek <- c(NA,NA,NA)

Pocty <- data.matrix(pocty_final)

Pismena <- pism_final

# preparing data for drawing heatmap in gglot
druhy_poradi <- rev(rownames(Pocty))
druhy_poradi <- chartr("\\.","  ",druhy_poradi)

Pocty <- cbind(Pocty[,3],Pocty[,2],Pocty[,1])

colnames(Pocty) <- paste("Col",1:3)
rownames(Pocty) <- paste("Row",1:42)

P <- c(rev(Pismena[,1]),rev(Pismena[,2]),rev(Pismena[,3]))

df <- melt(Pocty)
colnames(df) <- c("x", "y", "value")


my.lines.stejne <- data.frame(x=c(0.5,0.5,0.5,3.5), y=c(7.5,21.5,7.5,7.5), 
    xend=c(3.5,3.5,0.5,3.5), yend=c(7.5,21.5,21.5,21.5))

my.lines.2F <- data.frame(x=c(1.5,1.5,1.5,2.5), y=c(29.5,41.5,29.5,29.5), 
    xend=c(2.5,2.5,1.5,2.5), yend=c(29.5,41.5,41.5,41.5))

my.lines.2D <- data.frame(x=c(0.5,0.5,0.5,1.5), y=c(41.5,42.5,41.5,41.5), 
    xend=c(1.5,1.5,0.5,1.5), yend=c(41.5,42.5,42.5,42.5))

my.lines.10 <- data.frame(x=c(1.5,1.5,1.5,3.5), y=c(22.5,29.5,22.5,22.5), 
    xend=c(3.5,3.5,1.5,3.5), yend=c(22.5,29.5,29.5,29.5))

# drawing the figure
p <- ggplot(df, aes(x = y, y = x, fill = rev(log(value+1)))) +
  geom_tile() + scale_fill_gradient(low="White",high="#2171B5", breaks = c(0,2.5,5,7.5),labels= c("0","15","150","1500")) +
  guides(fill = guide_colorbar(title = "Number of \n seedlings")) +
  geom_text(aes(label= P), size = 3) +
  geom_segment(data=my.lines.stejne, aes(x,y,xend=xend,yend=yend), size=0.5, inherit.aes=F) +
  geom_segment(data=my.lines.2F, aes(x,y,xend=xend,yend=yend), size=0.5, inherit.aes=F) +
  geom_segment(data=my.lines.2D, aes(x,y,xend=xend,yend=yend), size=0.5, inherit.aes=F) +
  geom_segment(data=my.lines.10, aes(x,y,xend=xend,yend=yend), size=0.5, inherit.aes=F) 

p <- p + scale_x_discrete(label = c("Flooded \n 40cm","Drainage \n in June","Drainage \n in August"))

p <- p + scale_y_discrete(label = druhy_poradi)

p <- p + labs(x="Hydrological treatments", y="Species")

png(filename="Heatmap_drainage.png",width=14,height=21,
	units="cm",pointsize=10,res=120)
p
dev.off()
