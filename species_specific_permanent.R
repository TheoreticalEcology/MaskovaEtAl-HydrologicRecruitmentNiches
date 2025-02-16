# used libraries
library(FSA)
library(rcompanion)
library(fields)
library(reshape2)
library(ggplot2) 

# import the data
data <- read.table("data.csv", header=T, sep=",")
head(data)

# 4 treatments (Dry, Wet, Flooded10 and Flooded40; two drainage treatments are excluded)
# only species with more than 10 seedling in total
# sum up seedlings in each treatment -> result is in pocty_final
# sum up seedling in each treatment and divide them by total number of seedlings -> result is in pocty_final_st

d1 <- data[data$hydrological_variant!="DrainedJune" & data$hydrological_variant!="DrainedAugust",]

regime <- d1$hydrological_variant
site <- d1$area

regime <- as.factor(regime)
site = as.factor(site)

d1 <- d1[,-c(1,2)]
d1 <- d1[,-which(colSums(d1)<10)]

pocty_final <- NULL
for(i in 1:ncol(d1)){
    druh <- d1[,i]

    pocty <- aggregate(druh, list(regime), sum)

    pocty_st <- c(pocty[1,2],pocty[4,2], pocty[2,2], pocty[3,2])

    pocty_final <- rbind(pocty_final,pocty_st)
    rownames(pocty_final)[i] <-names(d1)[i]
}


pocty_final_st <- NULL
for(i in 1:ncol(d1)){
    druh <- d1[,i]

    pocty <- aggregate(druh, list(regime), sum)

    pocty_st <- c(pocty[1,2],pocty[4,2], pocty[2,2], pocty[3,2])/sum(pocty[,2])

    pocty_final_st <- rbind(pocty_final_st,pocty_st)
    rownames(pocty_final_st)[i] <-names(d1)[i]
}


#save new dataset
P <- data.frame(rownames(pocty_final_st), pocty_final_st)
rownames(P) <- NULL
colnames(P) <- c("Species","Dry","Wet","F10","F40")

write.table(P,"data1.csv", sep=",", row.names=F)


# models for individual species
pismena <- NULL
for(j in 1:ncol(d1)){
    druh <- d1[,j]
    site_p <- aggregate(druh, list(site), sum)
    d = data.frame(druh=druh,regime=regime, site=site, siteSum = site_p$x[as.integer(site)]
)
    kruskal_ph <- dunnTest(druh ~ regime, data = d, method="bonferroni")
    vysl <- kruskal_ph$res
    pp <- cldList(P.adj ~ Comparison, data = vysl, threshold = 0.05)$Letter    
    pism <- c(pp[1],pp[4],pp[2],pp[3])
    pismena <- rbind(pismena, pism)
}

# reordering species according to results for drawing the figure
pism_stejna <- NULL
pocty_stejna <- NULL

pism_1 <- NULL
pocty_1 <- NULL


pism_2 <- NULL
pocty_2 <- NULL


pism_ruzna <- NULL
pocty_ruzna <- NULL

for(i in 1:nrow(pismena)){

    if(length(unique(pismena[i,])) == 1){
        pism_stejna <- rbind(pism_stejna,pismena[i,])
        pocty_stejna <- rbind(pocty_stejna,as.data.frame(pocty_final)[i,])
        
        }else{
            if(sum(table(pismena[i,]) == 1) == 1){
                pism_1 <- rbind(pism_1,pismena[i,])
                pocty_1 <- rbind(pocty_1,as.data.frame(pocty_final)[i,])
            }else{
                if(sum(table(pismena[i,]) == 2) == 2){
                    pism_2 <- rbind(pism_2,pismena[i,])
                    pocty_2 <- rbind(pocty_2,as.data.frame(pocty_final)[i,])
                }else{
                    pism_ruzna <- rbind(pism_ruzna,pismena[i,])
                    pocty_ruzna <- rbind(pocty_ruzna,as.data.frame(pocty_final)[i,]) 
                }
            }
        }

}

nuly <- which(pocty_2[,1]==0)

pocty_2final <- rbind(pocty_2[-nuly,],pocty_2[nuly,])
pism_2final <- rbind(pism_2[-nuly,],pism_2[nuly,])

pocty_2final <- rbind(pocty_ruzna[3,],pocty_2final[1:2,],pocty_ruzna[4,],pocty_2final[3:9,],pocty_ruzna[8,],pocty_2final[10:11,],pocty_ruzna[11,],pocty_2final[12:13,],pocty_ruzna[17,],pocty_2final[14:16,],pocty_2final[18:25,])
pism_2final <- rbind(pism_ruzna[3,],pism_2final[1:2,],pism_ruzna[4,],pism_2final[3:9,],pism_ruzna[8,],pism_2final[10:11,],pism_ruzna[11,],pism_2final[12:13,],pism_ruzna[17,],pism_2final[14:16,],pism_2final[18:25,])


pocty_1final <- rbind(pocty_1[-c(1,2,4,7,8,12),], pocty_1[c(7,12),], pocty_1[c(8),])
pism_1final <- rbind(pism_1[-c(1,2,4,7,8,12),], pism_1[c(7,12),], pism_1[c(8),])

pocty_stejna <- pocty_stejna[-13,]
pism_stejna <- pism_stejna[-13,]

pocty_ruzna <- pocty_ruzna[-c(3,4,8,11,17),]
pism_ruzna <- pism_ruzna[-c(3,4,8,11,17),]


prazdny_radek <- c(0,0,0,0)

pr_radek <- c(NA,NA,NA,NA)

Pocty <- data.matrix(rbind(pocty_1final, pocty_2final, pocty_stejna, pocty_ruzna))

Pismena <- rbind(pism_1final, pism_2final, pism_stejna, pism_ruzna)

# preparing data for drawing heatmap in gglot
druhy_poradi <- rev(rownames(Pocty))
druhy_poradi <- chartr("\\.","  ",druhy_poradi)

Pocty <- cbind(Pocty[,4],Pocty[,3],Pocty[,2],Pocty[,1])

colnames(Pocty) <- paste("Col",1:4)
rownames(Pocty) <- paste("Row",1:66)

P <- c(rev(Pismena[,1]),rev(Pismena[,2]),rev(Pismena[,3]),rev(Pismena[,4]))

df <- melt(Pocty)
colnames(df) <- c("x", "y", "value")

my.lines.stejne <- data.frame(x=c(0.5,0.5,0.5,4.5), y=c(15.5,27.5,15.5,15.5), 
    xend=c(4.5,4.5,0.5,4.5), yend=c(15.5,27.5,27.5,27.5))

my.lines.2F <- data.frame(x=c(2.5,2.5,2.5,4.5), y=c(27.5,35.5,27.5,27.5), 
    xend=c(4.5,4.5,2.5,4.5), yend=c(27.5,35.5,35.5,35.5))

my.lines.2D <- data.frame(x=c(0.5,0.5,0.5,2.5), y=c(35.5,56.5,35.5,35.5), 
    xend=c(2.5,2.5,0.5,2.5), yend=c(35.5,56.5,56.5,56.5))

my.lines.10 <- data.frame(x=c(2.5,2.5,2.5,3.5), y=c(56.5,57.5,56.5,56.5), 
    xend=c(3.5,3.5,2.5,3.5), yend=c(56.5,57.5,57.5,57.5))

my.lines.M <- data.frame(x=c(1.5,1.5,1.5,2.5), y=c(57.5,59.5,57.5,57.5), 
    xend=c(2.5,2.5,1.5,2.5), yend=c(57.5,59.5,59.5,59.5))

my.lines.D <- data.frame(x=c(0.5,0.5,0.5,1.5), y=c(59.5,66.5,59.5,59.5), 
    xend=c(1.5,1.5,0.5,1.5), yend=c(59.5,66.5,66.5,66.5))

# drawing the figure
p <- ggplot(df, aes(x = y, y = x, fill = rev(log(value+1)))) +
  geom_tile() + scale_fill_gradient(low="White",high="#2171B5", breaks = c(0,2.5,5,7.5),labels= c("0","15","150","1500")) +
  guides(fill = guide_colorbar(title = "Number of \n seedlings")) +
  geom_text(aes(label= P), size = 3) +
  geom_segment(data=my.lines.stejne, aes(x,y,xend=xend,yend=yend), size=0.5, inherit.aes=F) +
  geom_segment(data=my.lines.2F, aes(x,y,xend=xend,yend=yend), size=0.5, inherit.aes=F) +
  geom_segment(data=my.lines.2D, aes(x,y,xend=xend,yend=yend), size=0.5, inherit.aes=F) +
  geom_segment(data=my.lines.10, aes(x,y,xend=xend,yend=yend), size=0.5, inherit.aes=F) +
  geom_segment(data=my.lines.M, aes(x,y,xend=xend,yend=yend), size=0.5, inherit.aes=F) +
  geom_segment(data=my.lines.D, aes(x,y,xend=xend,yend=yend), size=0.5, inherit.aes=F)

p <- p + scale_x_discrete(label = c("Dry to \n moist","Wet","Flooded \n 10cm","Flooded \n 40cm"))

p <- p + scale_y_discrete(label = druhy_poradi)

p <- p + labs(x="Permanent hydrological treatments", y="Species")

png(filename="Heatmap_water_regime.png",width=14,height=21,
	units="cm",pointsize=10,res=120)
p
dev.off()
