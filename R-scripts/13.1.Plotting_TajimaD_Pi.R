# Anik Dutta

# Script to plot TajimaD and Pi output from vcftools

library("plyr")
library("RColorBrewer")
library("ggplot2")
group = c("1958","1963","2747","2931","2932","3111","4107","4117","4330")
a<-read.table("1963.10kb.TajD.Tajima.D", header = T)
names(a)
data <- c()
for(i in 1:length(group)){
    data_temp <- read.table(paste(group[i],"_HSP3.1bp.Tajima.D",sep=""),header=TRUE,stringsAsFactors=F)
    data_temp <- cbind(data_temp,rep(paste(group[i],sep=""),dim(data_temp)[1]))
    colnames(data_temp) <- c("CHROM","BIN_START","N_SNPS","TajimaD","Group")
    data <- rbind(data,data_temp)
  }
data <- na.omit(data)
mypal <- brewer.pal(9,"Spectral")
mypal2 <- c(mypal[5],mypal[5],mypal[5],mypal[5],mypal[3],mypal[5],mypal[3],mypal[7],mypal[7])

pdf(file="Tajima_Pop.pdf")
ggplot(data,aes(x=reorder(factor(Group),TajimaD,FUN=mean),y=TajimaD)) +
  #geom_dotplot(binaxis = 'y', binwidth=1/20,stackdir = 'center', position = position_dodge(), aes(fill=factor(Group)))+
  geom_violin(aes(fill = factor(Group)))+
  stat_summary(fun=mean, geom="point", size=2, color="black") + 
  scale_fill_manual(values=mypal2) + #ylab(expression(pi)) + 
  theme(axis.title.y = element_text(size = rel(1.8)))+ theme(legend.position="none") +
  xlab("Population") +ylab("Tajima's D") + theme(axis.title = element_text(size = 16),axis.text.x = element_text(size = 12, angle = 45, hjust = 1),axis.text.y = element_text(size = 12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(),panel.border = element_rect(color = "black",
                                                                     fill = NA,
                                                                     size = 1))
dev.off()

x = reorder(Population, Singleton, FUN = mean)
logpi_mean <- aggregate(log(data$PI,base=10),list(data$Group),mean, na.rm=T)
logpi_se <- aggregate(log(data$PI,base=10),list(data$Group),sd,na.rm=T)
ci_logpi <- logpi_se$x/221157*1.96
logpi_summary <- cbind(logpi_mean$x,logpi_mean$x - ci_logpi,logpi_mean$x + ci_logpi)
rownames(logpi_summary) <- logpi_mean$Group.1
colnames(logpi_summary) <- c(paste("mean log",expression(TajimaD)),"L_CI 95%","H_CI 95%")
write.table(logpi_summary,file="logpi.summary.group.tab",quote=FALSE,row.names=TRUE,col.names=TRUE)

## This code can estimate average Fst from pixy data
pi_mean <- aggregate(data$TajimaD,list(data$Group),mean,na.rm=T)
pi_se <- aggregate(data$TajimaD,list(data$Group),sd,na.rm=T)
DXY_sd <- aggregate(filt_Haps_DXY$avg_dxy,list(filt_Haps_DXY$pop1,filt_Haps_DXY$pop2),sd,na.rm=T)
ci_pi <- pi_se$x/221157*1.96
pi_summary <- cbind(pi_mean$x,pi_mean$x - ci_pi,pi_mean$x + ci_pi)
rownames(pi_summary) <- pi_mean$Group.1
colnames(pi_summary) <- c(paste(expression(pi)),"L_CI 95%","H_CI 95%")
write.table(pi_summary,file="pi.summary.group.tab",quote=FALSE,row.names=TRUE,col.names=TRUE)

log.glm <- glm(PI ~ Group, family=gaussian, data=data)
summary(log.glm)
