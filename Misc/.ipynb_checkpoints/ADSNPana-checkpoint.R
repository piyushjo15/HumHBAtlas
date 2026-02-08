library(ggplot2)
setwd("~/DeepHB/Outs/")
df <- read.delim("AD_SNP/pred_AD_lstmv2_57MP_73.txt",header = FALSE) ##predicted similar to fetal brain
#df <- read.delim("AD_SNP/pred_AD_cnn_194.txt",header = FALSE)

head(df)

beds <- read.delim("AD_SNP/CREs_AD_SNP.txt")
cres <- row.names(beds) <- beds$X
head(beds)

## MP classes
mp_cls <- readLines("~/DeepHB/TopicMPs.txt")
colnames(df) <- mp_cls
row.names(df) <- row.names(beds)
df[1:4,1:4]

thres_mp <- read.delim("MP_threshold_summary_deephb_lstmv2_57MP_73.txt")

row.names(thres_mp) <- thres_mp$MP
thres_mp <- thres_mp$Best.Threshold
names(thres_mp) <- mp_cls

beds$ID <- NA
beds$Score <- 0

## create a matrix of 0 and 1
df2 <- matrix(0,dim(df)[1],dim(df)[2])
colnames(df2)<- colnames(df)
row.names(df2)<- row.names(df)
for (i in 1:ncol(df)) {
  df2[, i] <- ifelse(df[, i] >= thres_mp[i], 1,0)
}

rs <- rowSums(df2)
# ## find which cres didn't get assigned any id and remove them
# keep <- rs==0
# cres_rem <- cres[keep]
## which CREs have only one prediction
keep <- rs==1
cres_1c<- cres[keep]
sel_ind <- max.col(df2[cres_1c,])
prob <- apply(df[cres_1c,],1,max)
beds[cres_1c,"ID"] <-mp_cls[sel_ind]
beds[cres_1c,"Score"] <-prob

## now for CREs with more than one prediction,
## find the prediction that passes the cut-off and is max
keep <- rs>1
cres_mc<- cres[keep]

## make below cut-off as 0
df2x <- df2[cres_mc,]
dfx <- df[cres_mc,]
keep <-df2x==0
dfx[keep]<- 0
sel_ind <- max.col(dfx)
prob <- apply(dfx,1,max)

beds[cres_mc,"ID"] <-mp_cls[sel_ind]
beds[cres_mc,"Score"] <-prob


ll <- beds$ID

mdf <- data.frame(table(ll))
mdf <- mdf[order(-mdf$Freq), ]
mdf$cum_freq <- cumsum(mdf$Freq)
mdf$cum_percent <- 100 * mdf$cum_freq / sum(mdf$Freq)
save(mdf,beds,ll, file = "AD_SNP/Prediction_lstmv2_57MP_73.RData")
load("AD_SNP/Prediction_lstmv2_57MP_73.RData")

p <- ggplot(mdf, aes(x = reorder(ll, -Freq))) +
  geom_bar(aes(y = Freq), stat = "identity", fill = "steelblue") +
  geom_line(aes(y = cum_percent * max(Freq) / 100, group = 1), color = "red", linewidth = 1) +
  geom_point(aes(y = cum_percent * max(Freq) / 100), color = "red", size = 2) +
  theme_bw()+
  theme(axis.line = element_line(colour = 'black',linewidth=.5),
        axis.ticks = element_line(colour = 'black',linewidth=.5),
        axis.text.x = element_text(face="bold",colour = "black",size=8, angle = 90, vjust = 0.5,hjust=1),
        axis.text.y = element_text(face="bold",colour = "black",size=10),
        axis.title=element_blank(),
        legend.position = "None")+
  scale_y_continuous(
    name = "Frequency",
    sec.axis = sec_axis(~ . * 100 / max(mdf$Freq), name = "Cumulative %")
  ) +
  labs(title = "Pareto Plot", x = "Name") 

pdf("AD_SNP/Prediction_lstmv2_57MP_73.pdf",  width=8, height=3, pointsize = 10)
print(p)
dev.off()

##add rsid
md <- read.delim("AD_SNP/AD_SNP_aug.bed",header = FALSE) ##augment bed
row.names(md)<- paste0(md$V1,":",md$V2,"-",md$V3)
table(row.names(beds) %in% row.names(md))
beds$Peak <- md[row.names(beds),"V4"]


## now i want to find which SNP is associated which class and which peak region
## is representative. so i want to first find what is the SNP id for each of the
## peak. so all the augmented peaks for a SNP. then I want to find which particular
## class  appeared the most and then select the peak with maximum score that class.
## if there is draw between classes, i will chose the one with max score.
beds2 <- beds %>% separate(Peak,c("RSID","Adj"),sep="_")
beds2[is.na(beds2$ID),"ID"] <- "ND"
beds2$Adj <- as.numeric(beds2$Adj)
beds2$Peak <- beds$Peak
head(beds2)
beds2 <- beds2[beds2$ID!="ND",]
rsids <- unique(beds2$RSID)
sel_df <- c()

for(x in rsids ){
  del <- beds2[beds2$RSID==x,]
  if(dim(del)[1]==1){
    del_df <- data.frame(RSID=x,MP=del$ID,Score=del$Score,
                         Peak=del$Peak,Adj=del$Adj)
  } else {
    mpss <- del$ID
    scrs <- del$Score
    tbs <- data.frame(table(mpss))
    mpss_sel <- as.character(tbs$mpss)[tbs$Freq==max(tbs$Freq)]
    keep <- mpss %in% mpss_sel
    mpss <- mpss[keep]
    scrs <- scrs[keep]
    id <- which(scrs==max(scrs))
    ## if the length of id is more than 1, then I want to select the one which is
    ## close to center, if both are same distance just select on positive side
    if(length(id)>1){
      deldel <- del[id,]
      deldel$Adj2 <- abs(deldel$Adj-50)
      id <- id[which(deldel$Adj2==min(deldel$Adj2))]
      rm(deldel)
    }
    mpss <- mpss[id]
    scrs <- scrs[id]
    
    del_df <- data.frame(RSID=x,MP=mpss,Score=scrs,
                         Peak=del[id,"Peak"],Adj=del[id,"Adj"])
    rm(mpss,keep,scrs,tbs,mpss_sel,id)
    
  }
  
  sel_df <- rbind(sel_df,del_df)
  rm(del,del_df)
}

head(sel_df)
table(duplicated(sel_df$RSID))
ll <- sel_df$MP

mdf <- data.frame(table(ll))
mdf <- mdf[order(-mdf$Freq), ]
mdf$cum_freq <- cumsum(mdf$Freq)
mdf$cum_percent <- 100 * mdf$cum_freq / sum(mdf$Freq)
save(sel_df,ll,beds2, file = "AD_SNP/Prediction_lstmv2_57MP_73_v2.RData")

p <- ggplot(mdf, aes(x = reorder(ll, -Freq))) +
  geom_bar(aes(y = Freq), stat = "identity", fill = "steelblue") +
  geom_line(aes(y = cum_percent * max(Freq) / 100, group = 1), color = "red", linewidth = 1) +
  geom_point(aes(y = cum_percent * max(Freq) / 100), color = "red", size = 2) +
  theme_bw()+
  theme(axis.line = element_line(colour = 'black',linewidth=.5),
        axis.ticks = element_line(colour = 'black',linewidth=.5),
        axis.text.x = element_text(face="bold",colour = "black",size=8, angle = 90, vjust = 0.5,hjust=1),
        axis.text.y = element_text(face="bold",colour = "black",size=10),
        axis.title=element_blank(),
        legend.position = "None")+
  scale_y_continuous(
    name = "Frequency",
    sec.axis = sec_axis(~ . * 100 / max(mdf$Freq), name = "Cumulative %")
  ) +
  labs(title = "Pareto Plot", x = "Name") 


