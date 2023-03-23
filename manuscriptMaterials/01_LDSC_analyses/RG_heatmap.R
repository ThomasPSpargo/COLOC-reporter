library(ggplot2)
library(dplyr)

args<- commandArgs(trailingOnly = T)
path <- dirname(args[1])
setwd(path)

scor = read.table(args[1],header=T)               # read in
#scor = scor[,c("p1","p2","rg","p")]               # retain key headers
scor$p1 = gsub("_munge.sumstats.gz","",scor$p1)   # assuming the munged files have format [phenotype]_munge.sumstats.gz
scor$p2 = gsub("_munge.sumstats.gz","",scor$p2)   # (adapt as necessary)

scor$p1 = gsub("\\./","",scor$p1) #Drop trait prefixes in p1 and p2
scor$p2 = gsub("\\./","",scor$p2)

#Generate an empty dataset for plot diagonal
traits <- unique(c(scor$p1,scor$p2)) 
newscor<- as.data.frame(matrix(NA,nrow=length(traits),ncol=ncol(scor),dimnames=list(NULL,names(scor))))
newscor$p1 <- newscor$p2 <- traits

#Create entries with flipped trait order to generate symmetrical matrix
scor.flip <- scor
colnames(scor.flip)[1:2] <- c("p2","p1")

#Append additional pairings
scor <- rbind(scor,newscor,scor.flip) #Append to scor

#Prepare plot labels
scor <- scor %>%            
  mutate(p1=as.factor(p1),
         p2=as.factor(p2),
         diag=if_else(p1==p2,as.character(p1),""),
         rg_label=if_else(as.numeric(p1)<=as.numeric(p2),NA_character_,
                          paste0("atop(r[g]=='",formatC(rg,format = "f",digits=2)," (",formatC(se,format = "f",digits=2),");',p=='",formatC(p,format = "f",digits=3),"')"),
         ),
         rg=if_else(as.numeric(p1)<=as.numeric(p2),NA_real_,rg) #Drop the 'upper tri' for the figure
  )


rgheat<- ggplot(data=scor,aes(x=p1,y=p2,fill=rg))+
  geom_tile()+
  geom_tile(data=scor[!is.na(scor$rg),],colour="black")+ #Generate borders specifically across the plotted datapoints
  geom_text(aes(label=diag))+
  geom_text(aes(label = rg_label), color = "black", size = 3,na.rm = TRUE,parse = TRUE) +
  scale_fill_distiller(palette = "RdBu", guide="colourbar",
                       limit = c(-1,1),
                          na.value = "white",
                       name=bquote(r[g]))+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = c(0.14,0.55),
        panel.grid = element_blank(),
        panel.border = element_blank()
  )#+

ggsave(file.path(path,"heatmap_global_rg.pdf")
       , rgheat, device = "pdf",
       units="mm",width=150,height=100)
