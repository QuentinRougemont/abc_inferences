source('00-scripts/rscript/cv4abc.R')
if("data.table" %in% rownames(installed.packages()) == FALSE) 
{install.packages("qqman", repos="https://cloud.r-project.org") 
    print("installing packages data.table..." ) }

library(data.table)

target=as.numeric(read.table("00.data/OBS.ABC.stat.txt",skip=2,h=F))

#load simulations
IM1=fread("zcat 00.data/im.homom.homon.ABC.stat.txt.gz"    )
IM2=fread("zcat 00.data/im.heterom.homon.ABC.stat.txt.gz"  )
IM3=fread("zcat 00.data/im.homom.heteron.ABC.stat.txt.gz"  )
IM4=fread("zcat 00.data/im.heterom.heteron.ABC.stat.txt.gz")
AM1=fread("zcat 00.data/am.homom.homon.ABC.stat.txt.gz"    )
AM2=fread("zcat 00.data/am.heterom.homon.ABC.stat.txt.gz"  )
AM3=fread("zcat 00.data/am.homom.heteron.ABC.stat.txt.gz"  )
AM4=fread("zcat 00.data/am.heterom.heteron.ABC.stat.txt.gz")
SC1=fread("zcat 00.data/sc.homom.homon.ABC.stat.txt.gz"    )
SC2=fread("zcat 00.data/sc.heterom.homon.ABC.stat.txt.gz"  )
SC3=fread("zcat 00.data/sc.homom.heteron.ABC.stat.txt.gz"  )
SC4=fread("zcat 00.data/sc.heterom.heteron.ABC.stat.txt.gz")
SI1=fread("zcat 00.data/si.homom.homon.ABC.stat.txt.gz"    )
SI3=fread("zcat 00.data/si.homom.heteron.ABC.stat.txt.gz"  )

f <- function(x){
    m <- mean(x, na.rm = TRUE)
    x[is.na(x)] <- m
    x
}

SI1 <- apply(SI1, 2, f)
SI3 <- apply(SI3, 2, f)
SC1 <- apply(SC1, 2, f)
SC2 <- apply(SC2, 2, f)
SC3 <- apply(SC3, 2, f)
SC4 <- apply(SC4, 2, f)
IM1 <- apply(IM1, 2, f)
IM2 <- apply(IM2, 2, f)
IM3 <- apply(IM3, 2, f)
IM4 <- apply(IM4, 2, f)
AM1 <- apply(AM1, 2, f)
AM2 <- apply(AM2, 2, f)
AM3 <- apply(AM3, 2, f)
AM4 <- apply(AM4, 2, f)

nlinesFul=min(nrow(IM1), nrow(IM2), nrow(IM3), nrow(IM4),
              nrow(SI1), nrow(SI3),
              nrow(AM1), nrow(AM2), nrow(AM3), nrow(AM4),
              nrow(SC1), nrow(SC2), nrow(SC3), nrow(SC4))

IM1b=IM1[c(1:nlinesFul),-c(1:3,12,13,18:25)]
IM2b=IM2[c(1:nlinesFul),-c(1:3,12,13,18:25)]
IM3b=IM3[c(1:nlinesFul),-c(1:3,12,13,18:25)]
IM4b=IM4[c(1:nlinesFul),-c(1:3,12,13,18:25)]
SC1b=SC1[c(1:nlinesFul),-c(1:3,12,13,18:25)]
SC2b=SC2[c(1:nlinesFul),-c(1:3,12,13,18:25)]
SC3b=SC3[c(1:nlinesFul),-c(1:3,12,13,18:25)]
SC4b=SC4[c(1:nlinesFul),-c(1:3,12,13,18:25)]
AM1b=AM1[c(1:nlinesFul),-c(1:3,12,13,18:25)]
AM2b=AM2[c(1:nlinesFul),-c(1:3,12,13,18:25)]
AM3b=AM3[c(1:nlinesFul),-c(1:3,12,13,18:25)]
AM4b=AM4[c(1:nlinesFul),-c(1:3,12,13,18:25)]
SI1b=SI1[c(1:nlinesFul),-c(1:3,12,13,18:25)]
SI3b=SI3[c(1:nlinesFul),-c(1:3,12,13,18:25)]

target=target[-c(1:3,12,13,18:25)]
obs=matrix(rep(target[1:19],20), byrow=T, nrow=20)

#Comparaison All:
x=as.factor(c(rep(1:14,each=nlinesFul)))

z=rbind(SI1b,SI3b,
        IM1b,IM2b,IM3b,IM4b,
        AM1b,AM2b,AM3b,AM4b,
        SC1b,SC2b,SC3b,SC4b)

res5=model_selection_abc_nnet(target=obs, 
                              x=x, 
                              sumstat=z, 
                              tol=0.00025,
                              noweight=F,
                              rejmethod=F,
                              nb.nnet=50,
                              size.nnet=15,
                              output="ALL.MODELS.29.06.16")

mean_all=apply(res5,2,mean)
sd_all=apply(res5,2,sd)
write.table(mean_all,
            "mean_all",
            quote=F,
            row.names=c("SI1","SI2",
                        "IM1","IM2","IM3","IM4",
                        "AM1","AM2","AM3","AM4",
                        "SC1","SC2","SC3","SC4"),
            col.names=F)
