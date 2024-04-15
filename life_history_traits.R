library(Boruta)
library(ggplot2)
library(glmnet)
library(caret)
library(foreach)
library(doParallel)

cores <- 2
registerDoParallel(cores)

lht.dt <- fread("./data/lht.csv.gz")
lht.dt <- lht.dt[,apply(lht.dt, 2, function(col) { length(unique(col)) > 1 }), with=F]
lht.dt[,org:=gsub("GCF","GCA",str_extract(org, "GC(A|F)?_[0-9]+\\.1"))] 
lht.dt <- lht.dt[order(org)]

meta.dat <- fread("./data/lht_meta-data.csv")

all.equal(lht.dt$org, meta.dat$org)

covariats <- lht.dt[,-"org"]
response.source  <- factor(meta.dat$Origin2)
response.growth  <- meta.dat$mumax_mean

idx.missing <- which(meta.dat$strain %in% c("PA14", "PAO1"))
response.virulence   <- factor(meta.dat$`in mouse`[-idx.missing], 
                               levels=c("very low","medium","very high"), 
                               labels=c("low", "medium", "high"), ordered=T)
covariats.virulence <- covariats[-idx.missing,]

# helper functions
printPerf <- function(caret.train, n.round=2){
	caret.train.perf <- caret.train$results[rownames(caret.train$bestTune),] # 
	perf <- apply(apply(caret.train.perf, 1, function(x) 
	  {n <- names(caret.train.perf); paste0(paste(n,round(as.numeric(x),n.round), sep = ":"))}), 
	  2, paste0, collapse = ",")
	return(paste(perf,"(LOOCV)"))
}
addNameDT <- function(feat.dt){
	lht.meta.dt <- fread("./data/lht_strategies.tsv")
  feat.dt[,name:=str_remove(rn, ".*(?<=_)")]
	feat.dt[!is.na(match(rn, lht.meta.dt$id)), name:=lht.meta.dt$name[match(rn, lht.meta.dt$id)]]
}
mycombine <- function(x, ...) { mapply(rbind,x,...,SIMPLIFY=FALSE) }
# filter selected features by stability (occurence in 50% of the runs)
fsel.stability <- function(fsel.dt, cutoff=0.5){
    min.n <- 1/cutoff
    if(nrow(fsel.dt)==0) return(NULL)
    fsel.dt[,.N,by=.(rn)][N>fsel.dt[,length(unique(i))]/min.n, rn] 
}

# feature selection function
feat.selection <- function(covariats, response, lasso.family, n=100, skip.lasso=F){
    boruta.fsel.dt <- data.table()
    lasso.fsel.dt  <- data.table()
    fsel.lst <- foreach(i=1:n, .combine="mycombine", .multicombine=TRUE) %dopar%{	
        boruta.fsel <- Boruta(x=covariats, y=response, maxRuns=1000)	
        boruta.fsel <- TentativeRoughFix(boruta.fsel)
        boruta.tmp  <- data.table(i, attStats(boruta.fsel), keep.rownames=T)[decision=="Confirmed"]
        if(!skip.lasso){
            lasso.cv    <- cv.glmnet(x=data.matrix(covariats), 
                                     y=response, alpha=1, family=lasso.family)
            lasso.fsel  <- glmnet(x=data.matrix(covariats), 
                                  y=response, alpha = 1, lambda = lasso.cv$lambda.min, 
                                  family=lasso.family)
            if(sum(coef(lasso.fsel)!=0)>1) 
              lasso.tmp <- data.table(i, 
                                      rn=names(coef(lasso.fsel)[which(coef(lasso.fsel)[,1]!=0),][-1]), 
                                      coef=unname(coef(lasso.fsel)[which(coef(lasso.fsel)[,1]!=0),][-1])) 
            else lasso.tmp <- data.table()
        } else lasso.tmp <- data.table()
        list(boruta.tmp, lasso.tmp)
    }
    return(fsel.lst)
}

# cross-validation
feat.cv <- function(boruta.feat, lasso.feat, covariats, response, n.tune=1000){
    caret.boruta <- train(x=covariats[,boruta.feat,with=F], y=response, 
                          method = "rf", 
                          trControl = trainControl(method="LOOCV", number = 10, search="random"), 
                          tuneLength=n.tune)
    cat("boruta:", printPerf(caret.boruta, 3),"\n")
    if(length(lasso.feat)>0){
        if(length(lasso.feat)>2) lin.method = "glmnet" else lin.method="lm"
        caret.lasso <- train(x=covariats[,lasso.feat,with=F], y=response, 
                             method = lin.method, 
                             trControl = trainControl(method="LOOCV", number = 10, search="random"), 
                             tuneLength=n.tune, trControl = trainControl(search = "grid"), 
                             tuneGrid = data.frame(alpha = 1, lambda = 2^runif(n.tune, min = -10, 3)))
        cat("lasso: ", printPerf(caret.lasso,  3),"\n")
    } else caret.lasso <- NULL
    return(list(caret.boruta, caret.lasso))
}


plot.boruta <- function(fsel.boruta, boruta.cv, covariats, response, boruta.feat){
    cor.dt <- covariats[,list(rn=colnames(covariats[,.SD,.SDcols=is.numeric]), 
                              cor=unlist(
                                lapply(.SD, function(x) 
                                  cor.test(x,as.numeric(response), method="spearman")$estimate))),.SDcols=is.numeric]
    boruta.dt <- merge(fsel.boruta, cor.dt, by="rn")
    g <- ggplot(addNameDT(boruta.dt[rn %in% boruta.feat])) + 
      geom_boxplot(aes(x=name, y=meanImp, fill=cor)) + 
      theme_minimal(base_size=14) + 
      coord_flip() + 
      xlab("") + 
      ylab("Importance of feature (random forest)") + 
      theme(axis.title.x = element_text(hjust=+0.5), 
            plot.caption = element_text(hjust = 0, face= "italic"), 
            plot.caption.position = "plot") + 
      labs(x=NULL, caption=printPerf(boruta.cv)) +  
      scale_x_discrete(labels = scales::label_wrap(40)) + 
      scale_fill_gradient2(midpoint=0, low="deepskyblue", mid="white", high="deeppink", space ="Lab" ) + 
      labs(color="Correlation")
    return(g)
}

plot.lasso <- function(fsel.lasso, lasso.cv, lasso.feat){
    g <- ggplot(addNameDT(fsel.lasso[rn %in% lasso.feat])) + 
      geom_boxplot(aes(x=name, y=coef)) + 
      theme_minimal(base_size=14) + 
      coord_flip() + 
      xlab("") + 
      ylab("Regression coefficient of feature (lasso)") + 
      theme(axis.title.x = element_text(hjust=+0.5), 
            plot.caption = element_text(hjust = 0, face= "italic"), 
            plot.caption.position = "plot") + 
      labs(x=NULL, caption=printPerf(lasso.cv)) + 
      scale_x_discrete(labels = scales::label_wrap(40)) + 
      geom_hline(yintercept=0, linetype="dashed", color = "red")
    return(g)
}


#################
# strain source #
#################
# compute
fsel.source      <- feat.selection(covariats, response.source, lasso.family="binomial")
cv.source        <- feat.cv(fsel.stability(fsel.source[[1]]) , 
                            fsel.stability(fsel.source[[2]]), 
                            covariats, 
                            response.source)
saveRDS(list(fsel.source, cv.source), "./data/hannover_lht-source.RDS")
# plot
source.tmp <- readRDS("./data/hannover_lht-source.RDS"); 
fsel.source <- source.tmp[[1]]; cv.source <- source.tmp[[2]]
plot.boruta(fsel.source[[1]], cv.source[[1]], covariats, response.source, fsel.stability(fsel.source[[1]])) + 
  scale_fill_gradient2(midpoint=0, low="#F8766D", mid="white", high="#00BFC4", space ="Lab" ) 
ggsave("./figures/hannover-lht_boruta-source.pdf", width=5.5, height=4)
plot.lasso(fsel.source[[2]], cv.source[[2]], fsel.stability(fsel.source[[2]]))
ggsave("./figures/hannover-lht_lasso-source.pdf", width=4.5, height=4)


#################
# growth        #
#################
# compute
fsel.growth      <- feat.selection(covariats, response.growth, lasso.family="gaussian")
feat.growth      <- c(fsel.stability(fsel.growth[[1]], cutoff=0.1), 
                      fsel.stability(fsel.growth[[2]], cutoff=0.1)) 
cv.growth        <- feat.cv(fsel.stability(fsel.growth[[1]], cutoff=0.1), 
                            fsel.stability(fsel.growth[[2]], cutoff=0.1), covariats, response.growth)
saveRDS(list(fsel.growth, cv.growth), "./data/hannover_lht-growth.RDS")
# plot
growth.tmp <- readRDS("./data/hannover_lht-growth.RDS"); 
fsel.growth <- growth.tmp[[1]]; 
cv.growth <- growth.tmp[[2]]
plot.boruta(fsel.growth[[1]], cv.growth[[1]], covariats, response.growth, fsel.stability(fsel.growth[[1]], cutoff=0.1))
ggsave("./figures/hannover-lht_boruta-growth.pdf", width=5, height=2.3)
plot.lasso(fsel.growth[[2]], cv.growth[[2]], fsel.stability(fsel.growth[[2]], cutoff=0.1))
# plot prediction vs experimental data
boruta.growth.pred.dt <- data.table(pred=predict(cv.growth[[1]]), 
                                    exp=response.growth, 
                                    org=meta.dat$strain, 
                                    source=meta.dat$Origin2)											
ggplot(boruta.growth.pred.dt, aes(x=exp, y=pred)) + 
  geom_abline(slope=1, intercept=0, linetype="dashed", color = "red") + geom_point(size=3) + 
  theme_minimal(base_size=14) + 
  coord_flip() + 
  xlab("Experimental growth rate") + 
  ylab("Predicted growth rate") + 
  ggrepel::geom_label_repel(aes(label=org, fill=source), show.legend = T) + 
  ggpubr::stat_cor(method = "spearman", label.x.npc=1, label.y.npc=0)
ggsave("./figures/hannover-lht_boruta-growth-prediction.pdf", width=5.5, height=4)


#################
# virulence     #
#################
# compute
fsel.virulence      <- feat.selection(covariats.virulence, 
                                      response.virulence, 
                                      lasso.family="multinomial", 
                                      skip.lasso=T)
cv.virulence        <- feat.cv(fsel.stability(fsel.virulence[[1]]) , 
                               fsel.stability(fsel.virulence[[2]]), 
                               covariats.virulence, 
                               response.virulence)
saveRDS(list(fsel.virulence, cv.virulence), "./data/hannover_lht-virulence.RDS")
# plot
virulence.tmp <- readRDS("./data/hannover_lht-virulence.RDS"); 
fsel.virulence <- virulence.tmp[[1]]; 
cv.virulence <- virulence.tmp[[2]]
plot.boruta(fsel.virulence[[1]], cv.virulence[[1]], covariats.virulence, response.virulence, fsel.stability(fsel.virulence[[1]]))
ggsave("./figures/hannover-lht_boruta-vir.pdf", width=4.5, height=1.7)
#plot.lasso(fsel.virulence[[2]], cv.virulence[[2]], fsel.stability(fsel.virulence[[2]]))
