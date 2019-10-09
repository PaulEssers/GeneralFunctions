norm_01<-function(x){(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))}

# needed to rename signature names to R dimnames standards (replace all non-numeric characters with "_")
is.letter <- function(x) grepl("[[:alpha:]]", x)
alphaNumeric<-function(s){
    first<-sapply(s,function(x){strsplit(x,"")[[1]][1]})
    new<-sapply(1:length(s),function(i){ ifelse(is.letter(first[i]),s[i],paste0("X",s[i]))})
    return(gsub("[^[:alnum:]]", "_", new))
}

signatureListFromFile <- function(fl, name=1, start=3){
    #make a list containing all the signatures from a text file
    # if(!grepl(".txt", fl)){error("please use a txt file")}
    con <- file(fl)
    open(con);
    signatures.list <- list();

    while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
        vec<- unlist(strsplit(line, split="\t"))
        signatures.list[[vec[name]]]<-vec[start:length(vec)]
    }
    close(con)
    rm(con,line,vec)
    return(signatures.list)
}



calculateSignatureExpression<-function(geneExpr,signatures){
    library(GSVA)
    # RPKM values are the input
    sigExpression<-gsva(as.matrix(geneExpr), signatures, min.sz=5, max.sz=500, rnaseq=TRUE, method="gsva",parallel.sz=4)
    # change between algorithms with method parameter: method=c("gsva", "ssgsea", "zscore", "plage")
    signExpr<-sigExpression$es.obs
    row.names(signExpr)<-alphaNumeric(row.names(signExpr))
    return(data.frame(signExpr))
}

pairwise.var.test<-function(values,groups){
    #message("no p-value adjustments are implemented!")
    groups<-factor(groups)
    p.vals<-matrix(nrow=length(levels(groups))-1,ncol=length(levels(groups))-1)
    colnames(p.vals)<-paste0("col",seq(1,length(levels(groups))-1,1))
    row.names(p.vals)<-paste0("row",seq(1,length(levels(groups))-1,1))
    for(i in 1:(length(levels(groups))-1)){
        for(j in i+1:(length(levels(groups)))){
            if(j>length(levels(groups))){break}
            # message(paste0(i,"\t",j))
            x<-values[which(groups==as.character(levels(groups)[i]))]
            y<-values[which(groups==as.character(levels(groups)[j]))]
            p.vals[j-1,i]<-var.test(x,y)$p.value
            colnames(p.vals)[i]<-as.character(levels(groups)[i])
            row.names(p.vals)[j-1]<-as.character(levels(groups)[j])
        }

    }
    return(p.vals)
}

mergeByRowNames<-function(df1,df2,...){
    df.out<-merge(df1,df2,by="row.names",...)
    row.names(df.out)<-df.out$Row.names
    df.out$Row.names<-NULL
    return(df.out)
}

mergeListByRowNames <- function(dflist,...){
    finaldf = dflist[[1]]
    for(i in 2:length(dflist)){
        finaldf = mergeByRowNames(finaldf, dflist[[i]],...)
    }
    return(finaldf)
}


stderr<-function(x){sd(x,na.rm=T)/sqrt(length(x))}


setTheme<-function(){
    require(ggplot2)
    theme_set(theme_bw())
    theme_update(strip.background = element_rect(colour = "black"),
                 strip.text = element_text(size=12),
                 axis.text = element_text(size=12),
                 axis.text.x = element_text(angle=90,vjust=0.8,hjust=0.5,size=12),
                 axis.ticks.length = unit(0.1, "cm"),
                 #axis.text.x = element_blank(),
                 legend.key = element_rect(fill = "white", color = "white"),
                 #panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_rect(colour = "black"),
                 panel.background = element_blank(),
                 axis.line.x = element_line(colour = "black"),
                 axis.line.y = element_line(colour = "black"))

}


colfunc <- function(n,skew=1){
    f1<-colorRampPalette(c("blue","white","red"), bias=1)
    clr <- f1(floor(n/2))
    sk=ceiling(n* ((1 - 1/skew)/2))
    c(rep(clr[1],sk), clr, rep(clr[length(clr)],sk) )
}

round.pval <- function(p,digits=2){
    round(p,floor(-log10(p))+digits)
}

last <- function(x) { return( x[length(x)] ) }

showDuplicates = function(df, column){
    df[df[,column] %in% df[which(duplicated(df[,column])),column], ]
}

splitVector = function(string_to_split="", split="_", element=1, ...){
    sapply(string_to_split, function(x){strsplit(x,split, ...)[[1]][element]})
}

splitVectorAndPaste = function (string_to_split = "", split = "_", elements = 1, collapse="_", ...)
{
    unname(sapply(string_to_split, function(x) {
        paste(strsplit(x, split, ...)[[1]][elements], collapse=collapse)
    }))
}


pval2stars = function(p){
    ifelse(p<0.001,"***",ifelse(p<0.01, "**",ifelse(p<0.05, "*","")))
}

transparent<-function(someColor, alpha=100){
    newColor<-col2rgb(someColor)
    apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                                blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

addStarsToCorrplot = function(full=clinical, feature_names1=signames, feature_names2=feature_names1, star.cex=1.5, star.offset=0.22){
    # calculate significance of correlations, for adding stars to the plot
    cor_grid = expand.grid(feature1=feature_names1, feature2=feature_names2, stringsAsFactors = F)
    location_1 = data.frame(row.names=feature_names1, coordinate=seq(1,length(feature_names1),1))
    location_2 = data.frame(row.names=feature_names2, coordinate=seq(length(feature_names2),1,-1))
    cor_grid$coord.x = location_1[cor_grid$feature1,"coordinate"]
    cor_grid$coord.y = location_2[cor_grid$feature2,"coordinate"]
    cor_grid = cor_grid[cor_grid$feature1 != cor_grid$feature2,] # in case
    cor_grid$p.val = apply(cor_grid,1,function(x){ cor.test(full[,x[1] ],full[, x[2]])$p.val })

    # add the stars
    for(i in 1:nrow(cor_grid)){
        text(x=cor_grid$coord.x[i],y=cor_grid$coord.y[i]+star.offset,labels=pval2stars(cor_grid$p.val[i]), adj=c(0.5,0.5), cex=star.cex)
    }
}


