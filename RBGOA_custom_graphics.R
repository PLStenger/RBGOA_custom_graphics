# Custom by Pierre-Louis STENGER - Pierrelouis.stenger@gmail.com - phd student - 06-02-2019 

# All code is from https://github.com/z0on/GO_MWU for Rank-based Gene Ontology Analysis with Adaptive Clustering
# Here I just custom the gomwuPlot
# In order to: i) don't show the genes ("18/45" for example); ii) reverse the graph in other side (cluster tree at the right) in order to arrange later many graphics, iii) combine the both.

# So there 3 customized function from the incredible work of Mikhail V Matz (https://github.com/z0on/GO_MWU):
# gomwuPlot_reverse_without_genes()
# gomwuPlot_reverse()
# gomwuPlot_without_genes()


#####################################################################################################################################################################
######################################## gomwuPlot_reverse_without_genes ############################################################################################
#####################################################################################################################################################################
gomwuPlot_reverse_without_genes=function(inFile,goAnnotations,goDivision,level1=0.1,level2=0.05,level3=0.01,absValue=-log(0.05,10),adjusted=TRUE,txtsize=1,font.family="sans",treeHeight=0.5,colors=NULL) {
  require(ape)
  
  input=inFile
  in.mwu=paste("MWU",goDivision,input,sep="_")
  i
  
  level1=0.1
  level2=0.05
  level3=0.01
  absValue=-log(0.05,10)
  adjusted=TRUE
  txtsize=1
  font.family="sans"
  treeHeight=0.5
  colors=NULL
  
  #input=inFile
  in.mwu=paste("MWU",goDivision,input,sep="_")
  in.dissim=paste("dissim",goDivision,goAnnotations,sep="_")
  
  cutoff=-log(level1,10)
  pv=read.table(in.mwu,header=T)
  row.names(pv)=pv$term
  in.raw=paste(goDivision,input,sep="_")
  rsq=read.table(in.raw,sep="\t",header=T)
  rsq$term=as.factor(rsq$term)
  
  if (adjusted==TRUE) { pvals=pv$p.adj } else { pvals=pv$pval }
  heat=data.frame(cbind("pval"=pvals)) 
  row.names(heat)=pv$term
  heat$pval=-log(heat$pval+1e-15,10)
  heat$direction=0
  heat$direction[pv$delta.rank>0]=1
  if (cutoff>0) { 
    goods=subset(heat,pval>=cutoff) 
  } else {
    goods.names=unique(rsq$term[abs(rsq$value)>=absValue])
    goods=heat[row.names(heat) %in% goods.names,]
  }
  
  if (is.null(colors) | length(colors)<4 ) {
    colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral")
    if (sum(goods$direction)==nrow(goods) | sum(goods$direction)==0) { 
      colors=c("black","black","grey50","grey50")
    }
  }
  goods.names=row.names(goods)
  
  # reading and subsetting dissimilarity matrix
  diss=read.table(in.dissim,sep="\t",header=T,check.names=F)
  row.names(diss)=names(diss)
  diss.goods=diss[goods.names,goods.names]
  
  # how many genes out of what we started with we account for with our best categories?
  good.len=c();good.genes=c()
  for (g in goods.names) {
    sel=rsq[rsq$term==g,]	
    pass=abs(sel$value)>=absValue
    sel=sel[pass,]
    good.genes=append(good.genes,as.character(sel$seq))
    good.len=append(good.len,nrow(sel))
  }
  ngenes=length(unique(good.genes))
  
  ################### HERE TO DELETE GENES NUMBERS
  
  #hist(rsq$value)
  totSum=length(unique(rsq$seq[abs(rsq$value)>=absValue]))
  # row.names(goods)=paste(good.len,"/",pv[pv$term %in% goods.names,]$nseqs," ",pv[pv$term %in% goods.names,]$name,sep="")
  row.names(goods)=paste(pv[pv$term %in% goods.names,]$name,sep="") # modifier ############################################################################################# 
  #row.names(heat)=paste(good.len,"/",pv$nseqs," ",pv$name,sep="")
  row.names(heat)=paste(pv$name,sep="") # modifier ####################################
  #  row.names(diss.goods)=paste(good.len,"/",pv[pv$term %in% goods.names,]$nseqs," ",pv[pv$term %in% goods.names,]$name,sep="")
  row.names(diss.goods)=paste(pv[pv$term %in% goods.names,]$name,sep="")
  
  
  
  # clustering terms better than cutoff
  GO.categories=as.dist(diss.goods)
  cl.goods=hclust(GO.categories,method="average")
  labs=cl.goods$labels[cl.goods$order] # saving the labels to order the plot
  goods=goods[labs,]
  labs=sub(" activity","",labs)
  
  old.par <- par( no.readonly = TRUE )
  
  
  ########## ORGANISATION PLOT
  
  plots=layout(matrix(c(1,2,3),1,3,byrow=T),c(1,3,treeHeight),T)
  
  
  ########## P VALUE
  #par(mar = c(3,1,1,0))
  par(mar = c(3,0,1,0))
  
  
  plot(c(1:top)~c(1:top),type="n",axes=F,xlab="",ylab="")
  text(left,top-step*2,paste("p < ",level3,sep=""),font=2,cex=1* txtsize,adj=c(0,0),family=font.family)
  text(left,top-step*3,paste("p < ",level2,sep=""),font=1,cex=0.8* txtsize,adj=c(0,0),family=font.family)
  text(left,top-step*4,paste("p < ",10^(-cutoff),sep=""),font=3,col="grey50",cex=0.8* txtsize,adj=c(0,0),family=font.family)
  
  ########## GO CAT
  step=100
  left=1
  top=step*(2+length(labs))
  
  #par(mar = c(0,0,0.3,0))
  par(mar = c(0,0,0.3,0))
  
  
  plot(c(1:top)~c(1:top),type="n",axes=F,xlab="",ylab="")
  
  ii=1
  goods$color=1
  goods$color[goods$direction==1 & goods$pval>cutoff]=colors[4]
  goods$color[goods$direction==0 & goods$pval>cutoff]=colors[3]
  goods$color[goods$direction==1 & goods$pval>(-log(level2,10))]=colors[2]
  goods$color[goods$direction==0 & goods$pval>(-log(level2,10))]=colors[1]
  goods$color[goods$direction==1 & goods$pval>(-log(level3,10))]=colors[2]
  goods$color[goods$direction==0 & goods$pval>(-log(level3,10))]=colors[1]
  for (i in length(labs):1) {
    ypos=top-step*ii
    ii=ii+1
    if (goods$pval[i]> -log(level3,10)) { 
      text(2800,ypos,labs[i],font=2,cex=1*txtsize,col=goods$color[i],adj=c(0,0),family=font.family, pos=2)
    } else {
      if (goods$pval[i]>-log(level2,10)) { 
        text(2800,ypos,labs[i],font=1,cex=0.8* txtsize,col=goods$color[i],adj=c(0,0),family=font.family, pos=2)
      } else {
        #			if (goods$pval[i]>cutoff) { 
        #				text(left,ypos,labs[i],font=3,cex=0.8* txtsize,col=goods$color[i],adj=c(0,0),family=font.family)
        #		} else { 
        text(2800,ypos,labs[i],font=3,cex=0.8* txtsize,col=goods$color[i],adj=c(0,0),family=font.family, pos = 2) 
        #}
      }}}
  
  
  ########## CLUSTER
  
  #  par(mar = c(2,2,0.85,0))
  par(mar = c(2,0,0.85,0))
  
  # plot(as.phylo(cl.goods),show.tip.label=FALSE,cex=0.0000001 # normal plot
  plot(as.phylo(cl.goods),show.tip.label=FALSE,cex=0.0000001, x.lim=c(0.5, 0))
  

  
  cat(paste("GO terms dispayed: ",length(goods.names)),"\n")
  cat(paste("\"Good genes\" accounted for:  ", ngenes," out of ",totSum, " ( ",round(100*ngenes/totSum,0), "% )","\n",sep=""))
  par(old.par)	
  goods$pval=10^(-1*goods$pval)
  return(goods)
}




#####################################################################################################################################################################
######################################## gomwuPlot_reverse ##########################################################################################################
#####################################################################################################################################################################
gomwuPlot_reverse=function(inFile,goAnnotations,goDivision,level1=0.1,level2=0.05,level3=0.01,absValue=-log(0.05,10),adjusted=TRUE,txtsize=1,font.family="sans",treeHeight=0.5,colors=NULL) {
  require(ape)
  
  input=inFile
  in.mwu=paste("MWU",goDivision,input,sep="_")
  i
  
  level1=0.1
  level2=0.05
  level3=0.01
  absValue=-log(0.05,10)
  adjusted=TRUE
  txtsize=1
  font.family="sans"
  treeHeight=0.5
  colors=NULL
  
  #input=inFile
  in.mwu=paste("MWU",goDivision,input,sep="_")
  in.dissim=paste("dissim",goDivision,goAnnotations,sep="_")
  
  cutoff=-log(level1,10)
  pv=read.table(in.mwu,header=T)
  row.names(pv)=pv$term
  in.raw=paste(goDivision,input,sep="_")
  rsq=read.table(in.raw,sep="\t",header=T)
  rsq$term=as.factor(rsq$term)
  
  if (adjusted==TRUE) { pvals=pv$p.adj } else { pvals=pv$pval }
  heat=data.frame(cbind("pval"=pvals)) 
  row.names(heat)=pv$term
  heat$pval=-log(heat$pval+1e-15,10)
  heat$direction=0
  heat$direction[pv$delta.rank>0]=1
  if (cutoff>0) { 
    goods=subset(heat,pval>=cutoff) 
  } else {
    goods.names=unique(rsq$term[abs(rsq$value)>=absValue])
    goods=heat[row.names(heat) %in% goods.names,]
  }
  
  if (is.null(colors) | length(colors)<4 ) {
    colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral")
    if (sum(goods$direction)==nrow(goods) | sum(goods$direction)==0) { 
      colors=c("black","black","grey50","grey50")
    }
  }
  goods.names=row.names(goods)
  
  # reading and subsetting dissimilarity matrix
  diss=read.table(in.dissim,sep="\t",header=T,check.names=F)
  row.names(diss)=names(diss)
  diss.goods=diss[goods.names,goods.names]
  
  # how many genes out of what we started with we account for with our best categories?
  good.len=c();good.genes=c()
  for (g in goods.names) {
    sel=rsq[rsq$term==g,]	
    pass=abs(sel$value)>=absValue
    sel=sel[pass,]
    good.genes=append(good.genes,as.character(sel$seq))
    good.len=append(good.len,nrow(sel))
  }
  ngenes=length(unique(good.genes))
  
  ################### HERE TO DELETE GENES NUMBERS
  
  #hist(rsq$value)
  totSum=length(unique(rsq$seq[abs(rsq$value)>=absValue]))
  row.names(goods)=paste(good.len,"/",pv[pv$term %in% goods.names,]$nseqs," ",pv[pv$term %in% goods.names,]$name,sep="")
  row.names(heat)=paste(good.len,"/",pv$nseqs," ",pv$name,sep="")
  row.names(diss.goods)=paste(good.len,"/",pv[pv$term %in% goods.names,]$nseqs," ",pv[pv$term %in% goods.names,]$name,sep="")
  
  
  # clustering terms better than cutoff
  GO.categories=as.dist(diss.goods)
  cl.goods=hclust(GO.categories,method="average")
  labs=cl.goods$labels[cl.goods$order] # saving the labels to order the plot
  goods=goods[labs,]
  labs=sub(" activity","",labs)
  
  old.par <- par( no.readonly = TRUE )
  
  
  ########## ORGANISATION PLOT
  
  plots=layout(matrix(c(1,2,3),1,3,byrow=T),c(1,3,treeHeight),T)
  
  
  ########## P VALUE
  #par(mar = c(3,1,1,0))
  par(mar = c(3,0,1,0))
  
  
  plot(c(1:top)~c(1:top),type="n",axes=F,xlab="",ylab="")
  text(left,top-step*2,paste("p < ",level3,sep=""),font=2,cex=1* txtsize,adj=c(0,0),family=font.family)
  text(left,top-step*3,paste("p < ",level2,sep=""),font=1,cex=0.8* txtsize,adj=c(0,0),family=font.family)
  text(left,top-step*4,paste("p < ",10^(-cutoff),sep=""),font=3,col="grey50",cex=0.8* txtsize,adj=c(0,0),family=font.family)
  
  ########## GO CAT
  step=100
  left=1
  top=step*(2+length(labs))
  
  #par(mar = c(0,0,0.3,0))
  par(mar = c(0,0,0.3,0))
  
  
  plot(c(1:top)~c(1:top),type="n",axes=F,xlab="",ylab="")
  
  ii=1
  goods$color=1
  goods$color[goods$direction==1 & goods$pval>cutoff]=colors[4]
  goods$color[goods$direction==0 & goods$pval>cutoff]=colors[3]
  goods$color[goods$direction==1 & goods$pval>(-log(level2,10))]=colors[2]
  goods$color[goods$direction==0 & goods$pval>(-log(level2,10))]=colors[1]
  goods$color[goods$direction==1 & goods$pval>(-log(level3,10))]=colors[2]
  goods$color[goods$direction==0 & goods$pval>(-log(level3,10))]=colors[1]
  for (i in length(labs):1) {
    ypos=top-step*ii
    ii=ii+1
    if (goods$pval[i]> -log(level3,10)) { 
      text(2800,ypos,labs[i],font=2,cex=1*txtsize,col=goods$color[i],adj=c(0,0),family=font.family, pos=2)
    } else {
      if (goods$pval[i]>-log(level2,10)) { 
        text(2800,ypos,labs[i],font=1,cex=0.8* txtsize,col=goods$color[i],adj=c(0,0),family=font.family, pos=2)
      } else {
        #			if (goods$pval[i]>cutoff) { 
        #				text(left,ypos,labs[i],font=3,cex=0.8* txtsize,col=goods$color[i],adj=c(0,0),family=font.family)
        #		} else { 
        text(2800,ypos,labs[i],font=3,cex=0.8* txtsize,col=goods$color[i],adj=c(0,0),family=font.family, pos = 2) 
        #}
      }}}
  
  
  ########## CLUSTER
  
  #  par(mar = c(2,2,0.85,0))
  par(mar = c(2,0,0.85,0))
  
  # plot(as.phylo(cl.goods),show.tip.label=FALSE,cex=0.0000001 # normal plot
  plot(as.phylo(cl.goods),show.tip.label=FALSE,cex=0.0000001, x.lim=c(0.5, 0))
  
  
  
  cat(paste("GO terms dispayed: ",length(goods.names)),"\n")
  cat(paste("\"Good genes\" accounted for:  ", ngenes," out of ",totSum, " ( ",round(100*ngenes/totSum,0), "% )","\n",sep=""))
  par(old.par)	
  goods$pval=10^(-1*goods$pval)
  return(goods)
}




#####################################################################################################################################################################
######################################## gomwuPlot_without_genes ############################################################################################
#####################################################################################################################################################################
gomwuPlot_without_genes=function(inFile,goAnnotations,goDivision,level1=0.1,level2=0.05,level3=0.01,absValue=-log(0.05,10),adjusted=TRUE,txtsize=1,font.family="sans",treeHeight=0.5,colors=NULL) {
  require(ape)
  
  input=inFile
  in.mwu=paste("MWU",goDivision,input,sep="_")
  i
  
  level1=0.1
  level2=0.05
  level3=0.01
  absValue=-log(0.05,10)
  adjusted=TRUE
  txtsize=1
  font.family="sans"
  treeHeight=0.5
  colors=NULL
  
  #input=inFile
  in.mwu=paste("MWU",goDivision,input,sep="_")
  in.dissim=paste("dissim",goDivision,goAnnotations,sep="_")
  
  cutoff=-log(level1,10)
  pv=read.table(in.mwu,header=T)
  row.names(pv)=pv$term
  in.raw=paste(goDivision,input,sep="_")
  rsq=read.table(in.raw,sep="\t",header=T)
  rsq$term=as.factor(rsq$term)
  
  if (adjusted==TRUE) { pvals=pv$p.adj } else { pvals=pv$pval }
  heat=data.frame(cbind("pval"=pvals)) 
  row.names(heat)=pv$term
  heat$pval=-log(heat$pval+1e-15,10)
  heat$direction=0
  heat$direction[pv$delta.rank>0]=1
  if (cutoff>0) { 
    goods=subset(heat,pval>=cutoff) 
  } else {
    goods.names=unique(rsq$term[abs(rsq$value)>=absValue])
    goods=heat[row.names(heat) %in% goods.names,]
  }
  
  if (is.null(colors) | length(colors)<4 ) {
    colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral")
    if (sum(goods$direction)==nrow(goods) | sum(goods$direction)==0) { 
      colors=c("black","black","grey50","grey50")
    }
  }
  goods.names=row.names(goods)
  
  # reading and subsetting dissimilarity matrix
  diss=read.table(in.dissim,sep="\t",header=T,check.names=F)
  row.names(diss)=names(diss)
  diss.goods=diss[goods.names,goods.names]
  
  # how many genes out of what we started with we account for with our best categories?
  good.len=c();good.genes=c()
  for (g in goods.names) {
    sel=rsq[rsq$term==g,]	
    pass=abs(sel$value)>=absValue
    sel=sel[pass,]
    good.genes=append(good.genes,as.character(sel$seq))
    good.len=append(good.len,nrow(sel))
  }
  ngenes=length(unique(good.genes))
  
  ################### HERE TO DELETE GENES NUMBERS
  
  #hist(rsq$value)
  totSum=length(unique(rsq$seq[abs(rsq$value)>=absValue]))
  # row.names(goods)=paste(good.len,"/",pv[pv$term %in% goods.names,]$nseqs," ",pv[pv$term %in% goods.names,]$name,sep="")
  row.names(goods)=paste(pv[pv$term %in% goods.names,]$name,sep="") # modifier ############################################################################################# 
  #row.names(heat)=paste(good.len,"/",pv$nseqs," ",pv$name,sep="")
  row.names(heat)=paste(pv$name,sep="") # modifier ####################################
  #  row.names(diss.goods)=paste(good.len,"/",pv[pv$term %in% goods.names,]$nseqs," ",pv[pv$term %in% goods.names,]$name,sep="")
  row.names(diss.goods)=paste(pv[pv$term %in% goods.names,]$name,sep="")
  
  # clustering terms better than cutoff
  GO.categories=as.dist(diss.goods)
  cl.goods=hclust(GO.categories,method="average")
  labs=cl.goods$labels[cl.goods$order] # saving the labels to order the plot
  goods=goods[labs,]
  labs=sub(" activity","",labs)
  
  old.par <- par( no.readonly = TRUE )
  
  plots=layout(matrix(c(1,2,3),1,3,byrow=T),c(treeHeight,3,1),TRUE)
  
  par(mar = c(2,2,0.85,0))
  plot(as.phylo(cl.goods),show.tip.label=FALSE,cex=0.0000001)
  step=100
  left=1
  top=step*(2+length(labs))
  
  par(mar = c(0,0,0.3,0))
  plot(c(1:top)~c(1:top),type="n",axes=F,xlab="",ylab="")
  ii=1
  goods$color=1
  goods$color[goods$direction==1 & goods$pval>cutoff]=colors[4]
  goods$color[goods$direction==0 & goods$pval>cutoff]=colors[3]
  goods$color[goods$direction==1 & goods$pval>(-log(level2,10))]=colors[2]
  goods$color[goods$direction==0 & goods$pval>(-log(level2,10))]=colors[1]
  goods$color[goods$direction==1 & goods$pval>(-log(level3,10))]=colors[2]
  goods$color[goods$direction==0 & goods$pval>(-log(level3,10))]=colors[1]
  for (i in length(labs):1) {
    ypos=top-step*ii
    ii=ii+1
    if (goods$pval[i]> -log(level3,10)) { 
      text(left,ypos,labs[i],font=2,cex=1*txtsize,col=goods$color[i],adj=c(0,0),family=font.family) 
    } else {
      if (goods$pval[i]>-log(level2,10)) { 
        text(left,ypos,labs[i],font=1,cex=0.8* txtsize,col=goods$color[i],adj=c(0,0),family=font.family)
      } else {
        #			if (goods$pval[i]>cutoff) { 
        #				text(left,ypos,labs[i],font=3,cex=0.8* txtsize,col=goods$color[i],adj=c(0,0),family=font.family)
        #		} else { 
        text(left,ypos,labs[i],font=3,cex=0.8* txtsize,col=goods$color[i],adj=c(0,0),family=font.family) 
        #}
      }
    }
  }
  
  par(mar = c(3,1,1,0))
  
  plot(c(1:top)~c(1:top),type="n",axes=F,xlab="",ylab="")
  text(left,top-step*2,paste("p < ",level3,sep=""),font=2,cex=1* txtsize,adj=c(0,0),family=font.family)
  text(left,top-step*3,paste("p < ",level2,sep=""),font=1,cex=0.8* txtsize,adj=c(0,0),family=font.family)
  text(left,top-step*4,paste("p < ",10^(-cutoff),sep=""),font=3,col="grey50",cex=0.8* txtsize,adj=c(0,0),family=font.family)
  
  cat(paste("GO terms dispayed: ",length(goods.names)),"\n")
  cat(paste("\"Good genes\" accounted for:  ", ngenes," out of ",totSum, " ( ",round(100*ngenes/totSum,0), "% )","\n",sep=""))
  par(old.par)	
  goods$pval=10^(-1*goods$pval)
  return(goods)
}

