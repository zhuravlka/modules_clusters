library('RColorBrewer')
library('gnuplot')
library('gplots')


count.correlation.count.array <-function(result.data.frame, correlation.in.interest){
  #result.data.frame<-correlation_result.module_cluster_brain
  correlation.in.interest<-correlation.in.interest
  #create result.correlation.count.array
  #result.correlation.count.array <-array(list(NULL), dim=c(nrow(result.data.frame),ncol(result.data.frame)), dimnames= list(rownames(result.data.frame),colnames(result.data.frame)))
  
  result.correlation.count.matrix <-matrix(data= NA, nrow=nrow(result.data.frame), ncol=ncol(result.data.frame), byrow=FALSE,dimnames= list(rownames(result.data.frame),colnames(result.data.frame)))
  
  
  #substrate correlation counts and fill the array with them
  for (module.name in rownames(result.data.frame)) {
    for (cluster.name in colnames(result.data.frame)) {
      
      #correlation.in.interest<-'relative.distances.ecdf.area.correlation'
      
      #correlation.in.interest<-'jaccard.measure'
      
      #correlation.in.interest<-'scaled.absolute.min.distance.sum.p.value'
      
      #correlation.in.interest<-'projection.test.p.value'
      #true.false.list<- list("scaled.absolute.min.distance.sum.p.value", "projection.test.p.value")
      
      #if(correlation.in.interest %in% true.false.list){
      #if (correlation.in.interest == "scaled.absolute.min.distance.sum.p.value" ){
      #result.data.frame[unname(module.name), unname(cluster.name)][[1]]$awhole[correlation.in.interest][[1]]
      # }
      #if (result.data.frame[unname(module.name), unname(cluster.name)][[1]]$awhole$scaled.absolute.min.distance.sum.lower.tail)    
      
      #}
      
      
      
      if("awhole" %in% names(result.data.frame[unname(module.name), unname(cluster.name)][[1]])){
        #print('it is true')
        current.correlation.count <- result.data.frame[unname(module.name), unname(cluster.name)][[1]]$awhole[correlation.in.interest][[1]]
      
        if (is.na(as.numeric(current.correlation.count))){
          result.correlation.count.matrix[unname(module.name), unname(cluster.name)]<-0
        }
        else{
          result.correlation.count.matrix[unname(module.name), unname(cluster.name)]<-as.numeric(current.correlation.count)
        }
        
      }
      else {
        #print('false')
        #result.correlation.count.matrix[unname(module.name), unname(cluster.name)]<-c(0)
        result.correlation.count.matrix[module.name, cluster.name]<-c(0)
      }
      
    }
  }
  return(result.correlation.count.matrix)
}



heatmap.and.save<-function(matrix.brain,matrix.all, matrix.active, correlation.in.interest){
  name<-paste(correlation.in.interest,"matrix.brain", sep = "_")
  png(paste(name, ".png", sep=""), pointsize = 10,width = 600, height = 600)
  col.na.vector<-c()
  for (element in strsplit(colnames(matrix.brain), "\\.")){
    col.na.vector<-c(col.na.vector, paste(element[1], element[2], sep="."))
  }
  heatmap.2(matrix.brain,trace='none', col=bluered ,dendrogram="none", Colv=NA, Rowv=NA, scale="none", main=paste("brain",correlation.in.interest,  sep="_"), labRow=rownames(matrix.brain), labCol=col.na.vector)
  dev.off()
  
  
  name<-paste(correlation.in.interest,"matrix.all", sep = "_")
  png(paste(name, ".png", sep=""),pointsize = 10,width = 600, height = 600)
  col.na.vector<-c()
  for (element in strsplit(colnames(matrix.all), "\\.")){
    col.na.vector<-c(col.na.vector, paste(element[1], element[2], sep="."))
  }
  heatmap.2(matrix.all,col=bluered,trace='none',dendrogram="none",Colv=NA, Rowv=NA,   scale="none", main=paste("all",correlation.in.interest,  sep="_"), labRow=rownames(matrix.all), labCol=col.na.vector)
  dev.off()
  
  
  name<-paste(correlation.in.interest,"matrix.active", sep = "_")
  png(paste(name, ".png", sep=""),pointsize = 10,width = 600, height = 600)
  col.na.vector<-c()
  for (element in strsplit(colnames(matrix.active), "\\.")){
    col.na.vector<-c(col.na.vector, paste(element[1], element[2], sep="."))
  }
  heatmap.2(matrix.active,col=bluered,trace='none',dendrogram="none",Colv=NA, Rowv=NA, scale="none", main=paste("active",correlation.in.interest,  sep="_"), labRow=rownames(matrix.active), labCol=col.na.vector)
  dev.off()
  
  
}


#remake.matrix(matrix.brain,matrix.brain.lowertail)

remake.matrix <- function(pvalues, lowertails){
  
  for (current.row in rownames(pvalues)){
    for(current.column in colnames(pvalues)){
      #print(lowertails[current.row, current.column])
    
      if(lowertails[current.row, current.column] == 0){
        pvalues[current.row, current.column]<- 1-pvalues[current.row, current.column]
      }
      else {
        pvalues[current.row, current.column]<-pvalues[current.row, current.column]
      }
    }
  }
  return(pvalues)

}


#heatmap.and.save.pvalue2(new.matrix.brain,new.matrix.all, new.matrix.active, correlation.in.interest)

heatmap.and.save.pvalue2<-function(new.matrix.brain,new.matrix.all, new.matrix.active,  correlation.in.interest){

#heatmap.and.save.pvalue2<-function(matrix.brain,matrix.all, matrix.active,matrix.brain.lowertail,matrix.all.lowertail, matrix.active.lowertail, correlation.in.interest, correlation.lowertails){
  
  name<-paste(correlation.in.interest,"matrix.brain", sep = "_")
  png(paste(name, ".png", sep=""), pointsize = 10,width = 600, height = 600)
  col.na.vector<-c()
  for (element in strsplit(colnames(new.matrix.brain), "\\.")){
    col.na.vector<-c(col.na.vector, paste(element[1], element[2], sep="."))
  }
  

  #breaks = seq(0,max(matrix.brain),by=0.01)
  #breaks2 <- seq(0.05,1,length.out=5)
  #breaks1 <- seq(0, 0.05,length.out=5)
  #gradient1 <- colorpanel( 4, "white", "blue")
  #gradient2 <- colorpanel( 5, "white","red"  )
  #hm.colors <- c(rev(gradient1),gradient2[-1])
  #breaks <-c(breaks1, breaks2[-1])
  
  breaks <- c(seq(0, 0.05, length.out=5),seq(0.05, 1, length.out=6)[-1])
  hm.colors <- colorpanel(n=9,low="blue",mid="white",high="red")
  
  heatmap.2(new.matrix.brain,trace='none',breaks=breaks,col=hm.colors,dendrogram="none", Colv=NA, Rowv=NA, scale="none", main=paste("brain",correlation.in.interest,  sep="_"), labRow=rownames(new.matrix.brain), labCol=col.na.vector)
  dev.off()
  
  
  name<-paste(correlation.in.interest,"matrix.all", sep = "_")
  png(paste(name, ".png", sep=""),pointsize = 10,width = 600, height = 600)
  col.na.vector<-c()
  for (element in strsplit(colnames(new.matrix.all), "\\.")){
    col.na.vector<-c(col.na.vector, paste(element[1], element[2], sep="."))
  }
  
  #breaks = seq(0,max(matrix.all),by=0.01)
  breaks <- c(seq(0, 0.05, length.out=5),seq(0.05, 1, length.out=6)[-1])
  hm.colors <- colorpanel(n=9,low="blue",mid="white",high="red")
  
  #heatmap.2(matrix.brain,trace='none',breaks=breaks,col=hm.colors,dendrogram="none", Colv=NA, Rowv=NA, scale="none", main=paste("brain",correlation.in.interest,  sep="_"), labRow=rownames(matrix.brain), labCol=col.na.vector)
  
  heatmap.2(new.matrix.all,trace='none',breaks=breaks,col=hm.colors,dendrogram="none",Colv=NA, Rowv=NA,   scale="none", main=paste("all",correlation.in.interest,  sep="_"), labRow=rownames(new.matrix.all), labCol=col.na.vector)
  dev.off()
  
  
  name<-paste(correlation.in.interest,"matrix.active", sep = "_")
  png(paste(name, ".png", sep=""),pointsize = 10,width = 600, height = 600)
  col.na.vector<-c()
  for (element in strsplit(colnames(new.matrix.active), "\\.")){
    col.na.vector<-c(col.na.vector, paste(element[1], element[2], sep="."))
  }
  
  #breaks = seq(0,max(matrix.active),by=0.01)
  
  breaks <- c(seq(0, 0.05, length.out=5),seq(0.05, 1, length.out=6)[-1])
  hm.colors <- colorpanel(n=9,low="blue",mid="white",high="red")
  
  heatmap.2(new.matrix.active,trace='none',breaks=breaks,col=hm.colors,dendrogram="none",Colv=NA, Rowv=NA, scale="none", main=paste("active",correlation.in.interest,  sep="_"), labRow=rownames(new.matrix.active), labCol=col.na.vector)
  dev.off()
   
}

#########_________MAIN_____#############
#setwd("F:/work/Favorov_data/favorov_data")

#load('ttt_all_1.RData')
#load('ttt_active_2.RData')
#load('ttt_brain_3.RData')

#load('ttt_all_1_10000.RData')
##load('ttt_active_2_10000.RData')
#load('ttt_brain_3_10000.RData')

#load('new_cl_1_1_10000_10000.RData')
#load('new_cl_1_2_10000_10000.RData')
#load('new_cl_1_3_10000_10000.RData')

load('gener_1.RData')
load('gener_2.RData')
load('gener_3.RData')

#load('new_cl_1_1.RData')
#load('new_cl_1_2.RData')
#load('new_cl_1_3.RData')

correlation.in.interest<-'jaccard.measure'
matrix.brain<-count.correlation.count.array(correlation_result.module_cluster_brain, correlation.in.interest)
matrix.all<-count.correlation.count.array(correlation_result.module_cluster_all, correlation.in.interest)
matrix.active<-count.correlation.count.array(correlation_result.module_cluster_active,correlation.in.interest)
heatmap.and.save(matrix.brain,matrix.all, matrix.active, correlation.in.interest)

correlation.in.interest<-'relative.distances.ecdf.area.correlation'
matrix.brain<-count.correlation.count.array(correlation_result.module_cluster_brain, correlation.in.interest)
matrix.all<-count.correlation.count.array(correlation_result.module_cluster_all, correlation.in.interest)
matrix.active<-count.correlation.count.array(correlation_result.module_cluster_active,correlation.in.interest)
heatmap.and.save(matrix.brain,matrix.all, matrix.active, correlation.in.interest)

############################___pvalue_statistics

correlation.in.interest<-'jaccard.measure.p.value'
correlation.lowertails<-'jaccard.measure.lower.tail'



matrix.brain<-count.correlation.count.array(correlation_result.module_cluster_brain, correlation.in.interest)
matrix.all<-count.correlation.count.array(correlation_result.module_cluster_all, correlation.in.interest)
matrix.active<-count.correlation.count.array(correlation_result.module_cluster_active,correlation.in.interest)

matrix.brain.lowertail<- count.correlation.count.array(correlation_result.module_cluster_brain, correlation.lowertails)
matrix.all.lowertail<-count.correlation.count.array(correlation_result.module_cluster_all, correlation.lowertails)
matrix.active.lowertail<-count.correlation.count.array(correlation_result.module_cluster_active, correlation.lowertails)

new.matrix.brain<-remake.matrix(matrix.brain,matrix.brain.lowertail)
new.matrix.all<-remake.matrix(matrix.all,matrix.all.lowertail)
new.matrix.active<-remake.matrix(matrix.active,matrix.active.lowertail)

heatmap.and.save.pvalue2(new.matrix.brain,new.matrix.all, new.matrix.active, correlation.in.interest)



correlation.in.interest<-'scaled.absolute.min.distance.sum.p.value'
correlation.lowertails<-'scaled.absolute.min.distance.sum.lower.tail'
matrix.brain<-count.correlation.count.array(correlation_result.module_cluster_brain, correlation.in.interest)
matrix.all<-count.correlation.count.array(correlation_result.module_cluster_all, correlation.in.interest)
matrix.active<-count.correlation.count.array(correlation_result.module_cluster_active,correlation.in.interest)

matrix.brain.lowertail<- count.correlation.count.array(correlation_result.module_cluster_brain, correlation.lowertails)
matrix.all.lowertail<-count.correlation.count.array(correlation_result.module_cluster_all, correlation.lowertails)
matrix.active.lowertail<-count.correlation.count.array(correlation_result.module_cluster_active, correlation.lowertails)

new.matrix.brain<-remake.matrix(matrix.brain,matrix.brain.lowertail)
new.matrix.all<-remake.matrix(matrix.all,matrix.all.lowertail)
new.matrix.active<-remake.matrix(matrix.active,matrix.active.lowertail)

heatmap.and.save.pvalue2(new.matrix.brain,new.matrix.all, new.matrix.active, correlation.in.interest)


correlation.in.interest<-'projection.test.p.value'
correlation.lowertails<-'projection.test.lower.tail'
matrix.brain<-count.correlation.count.array(correlation_result.module_cluster_brain, correlation.in.interest)
matrix.all<-count.correlation.count.array(correlation_result.module_cluster_all, correlation.in.interest)
matrix.active<-count.correlation.count.array(correlation_result.module_cluster_active,correlation.in.interest)

matrix.brain.lowertail<- count.correlation.count.array(correlation_result.module_cluster_brain, correlation.lowertails)
matrix.all.lowertail<-count.correlation.count.array(correlation_result.module_cluster_all, correlation.lowertails)
matrix.active.lowertail<-count.correlation.count.array(correlation_result.module_cluster_active, correlation.lowertails)

new.matrix.brain<-remake.matrix(matrix.brain,matrix.brain.lowertail)
new.matrix.all<-remake.matrix(matrix.all,matrix.all.lowertail)
new.matrix.active<-remake.matrix(matrix.active,matrix.active.lowertail)

heatmap.and.save.pvalue2(new.matrix.brain,new.matrix.all, new.matrix.active, correlation.in.interest)

################_pvalue_unsigned

correlation.in.interest<-'relative.distances.ecdf.deviation.area.p.value'
matrix.brain<-count.correlation.count.array(correlation_result.module_cluster_brain, correlation.in.interest)
matrix.all<-count.correlation.count.array(correlation_result.module_cluster_all, correlation.in.interest)
matrix.active<-count.correlation.count.array(correlation_result.module_cluster_active,correlation.in.interest)
heatmap.and.save(matrix.brain,matrix.all, matrix.active, correlation.in.interest)

correlation.in.interest<-'relative.distances.ks.p.value'
matrix.brain<-count.correlation.count.array(correlation_result.module_cluster_brain, correlation.in.interest)
matrix.all<-count.correlation.count.array(correlation_result.module_cluster_all, correlation.in.interest)
matrix.active<-count.correlation.count.array(correlation_result.module_cluster_active,correlation.in.interest)
heatmap.and.save(matrix.brain,matrix.all, matrix.active, correlation.in.interest)







