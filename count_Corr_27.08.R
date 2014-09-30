library("rtracklayer")
library("GenometriCorr")

library('doParallel')
#cl <- makeCluster(1)
registerDoParallel(cores=18)

###
#___how_to_start:
#R CMD BATCH --no-save --no-restore '--args RData_file_name_with_result_matrix_with_corr_stat name_of_wich_cluster_tracks_to_conside name_of_matrix_to_wich_correl_res_write' count_Corr_27.08.R
###

load.all.data<-function(input.file){
  
  con <- file(input.file, 'r') 
  modules_names <- readLines(con)
  close(con)
  
  modules.tracks<- c()
  cluster.active.tracks<- c()
  cluster.all.tracks<- c()
  cluster.brain.tracks<- c()
  
  
  
  for (i in 1:length(modules_names)){
    module_name <- as.character(modules_names[i])
    #data <- import(name)
    current_module.track <- as(import(paste("C",module_name,".tss.nr.bed",sep = "")), "RangedData")
    
    current_cluster.active.track <- as(import(paste("C",module_name,".clust.ge_2.active.bed",sep = "")) , "RangedData")   
    current_cluster.all.track <- as(import(paste("C",module_name,".clust.ge_2.all.bed",sep = "")), "RangedData")
    current_cluster.brain.track <- as(import(paste("C",module_name,".clust.ge_2.brain.bed",sep = "")) , "RangedData")
    
    modules.tracks<- c(modules.tracks, current_module.track)
    cluster.active.tracks<- c(cluster.active.tracks, current_cluster.active.track)
    cluster.all.tracks<- c(cluster.all.tracks, current_cluster.all.track)
    cluster.brain.tracks<-c(cluster.brain.tracks, current_cluster.brain.track)
  } 
  
  names(modules.tracks)<-modules_names
  names(cluster.active.tracks)<-paste(modules_names,".cluster.active",sep = "")
  names(cluster.all.tracks)<-paste(modules_names,".cluster.all",sep = "")
  names(cluster.brain.tracks)<-paste(modules_names,".cluster.brain",sep = "")
  
  return(list(modules.tracks,cluster.active.tracks,cluster.all.tracks,cluster.brain.tracks))
  
}

"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}

countCorr<-function(tracks.query, tracks.reference,result_array){
  #add essential information
  human.chrom.length <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 59373566,155270560)
  names(human.chrom.length) <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrY", "chrX")
  pn.area <- 100
  pn.dist <- 100
  pn.jacc <- 100
  
  
  for (name.query in names(tracks.query)){
    data.query<- tracks.query[name.query]
    #count<- 0
    #count_prev<-0
    #if(count-count_prev>5){
    # save(result_array, file = paste(result_array,'_v',count,sep = "_"))
    # count_prev<- count
    #}
    
    
    x<-c()
    x<-foreach (j=1:length(tracks.reference)) %dopar% {
      #for (j in 1:length(tracks.reference)){
      name.reference <- names(tracks.reference)[j]
      data.reference <- tracks.reference[name.reference]
      
      
      correlation_result  <- suppressWarnings(GenometriCorrelation(data.query[[1]],data.reference[[1]], chromosomes.length = human.chrom.length, chromosomes.to.proceed = c("chr1", "chr2", "chr3","chr4", "chr5","chr6", "chr7","chr8", "chr9","chr10", "chr11", "chr12", "chr13","chr14", "chr15","chr16", "chr17","chr18", "chr19","chr20", "chr21","chr22"), ecdf.area.permut.number = pn.area, mean.distance.permut.number = pn.dist,  jaccard.measure.permut.number = pn.jacc, keep.distributions = FALSE, showProgressBar = FALSE))
      #correlation_result
      #x<-c(x,list(correlation_result))
      
      #correlation_result  <- GenometriCorrelation(data.first[[1]],data.second[[1]], chromosomes.length = human.chrom.length, chromosomes.to.proceed = c("chr1"), ecdf.area.permut.number = pn.area, mean.distance.permut.number = pn.dist,  jaccard.measure.permut.number = pn.jacc, keep.distributions = FALSE, showProgressBar = FALSE)
      #correlation_result  <- GenometriCorrelation(data.first[[name.first]],data.second[[name.second]], chromosomes.length = human.chrom.length, chromosomes.to.proceed = c("chr1"), ecdf.area.permut.number = pn.area, mean.distance.permut.number = pn.dist,  jaccard.measure.permut.number = pn.jacc, keep.distributions = TRUE, showProgressBar = FALSE)        
    }
    names(x)<-names(tracks.reference)
    for (name in names(x)){
      result_array[name.query,name]<-x[name]
    }
  }
  return(list(result_array))
  #return(array.mtt)
  
}

create.result.array<-function(tracks_1, tracks_2){
  result_array <-array(list(NULL), dim=c(length(tracks_1),length(tracks_2)), dimnames= list(names(tracks_1),names(tracks_2)))
  return(result_array)
}

change_modules_starts<-function(modules.tracks, upstream){
  new.starts<-lapply(lapply(modules.tracks, start) ,   function(x) x-upstream) 
  for (module.name in names(modules.tracks)){
    #print(module.name)
    #print(names(modules.tracks[module.name]))
    start(modules.tracks[module.name][[1]]) <-new.starts[module.name][[1]] 
  }
  return(modules.tracks)
}

change_modules_ends<-function(modules.tracks, downstream){
  new.ends<-lapply(lapply(modules.tracks, end) ,   function(x) x+downstream) 
  for (module.name in names(modules.tracks)){
    #print(module.name)
    #print(names(modules.tracks[module.name]))
    end(modules.tracks[module.name][[1]]) <-new.ends[module.name][[1]] 
  }
  return(modules.tracks)
}

load.modules<-function(input.modules.list){
  con <- file(input.modules.list, 'r') 
  modules_names <- readLines(con)
  close(con)
  modules.tracks<- c()
  
  for (i in 1:length(modules_names)){
    current_module.track <- as(import(paste(modules_names[i],".bed",sep = "")), "RangedData") 
    modules.tracks<- c(modules.tracks, current_module.track)
  }  
  names(modules.tracks)<-modules_names
  return(modules.tracks)
}

make.clusters.coords<-function(modules.tracks){
  clusters.all<-c()
  clusters.active<-c()
  clusters.brain<-c()
  for (module.name in names(modules.tracks)){   
    #module.name<-names(modules.tracks)[1]
    current.module<-modules.tracks[module.name]
    starts<-lapply(current.module, start)[[1]]   
    names_for<-space(current.module[[1]])
    #starts<-c(10, 20, 30)
    left_coords<-sapply(starts, function(x) sample(x:(x-30000), 1))
    rigth_coords<-sapply(left_coords, function(x) sample(x:(x+5000), 1))
    current_cluster_track_all<-RangedData(ranges = IRanges(start=left_coords, end=rigth_coords), space=names_for)
    
    left_coords<-sapply(starts, function(x) sample(x:(x-20000), 1))
    rigth_coords<-sapply(left_coords, function(x) sample(x:(x+5000), 1))
    current_cluster_track_active<-RangedData(ranges = IRanges(start=left_coords, end=rigth_coords), space=names_for)
    
    left_coords<-sapply(starts, function(x) sample(x:(x-10000), 1))
    rigth_coords<-sapply(left_coords, function(x) sample(x:(x+5000), 1))
    current_cluster_track_brain<-RangedData(ranges = IRanges(start=left_coords, end=rigth_coords), space=names_for)
    
    clusters.brain<-c(clusters.brain, current_cluster_track_brain)
    clusters.active<-c(clusters.active, current_cluster_track_active)
    clusters.all<-c(clusters.all, current_cluster_track_all)
  }
  
  names(clusters.brain)<-names(modules.tracks)
  names(clusters.active)<-names(modules.tracks)
  names(clusters.all)<-names(modules.tracks)
  
  return(list(clusters.brain, clusters.active, clusters.all)) 
}

#########____MAIN_for_simple_count_corr____#################

#tracks <- structure(NA,class="result")
#tracks[modules.tracks,cluster.active.tracks,cluster.all.tracks,cluster.brain.tracks]<- load.all.data('modules_list.txt')

#modules.tracks<- change_modules_starts(modules.tracks, 10000)
#modules.tracks<- change_modules_ends(modules.tracks, 10000)

#list <- structure(NA,class="result")

##create arrays for results
#correlation_result.module_cluster_all<-create.result.array(modules.tracks, cluster.all.tracks)
#correlation_result.module_cluster_active<-create.result.array(modules.tracks, cluster.active.tracks)
#correlation_result.module_cluster_brain<-create.result.array(modules.tracks, cluster.brain.tracks)

##count correlations and write results to RData files
###list[correlation_result.module_cluster_all]<- countCorr(modules.tracks, cluster.all.tracks,correlation_result.module_cluster_all)
###save(correlation_result.module_cluster_all, file = "ttt_all_1.RData")

###list[correlation_result.module_cluster_active]<- countCorr(modules.tracks, cluster.active.tracks,correlation_result.module_cluster_active)
###save(correlation_result.module_cluster_active, file = "ttt_active_2.RData")

#list[correlation_result.module_cluster_brain]<- countCorr(modules.tracks, cluster.brain.tracks,correlation_result.module_cluster_brain)
#save(correlation_result.module_cluster_brain, file = "new_cl_1_3_10000_10000.RData")

##___MAIN_for_generation_clusters_coord_and_check_count____####

args <- commandArgs(trailingOnly = TRUE)

#args<-list(modules_list.txt',test.RData','brain', '-u', '10000', '-d', '10000', '-g') 

if ("-g" %in% args){
  modules.tracks<-load.modules(args[1])
  tracks <- structure(NA,class="result")
  tracks[cluster.brain.tracks, cluster.active.tracks, cluster.all.tracks]<-make.clusters.coords(modules.tracks)
  
}
else{
  tracks <- structure(NA,class="result")
  tracks[modules.tracks,cluster.active.tracks,cluster.all.tracks,cluster.brain.tracks]<- load.all.data('modules_list.txt')  
}

if ("-u" %in% args){
  upstream_value_index<-match('-u', args)+1
  modules.tracks<- change_modules_starts(modules.tracks,args[upstream_value_index] )
  
}
if ("-d" %in% args){
  downstream_value_index<-match('-d', args)+1
  modules.tracks<- change_modules_ends(modules.tracks, args[downstream_value_index])
}

#args[1]='modules_list.txt'

#modules.tracks<- change_modules_starts(modules.tracks, 10000)
#modules.tracks<- change_modules_ends(modules.tracks, 10000)

list <- structure(NA,class="result")

correlation_result.module_cluster_all<-create.result.array(modules.tracks, cluster.all.tracks)
correlation_result.module_cluster_active<-create.result.array(modules.tracks, cluster.active.tracks)
correlation_result.module_cluster_brain<-create.result.array(modules.tracks, cluster.brain.tracks)




#args[3]='brain'
#paste('cluster.',args[3],'.tracks', sep="")
#paste('correlation_result.module_cluster_',args[3], sep="")
list[correlation_result.module_cluster_brain]<- countCorr(modules.tracks, get(paste('cluster.',args[3],'.tracks', sep="")),get(paste('correlation_result.module_cluster_',args[3], sep="")))
save(correlation_result.module_cluster_brain, file=args[2])
