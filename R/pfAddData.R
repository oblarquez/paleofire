pfAddData=function(files,metadata,type="NULL",Int=TRUE,first=NULL,last=NULL,yrInterp=NULL){

## Data part
  
## Read data
for (i in 1:length(files))
  assign(paste("data",i,sep=""),read.csv(paste(files[i])))

## Calculates Charcoal Accumulation rates
if(type=="CharAnalysis"){
  for (i in 1:length(files)){
    temp=pretreatment(get(paste("data",i,sep=""))[,1:5],
                      get(paste("data",i,sep=""))[,6],
                      Int=Int,first=first,last=last,yrInterp=yrInterp)
    assign(paste("dataI",i,sep=""),na.omit(as.data.frame(cbind(10000+i,temp$cmI, temp$ybpI, temp$accI))))
  }
  all_data=do.call(rbind,lapply(1:length(files),function(i) rbind(get(paste("dataI",i,sep="")))))
}

## If data already specified as Influx (i.e. 3 columns csv's with Depth, Age, Influx)
if(type=="NULL")
  all_data=do.call(rbind,lapply(1:length(files),function(i) rbind(cbind(1000+i,get(paste("data",i,sep=""))))))

## Set colnames
colnames(all_data)=c("ID_SITE","DEPTH","EST_AGE","QUANTITY")

## Metadata part
options(warn=-1)
assign("meta",read.csv(paste(metadata),header=T))
meta[,1]=as.character(meta[,1])
class(meta[,1])="character"
meta=cbind(unique(all_data[,1]),meta)
colnames(meta)[1]="ID_SITE"

## Output
output=structure(list(data=all_data,
                      metadata=meta))
return(output)
}




