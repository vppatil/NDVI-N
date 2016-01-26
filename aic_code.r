
aicc=function(aics,k,n)
{
	aicc= aics+ ((2*k*(k+1))/(n-k-1))
	return(aicc)
}

daic=function(aics){di=aics-min(aics); return(di)}

n.obs=function(x) sapply(x,nobs)

aic.weight=function(daics)
{
	weights=exp(-.5*daics)/sum(exp(-.5*daics))
	return(weights)
}

aic.table=function(aic.table,aics,k=NULL,n=NULL,aicc=FALSE)
{
	if(aicc)
	{
		aic.table$aicc=aicc(aic.table[[aics]],aic.table[[k]],aic.table[[n]])
		aic.table$daic=daic(aic.table$aicc)
		aic.table$use.aicc=1
	}
	else
	{
		aic.table$daic=daic(aic.table[[aics]])
		aic.table$w=aic.weight(aic.table$daic)
		aic.table$use.aicc=0
	}

	aic.table$w=aic.weight(aic.table$daic)
	##aic.table=aic.table[rev(order(aic.table$w)),] #don't reorder so models stay in consistent order.
	return(aic.table)
}

model.avg.prediction=function(model.list,aic.weights,dframe)
{
  weighteds=matrix(0,length(aic.weights),dim(dframe)[1])
  for(i in 1:dim(weighteds)[1])
  {
    weighteds[i,]=predict(model.list[[i]],dframe)*aic.weights[i]
  }
  weighteds=apply(weighteds,2,sum)
  return(weighteds)
}
