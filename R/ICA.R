ICA <-
function(cost,nvar,ncountries=80,nimp=10,maxiter=100,lb=-10,
                    ub=10,beta=2,P_revolve=0.3,zeta=0.02)
{
if(is.null(cost))
  {stop()}
if(is.null(nvar))
  {stop()}

ptm <- proc.time()
lb <- rep(lb,nvar)
ub <- rep(ub,nvar)
#---------------------- initial population algorithm

var <- c(rep(0,nvar))
fit=0
totalfit=0
y <- list(var=var,fit=fit)

colony <- list()
for(i in 1:ncountries)
{colony[[i]] <- y}

f <- c(rep(0,ncountries))
for(i in 1:ncountries)
{
  for(j in 1:nvar)
    {
     colony[[i]]$var[[j]] <- lb[j]+(runif(1)*(ub[j]-lb[j]))
     colony[[i]]$fit <- cost(colony[[i]]$var)
    }
}
for(i in 1:ncountries)
 { f[i] <- colony[[i]]$fit}
f1 <- sort(f)
k <- c(rep(0,ncountries))
for(i in 1:ncountries)
{
  for(j in 1:ncountries)
   {
      if(f1[i]==f[j])
      {
        k[i]=j
       break
      }
   }
  f[k]=0
}
colony2 <- colony
for(i in 1:ncountries)
{
  colony[[i]] <- colony2[[k[i]]]
}
#--------
d1 <- round((ncountries-nimp)/nimp)
d2 <- (ncountries-nimp)-(d1*(nimp-1))

imp <- list()
for(i in 1:(nimp-1))
imp[[i]] <- list(var=var,fit=fit,colonyi=list(),totalfit=totalfit)
imp[[nimp]] <- list(var=var,fit=fit,colonyi=list(),totalfit=totalfit)

for( i in 1:nimp)
{
 imp[[i]]$var <- colony[[i]]$var
 imp[[i]]$fit <- colony[[i]]$fit
}
#---------
colony3 <- list()

for(i in 1:(ncountries-nimp))
colony3[[i]] <- colony[[i+nimp]]

#---------
r <- 1:(ncountries-nimp)
r <- sample(r)
w=0
for(i in 1:(nimp-1))
{
  for(j in 1:d1)
   {
    w <- w+1
    imp[[i]]$colonyi[[j]] <- colony3[[r[w]]]
    }
}
for(i in 1:d2)
{
  w=w+1
  imp[[nimp]]$colonyi[[i]] <- colony3[[r[w]]]
}
#=========================================================Totla fitness      
nimp <- length(imp)
 k <- c(rep(0,nimp))
for(i in 1:(nimp-1))
{
  for(j in 1:d1)
  { k[i] <- k[i]+imp[[i]]$colonyi[[j]]$fit}
   k[i] <- k[i]/d1  
    imp[[i]]$totalfit=imp[[i]]$fit+zeta*k[i]
}
for(j in 1:d2)
  { k[nimp] <- k[nimp]+imp[[nimp]]$colonyi[[j]]$fit}
   k[nimp] <- k[nimp]/d2   
    imp[[nimp]]$totalfit <- imp[[nimp]]$fit+zeta*k[nimp]
#====================================================
p <- c(rep(0,nimp))
for(i in 1:nimp)
 {p[i] <- imp[[i]]$fit}
m <- min(p)
k=0
for(i in 1:nimp)
{
     if(m==imp[[i]]$fit) 
     {k <- i
      break
      }
}
gimp <- imp[[k]]
#==================================================== main loop algorithm
BEST <- matrix(rep(0,maxiter),ncol=1)
for(iter in 1:maxiter)
{
#=======================================Assimilation
nimp <- length(imp)
for(i in 1:nimp)
{
    ncolony <- length(imp[[i]]$colonyi)    
    for(j in 1:ncolony )
     {   
        d <- imp[[i]]$var-imp[[i]]$colonyi[[j]]$var
        d <- d*runif(nvar)*beta        
        imp[[i]]$colonyi[[j]]$var <- imp[[i]]$colonyi[[j]]$var+d
        for(k in 1:nvar)
        {
          imp[[i]]$colonyi[[j]]$var[k] <- max(imp[[i]]$colonyi[[j]]$var[k],lb[k])
          imp[[i]]$colonyi[[j]]$var[k] <- min(imp[[i]]$colonyi[[j]]$var[k],ub[k])
         }      
        imp[[i]]$colonyi[[j]]$fit <- cost(imp[[i]]$colonyi[[j]]$var)       
     }
}
#=======================================================Revolution
nimp <- length(imp)
for(i in 1:nimp)
{
    ncolony <- length(imp[[i]]$colonyi)
    for(j in 1:ncolony)
    { 
        if( runif(1,0,1)<P_revolve)
        {
         k <- sample(nvar,1)
         d <- ub[k]-lb[k]
         d <- 0.1*runif(1,-1:1)*d        
         imp[[i]]$colonyi[[j]]$var[k] <- imp[[i]]$colonyi[[j]]$var[k]+d 
         for(z in 1:nvar)
         {
          imp[[i]]$colonyi[[j]]$var[z] <- max(imp[[i]]$colonyi[[j]]$var[z],lb[z])
          imp[[i]]$colonyi[[j]]$var[z] <- min(imp[[i]]$colonyi[[j]]$var[z],ub[z])
         }        
         imp[[i]]$colonyi[[j]]$fit <- cost(imp[[i]]$colonyi[[j]]$var)
        }     
    }
}
#========================================================Exchange
nimp <- length(imp)
for(i in 1:nimp)
{
 l <- length(imp[[i]]$colonyi)
 p <- c(rep(0,l))
 for(j in 1:l)
 {
  p[j] <- imp[[i]]$colonyi[[j]]$fit
 }
 value1 <- min(p)
 index=0
  for(j in 1:l)
  {
   if(value1==(imp[[i]]$colonyi[[j]]$fit))  
    {
     index <- j
     break
    }
  }   
    if(value1<imp[[i]]$fit)
    {        
        bestcolony <- imp[[i]]$colonyi[[index]]        
        imp[[i]]$colonyi[[index]]$var <- imp[[i]]$var
        imp[[i]]$colonyi[[index]]$fit <- imp[[i]]$fit        
        imp[[i]]$var <- bestcolony$var
        imp[[i]]$fit <- bestcolony$fit       
    }   
}
#===============================================imperialistic competition 
nimp <- length(imp)
if(nimp>=2)
{
   p <- c(rep(0,nimp))
   for(j in 1:nimp)
     {p[j] <- imp[[j]]$totalfit}
   m <- max(p)
   index1=0
   for(j in 1:nimp)
   {
      if(m==imp[[j]]$totalfit)
      {
        index1 <- j
         break
      }
   }
   wimp <- imp[[index1]]
   l <- length(wimp$colonyi)
   p <- c(rep(0,l))
   for(j in 1:l)
      {p[j] <- wimp$colonyi[[j]]$fit}
   m <- max(p)
   index2=0
   for(j in 1:l)
   {
     if(m==wimp$colonyi[[j]]$fit)
      {
        index2 <- j
         break
      }
   }
   wcolony <- wimp$colonyi[[index2]]
   l <- length(imp[[index1]]$colonyi)
   if(index2==l)
    {length(imp[[index1]]$colonyi)=length(imp[[index1]]$colonyi)-1}
   if(index2!=l)
    {
      l <- l-1
      for(j in index2:l)
       {imp[[index1]]$colonyi[[j]] <- imp[[index1]]$colonyi[[j+1]]}
      length(imp[[index1]]$colonyi) <- length(imp[[index1]]$colonyi)-1
    }
   l2 <- length(imp)
   p <- c(rep(0,l2))
   for(j in 1:l2)
     {p[j] <- imp[[j]]$totalfit}
   p[index1] <- 0
   p <- max(p)-p-2.2204e-16   # eps=2.2204e-16
   p <- p/sum(p)
   p <- cumsum(p)
   d <- runif(1,0,1)
   for(i in 1:l2)
   {
      if((d<p[i]) || (d=p[i]))
      {
       k <- i
       break
      }
    }
   #------------------
   n1 <- length(imp[[k]]$colonyi)
   length(imp[[k]]$colonyi) <- length(imp[[k]]$colonyi)+1
   n1 <- n1+1
   imp[[k]]$colonyi[[n1]]$var <- wcolony$var
   imp[[k]]$colonyi[[n1]]$fit <- wcolony$fit
   n <- length(imp[[index1]]$colonyi)
   if (n==0)
   {
     m <- nimp
     if(index1==m)
       {length(imp) <- length(imp)-1}
     if(index1!=m)
       {
        m <- m-1
        for(j in index1:m)
          {imp[[j]] <- imp[[j+1]]}
        length(imp) <- length(imp)-1
       }
     l3 <- length(imp)
     p <- c(rep(0,l3))
     for(j in 1:l3)
       {p[j] <- imp[[j]]$totalfit}   
     p <- max(p)-p+2.2204e-16
     p <- p/sum(p)
     p <- cumsum(p)
     d <- runif(1,0,1)
     for(i in 1:l3)
     {
       if((d<p[i]) || (d=p[i]))
        {
         k <- i
         break
        }
      }
      #----------------
     n2 <- length(imp[[k]]$colonyi)
     length(imp[[k]]$colonyi) <- length(imp[[k]]$colonyi)+1
     n2 <- n2+1
     imp[[k]]$colonyi[[n2]]$var <- wimp$var
     imp[[k]]$colonyi[[n2]]$fit <- wimp$fit 
   }
   # n==0   
}
# nimp>=2
################################# End of imperialistic competition 
nimp <- length(imp)
p <- c(rep(0,nimp))
for(i in 1:nimp)
 {p[i] <- imp[[i]]$fit}
value2 <- min(p)
for(i in 1:nimp)
{ 
  if(value2==imp[[i]]$fit)
   {
    index <- i
    break
   }
}  
   if(value2<gimp$fit)
   {gimp <- imp[[index]]} 

BEST[iter,1] <- gimp$fit
nimp <- length(imp)
if (nimp==1)
    {
     break
     }

}
##### End
# results algorithm
#BEST
cat("Best Solution=",gimp$var,"\n")
cat("Best Fitness=",gimp$fit,sep = " ", "\n")
cat("Nimp = ",nimp,"\n")
plot(BEST,xlab="Iteration",ylab="Fitness",main="ICA",type="s",col="red")
proc.time() - ptm
}
