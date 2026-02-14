x1=runif(300)
x=matrix(x1,100)
gene.data1=function(x,nz,n,p)
{
      set.seed(x)
      beta<-1/sqrt(3)*c(1, 0, 0)
sigma2<-0.01             #true sigma2
A<-sqrt(3)/2-1.645/sqrt(12); C<-sqrt(3)/2+1.645/sqrt(12);
x1<-runif(p*n, min=0, max=1);
X<-matrix(x1,ncol=n)
index<-matrix(t(beta)%*%X,ncol=1,byrow=T);
epsilon<-matrix(rnorm(n,mean=0,sd=sqrt(sigma2)),ncol=1,byrow=T) #epsilon is 1*n
y<-matrix(sin(pi*(index-A)/(C-A))+epsilon,ncol=1,byrow=T)
      return(list(X=X,y=y))
}


 
x=gene.data1(x,3,100,3)$X
y=gene.data1(x,3,100,3)$y
