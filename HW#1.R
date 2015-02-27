myqueue=function(lambda,mu,n,p_num,c){#n:length of loop, p_num: which level you want to estimate, c: number you want for simulation
  w_p=c()
  for(i in 1:n){ # generate w
    t=rexp(100,lambda)
    s=rexp(100,mu)
    d=c()
    a=c()
    w=c()
    w[1]=0
    a[1]=t[1]
    d[1]=a[1]+w[1]+s[1]
    for (j in 2:100){
      a[j]=sum(t[1:j])
      w[j]=max(0,d[j-1]-a[j])
      d[j]=a[j]+w[j]+s[j]
    }
    w_p[i]=w[p_num] # generate w_p[n], iid of w_p
  }
  z=sum(w_p>c) 
  m_ba=z/n
  sigma=sqrt((sum((z-m_ba)^2))/(n-1))
  E=m_ba
  CI=c(m_ba-1.96*(sigma/sqrt(n)),m_ba+1.96*(sigma/sqrt(n)))
  return(list(E=E,CI=CI,w_p=w_p))
}



