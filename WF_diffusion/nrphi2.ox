
nrphi2(const g1n, const g2n, const g1o, const g2o, const sig, const x0, const eps)

{
  decl dif,f,df,x,it;

dif=1; x=x0; it=0;

while(dif>eps && it<1000)
{

f=dphi2(g1n,g2n,g1o,g2o,sig,x);
df=ddphi2(g1n,g2n,g1o,g2o,sig,x);

dif=f/df;

x=x-dif;

dif=fabs(dif);

it+=1;

}

  
return x;


}
