
nrphi(const g1, const g2, const sig, const x0, const eps)

{
  decl dif,f,df,x,it;


dif=1;
x=x0;
it=0;

while(dif>eps && it<1000)      // it<1000 guard (exact Newton pair converges fast)
{

f=dphi(g1,g2,sig,x);
df=ddphi(g1,g2,sig,x);

dif=f/df;

x=x-dif;

dif=fabs(dif);

it+=1;

}

  
return x;


}
