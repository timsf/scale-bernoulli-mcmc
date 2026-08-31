
leaf_sim_gamma(const p, const dtt, const X1, X2, const XX1, const g1o, const g2o, const g1n, const g2n, const sig, const xLrow, const xRrow, const scp_prb)
{

decl ind_stop,cont_lf,scp_ind,C1,out_lfs,i,cB,np,X11,up,indX2,X2i,xt,phix,k1,k2,Lval;

decl ccc=0;
cont_lf=0;
ind_stop=0;
decl n=rows(X1);
decl MMM=XX1[][2];

while(ccc==0)
{
decl cc=0;
cont_lf+=1;
scp_ind=ranbinomial(1,1,1,scp_prb);
if(scp_ind==1){ccc=1;ind_stop=1;out_lfs=0;}
else
{
  C1=ranbinomial(1,1,1,p);
   if(C1==1)
 {
  ccc=1; out_lfs=1;
 }
 else
 {
   i=1;
  cB=1;
  while (cB==1 && i<=n)	
  {
  np=ranpoisson(1,1,dtt*MMM[i-1]);
  if(np==0){i+=1;}
  else
  {
   X11=X1[i-1][:7];
   up=runif(np,2,zeros(np,2),(zeros(np,1)+dtt)~(zeros(np,1)+MMM[i-1]));
   up=sortbyc(up,0);
   indX2=vecindex(X2[][0],X11[0][0]);
   indX2=indX2[0]+vecindex(X2[indX2][1],X11[0][1]);
   X2i=X2[indX2][];
   [X2i,xt]=BBx(X11,X2i,up[][0],dtt);
   Lval=xLrow[i-1]+(xRrow[i-1]-xLrow[i-1])*(X1[i-1][1]*dtt+up[][0])+xt;
   phix=phi2(g1o,g2o,g1n,g2n,sig,Lval)-XX1[i-1][0];
   if(up[][1]>phix){i+=1;}  
   else{cB=0;}
   k1=rows(indX2);	
   k2=rows(X2i);
   X2=insertr(X2,indX2[0],k2);
   X2=dropr(X2,range((indX2[0]+k2),(indX2[0]+k2+k1-1),1));
   X2[indX2[0]:(indX2[0]+k2-1)][]=X2i;
   } //else np=0
  }	//while(cB==1)
if(cB==1){ccc=1;out_lfs=0;}
 } //else C1=1
}	// else scp_ind=0
} // while ccc

return {ind_stop,out_lfs,cont_lf,X2};
}
