
GBarker2(const XXX1, const XXX2, const dtt, const g1o, const g2o, const sig, const xLrow, const xRrow, const cpo, const s1, const s2, const MINI, const MAXI, const X0, const XT, const NS, const NNL, const ll, const scp_prb, const PR_g1_a,const PR_g1_b,const PR_g2_a,const PR_g2_b)

{
  decl S,gamao,gaman,g1n,g2n,n,ind0,ind1,cp,cphi,c1n,c2n,c3n,c4n,c1o,c2o,c3o,
  c4o,c1,c2,c3,c4,x,y,m,M,l1,l21,l22,mm,MM,l,u,C1,i,cB,np,X11,
  indm,sindm,up,phil,phiu,ind2,indX2,X2i,ind3,upp,xt,phix,k1,k2,gama,mm1,mm2,
  an,bn,ao,bo,cpp,a,b,dAn,dAo,cppp,temp;

serial decl X1=XXX1;
serial decl X2=XXX2;
serial decl NL=NNL;
  
gamao=g1o|g2o;


gaman=((g1o|g2o)-(s1|s2)+2*(s1|s2).*ranu(2,1));

g1n=gaman[0];
g2n=gaman[1];

n=rows(X1);

// A-endpoint coefficients (sigma-general: g1->g1/sig^2, g1g2->g1g2/sig^2)
an=-0.5*(g1n/sig^2-1);
bn=g1n*g2n/(2*sig^2)-1/4;

ao=-0.5*(g1o/sig^2-1);
bo=g1o*g2o/(2*sig^2)-1/4;

// phi-bound coefficients (sigma-general) for ind0/ind1
decl kmn,kpn,kmo,kpo,al,be,pp,qq,AA,BB;
al=g1n*g2n; be=g1n*(1-g2n); pp=al-be; qq=al+be-sig^2;
AA=(pp^2+qq^2)/(8*sig^2)-qq/4; BB=pp*qq/(4*sig^2)-pp/4; kpn=(AA+BB)/4; kmn=(AA-BB)/4;
al=g1o*g2o; be=g1o*(1-g2o); pp=al-be; qq=al+be-sig^2;
AA=(pp^2+qq^2)/(8*sig^2)-qq/4; BB=pp*qq/(4*sig^2)-pp/4; kpo=(AA+BB)/4; kmo=(AA-BB)/4;


ind0=0;
if( ( kmo - kmn )>0 ){ind0=1;}
ind1=0;
if( ( kpo - kpn )>0 ){ind1=1;}


if((ind0==0 && ind1==0) || (ind0==1 && ind1==1))
{
cp=nrphi2(g1o,g2o,g1n,g2n,sig,M_PI/(2*sig),10^(-10));
cphi=phi2(g1o,g2o,g1n,g2n,sig,cp);
}
else{cp=M_PI/(2*sig);cphi=0;}


serial decl XX1=Intii20(X1[][:7], dtt, g1o, g2o, g1n, g2n, sig, cp, cphi, xLrow, xRrow);


dAn=(1/NS)*(2*an*log(1/cos(XT/2)) + 2*bn*log(tan(XT/2)) - 2*an*log(1/cos(X0/2)) - 2*bn*log(tan(X0/2)));
dAo=(1/NS)*(2*ao*log(1/cos(XT/2)) + 2*bo*log(tan(XT/2)) - 2*ao*log(1/cos(X0/2)) - 2*bo*log(tan(X0/2)));

decl prn=(1/NS)*( (PR_g1_a-1)*log(g1n)-PR_g1_b*g1n + (PR_g2_a-1)*log(g2n)+(PR_g2_b-1)*log(1-g2n) );
decl pro=(1/NS)*( (PR_g1_a-1)*log(g1o)-PR_g1_b*g1o + (PR_g2_a-1)*log(g2o)+(PR_g2_b-1)*log(1-g2o) );


//serial decl MMM=XX1[][2];
decl s0=XX1[][3];

serial decl ind_lfs=rsample(NS,1,NL);

serial decl p_trms=dAo-dAn+s0+pro-prn;

// ---- [DIAGNOSTIC] factory balance (comment out this block for production) ----
//{
// decl DGintm=-sumc(XX1[][3]);
// decl DGpts =dtt*sumc(XX1[][2]);
// decl DGlogR=NS*(dAn-dAo+prn-pro);          // analytic log-ratio (new/old)
// decl DGsumc=sumc(p_trms);                  // = -DGlogR - DGintm
// decl DGleaf=DGsumc/NL;
// decl DGp   =1/(1+exp(DGleaf));
// decl DGll  =ceil(log(fabs(DGintm)+1e-12)/log(2));
 //println("  [GBarker2] g1 ",g1o," -> ",g1n,"   g2 ",g2o," -> ",g2n,
 //        "  | logRan=",DGlogR,"  int_m=",DGintm,"  sumc=",DGsumc);
 //println("             pts/pass=",DGpts,"  | NL=",NL,"  per-leaf expo=",DGleaf,
 //        "  1st-coin p=",DGp,"  | suggested ll=",DGll);
//}
// ---- [end DIAGNOSTIC] ----

serial decl ind_lfs_set = new array[NL];
serial decl ind_lfs_set2 = new array[NL];



parallel for (i=1;i<=NL;++i)
{
decl iind=vecindex(ind_lfs,i-1);
ind_lfs_set[i-1]=iind;

decl nn=rows(iind);
decl indX2=zeros(1,1); decl indX2_temp,j;
for(j=1;j<=nn;++j)
{
 indX2_temp=vecindex(X2[][0],X1[iind[j-1]][0]);
 indX2_temp=indX2_temp[0]+vecindex(X2[indX2_temp][1],X1[iind[j-1]][1]);
 indX2=indX2|indX2_temp;
}
indX2=indX2[1:];
ind_lfs_set2[i-1]=indX2;
}



serial decl out_lfs=zeros(NL,1);
decl I_atual;
serial decl cont_w=0;
serial decl cont_phi=0;
serial decl ind_stop=0;

serial decl XX2 = zeros(1,columns(X2));		  //println(NL);println(ind_lfs);println(ind_lfs_set);println(ind_lfs_set2);

parallel for(i=1;i<=NL;++i)
{
  decl iind2=ind_lfs_set2[i-1];
  decl X2_temp=X2[iind2][];
  if(ind_stop==0)
 {
  decl ind_stop_temp,cont_lf_temp,out_lfs_temp,X2_temp2;
  decl iind=ind_lfs_set[i-1];
  decl X1_temp=X1[iind][];
  decl XX1_temp=XX1[iind][];
  //decl MMM_temp=MMM[iind];
  decl p=p_trms[iind];
  p=1/(1 + exp(sumc(p)));
  [ind_stop_temp,out_lfs_temp,cont_lf_temp,X2_temp2] = leaf_sim_gamma(p,dtt,X1_temp,X2_temp,XX1_temp,g1o,g2o,g1n,g2n,sig,xLrow[iind],xRrow[iind],scp_prb);
  out_lfs[i-1]=out_lfs_temp;
  cont_phi=cont_phi+cont_lf_temp;
  cont_w=cont_w+1;
  X2_temp=X2_temp2;
  if(ind_stop_temp==1){ind_stop=1;}
 }

 serial{XX2=XX2|X2_temp;}
  
}

XX2=XX2[1:][];
XX2=sortbyc(XX2,<0,1,2>);

X2=XX2;

//println(X2);

serial decl cc=0;

while(cc==0)
{
//println("try");
if(ind_stop==1){cc=1;I_atual=0;cpp=cpo;}

else
{

if(out_lfs==0){I_atual=0;cc=1;cpp=cpo;}

else if(out_lfs==1){I_atual=1;cc=1;cpp=nrphi(g1n,g2n,sig,M_PI/(2*sig),10^(-10));}

else
{

serial decl lfs_sim=vecindex(DCBF(out_lfs)); //println(lfs_sim);
serial decl n_lfs_temp=rows(lfs_sim);


parallel for (i=1;i<=NL;++i)
{
decl iind=vecindex(ind_lfs,i-1);

decl nn=rows(iind);
decl indX2=zeros(1,1); decl indX2_temp,j;
for(j=1;j<=nn;++j)
{
 indX2_temp=vecindex(X2[][0],X1[iind[j-1]][0]);
 indX2_temp=indX2_temp[0]+vecindex(X2[indX2_temp][1],X1[iind[j-1]][1]);
 indX2=indX2|indX2_temp;
}
indX2=indX2[1:];
ind_lfs_set2[i-1]=indX2;
}


serial decl XX2 = zeros(1,columns(X2));

parallel for(i=1;i<=n_lfs_temp;++i)
{
  decl iind2=ind_lfs_set2[lfs_sim[i-1]];
  decl X2_temp=X2[iind2][];
  if(ind_stop==0)
 {
  decl ind_stop_temp,cont_lf_temp,out_lfs_temp,X2_temp2;
  decl iind=ind_lfs_set[lfs_sim[i-1]];
  decl X1_temp=X1[iind][];
  decl XX1_temp=XX1[iind][];
  //decl MMM_temp=MMM[iind];
  decl p=p_trms[iind];
  p=1/(1 + exp(sumc(p)));
  [ind_stop_temp,out_lfs_temp,cont_lf_temp,X2_temp2] = leaf_sim_gamma(p,dtt,X1_temp,X2_temp,XX1_temp,g1o,g2o,g1n,g2n,sig,xLrow[iind],xRrow[iind],scp_prb);
  out_lfs[lfs_sim[i-1]]=out_lfs_temp;
  cont_phi=cont_phi+cont_lf_temp;
  cont_w=cont_w+1;
  X2_temp=X2_temp2;
  if(ind_stop_temp==1){ind_stop=1;}
 }
  serial{XX2=XX2|X2_temp;}
}

if(n_lfs_temp<NL)
{
 serial decl ind_temp=dropr(range(0,NL-1)',lfs_sim);
 serial decl n_lfs_temp2=rows(ind_temp);
 parallel for(i=1;i<=n_lfs_temp2;++i)
 {
  decl iind2=ind_lfs_set2[ind_temp[i-1]];
  decl X2_temp=X2[iind2][];
  serial{XX2=XX2|X2_temp;}
 }
}

XX2=XX2[1:][];
XX2=sortbyc(XX2,<0,1,2>);

X2=XX2;
//println(out_lfs);

}  // else out_lfs not all 0 and not all 1
}  // else ind_stop=0


}  // while cc=0


if(I_atual==0){gama=gamao;}
else{gama=gaman;}

return {X2,gama,cpp,I_atual,ind_stop,cont_phi,cont_w};


}
