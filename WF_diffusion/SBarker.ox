
SBarker(const XXX1, const XXX2, const sigo, const g1, const g2, const s_sig, const zpairs, const x0, const xT, const dt, const dtt, const T, const NS, const NNL, const ll, const scp_prb, const PR_s_a,const PR_s_b)
{
  decl sign,cp,cpp,nn,ii,rr,al,be,Sz,anAo,bnAo,anAn,bnAn,Aendo,Aendn,logRan,m,I_atual,sigret,cpret;

serial decl X1=XXX1;
serial decl X2=XXX2;
serial decl NL=NNL;

nn=dt/dtt;

// ---- propose sigma' (uniform RW); reject immediately if non-positive ----
sign=sigo - s_sig + 2*s_sig*ranu(1,1);
if(sign<=0){ return {X1,X2,sigo,nrphi(g1,g2,sigo,M_PI/(2*sigo),10^(-10)),0}; }

cp =nrphi(g1,g2,sigo,M_PI/(2*sigo),10^(-10));
cpp=nrphi(g1,g2,sign,M_PI/(2*sign),10^(-10));

// per-row sigma=1 endpoints (for leaf_sim_sig)
serial decl zLrow=vec((zpairs[][0]*ones(1,nn))');
serial decl zRrow=vec((zpairs[][1]*ones(1,nn))');

// ---- sigma-difference bounds, interval by interval ----
serial decl XX1=zeros(NS,4);
for(ii=1;ii<=T;++ii)
{
 rr=((ii-1)*nn)|(ii*nn-1);
 XX1[rr[0]:rr[1]][]=Intii20sig(X1[rr[0]:rr[1]][4:7],dtt,g1,g2,sign,sigo,cpp,cp,zpairs[ii-1][0],zpairs[ii-1][1]);
}

// ---- analytic log-ratio  log pi(sign)/pi(sigo)  (IG prior + Jacobian + N-term + telescoped A) ----
m=T;
al=g1*g2; be=g1*(1-g2);
Sz=sumc( (zpairs[][1]-zpairs[][0]).^2 );
anAo=-0.5*(g1/sigo^2-1); bnAo=g1*g2/(2*sigo^2)-1/4;
anAn=-0.5*(g1/sign^2-1); bnAn=g1*g2/(2*sign^2)-1/4;
Aendo=2*anAo*(log(1/cos(xT/2))-log(1/cos(x0/2))) + 2*bnAo*(log(tan(xT/2))-log(tan(x0/2)));
Aendn=2*anAn*(log(1/cos(xT/2))-log(1/cos(x0/2))) + 2*bnAn*(log(tan(xT/2))-log(tan(x0/2)));
logRan = -(2*PR_s_a+1)*log(sign/sigo) - PR_s_b*(1/sign^2-1/sigo^2)
         - m*log(sign/sigo)
         - (Sz/2)*(1/sign^2-1/sigo^2)
         + (Aendn-Aendo);

// per-row first-coin log-terms:  -logRan/NS  distributed, plus the interval lower-bound part I=-dtt*l
serial decl s0=XX1[][3];
serial decl p_trms= constant(-logRan/NS,NS,1) + s0;

// ---- [DIAGNOSTIC] factory balance (comment out this block for production) ----
// int_m  : integral of the lower bound of Delta phi  (the term that saturates the 1st coin)
// sumc   : total first-coin exponent  = -logRan - int_m   (want |sumc/NL| ~ O(1))
// p_leaf : first-coin accept prob at the mean per-leaf load
// ll_rec : leaf depth needed so the per-leaf load is O(1)  (set ll >= ll_rec)
//{
// decl DGintm=-sumc(XX1[][3]);
// decl DGpts =dtt*sumc(XX1[][2]);
// decl DGsumc=sumc(p_trms);
// decl DGleaf=DGsumc/NL;
// decl DGp   =1/(1+exp(DGleaf));
// decl DGll  =ceil(log(fabs(DGintm)+1e-12)/log(2));
// println("  [SBarker]  sig ",sigo," -> ",sign,
//         "  | logRan=",logRan,"  int_m=",DGintm,"  sumc(=-logRan-int_m)=",DGsumc);
// println("             pts/pass=",DGpts,"  | NL=",NL,"  per-leaf expo=",DGleaf,
//         "  1st-coin p=",DGp,"  | suggested ll=",DGll);
//}
// ---- [end DIAGNOSTIC] ----

// ================= DCBF / leaf machinery (mirrors GBarker2) =================
serial decl ind_lfs=rsample(NS,1,NL);
serial decl ind_lfs_set  = new array[NL];
serial decl ind_lfs_set2 = new array[NL];

decl i;
parallel for (i=1;i<=NL;++i)
{
decl iind=vecindex(ind_lfs,i-1);
ind_lfs_set[i-1]=iind;
decl nnr=rows(iind);
decl indX2=zeros(1,1); decl indX2_temp,j;
for(j=1;j<=nnr;++j)
{
 indX2_temp=vecindex(X2[][0],X1[iind[j-1]][0]);
 indX2_temp=indX2_temp[0]+vecindex(X2[indX2_temp][1],X1[iind[j-1]][1]);
 indX2=indX2|indX2_temp;
}
indX2=indX2[1:];
ind_lfs_set2[i-1]=indX2;
}

serial decl out_lfs=zeros(NL,1);
serial decl ind_stop=0;

serial decl XX2 = zeros(1,columns(X2));
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
  decl p=p_trms[iind];
  p=1/(1 + exp(sumc(p)));
  [ind_stop_temp,out_lfs_temp,cont_lf_temp,X2_temp2] = leaf_sim_sig(p,dtt,X1_temp,X2_temp,XX1_temp,g1,g2,sign,sigo,zLrow[iind],zRrow[iind],scp_prb);
  out_lfs[i-1]=out_lfs_temp;
  X2_temp=X2_temp2;
  if(ind_stop_temp==1){ind_stop=1;}
 }
 serial{XX2=XX2|X2_temp;}
}
XX2=XX2[1:][];
XX2=sortbyc(XX2,<0,1,2>);
X2=XX2;

serial decl cc=0;
while(cc==0)
{
if(ind_stop==1){cc=1;I_atual=0;}
else
{
if(out_lfs==0){I_atual=0;cc=1;}
else if(out_lfs==1){I_atual=1;cc=1;}
else
{
serial decl lfs_sim=vecindex(DCBF(out_lfs));
serial decl n_lfs_temp=rows(lfs_sim);

parallel for (i=1;i<=NL;++i)
{
decl iind=vecindex(ind_lfs,i-1);
decl nnr=rows(iind);
decl indX2=zeros(1,1); decl indX2_temp,j;
for(j=1;j<=nnr;++j)
{
 indX2_temp=vecindex(X2[][0],X1[iind[j-1]][0]);
 indX2_temp=indX2_temp[0]+vecindex(X2[indX2_temp][1],X1[iind[j-1]][1]);
 indX2=indX2|indX2_temp;
}
indX2=indX2[1:];
ind_lfs_set2[i-1]=indX2;
}

serial decl XX2b = zeros(1,columns(X2));
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
  decl p=p_trms[iind];
  p=1/(1 + exp(sumc(p)));
  [ind_stop_temp,out_lfs_temp,cont_lf_temp,X2_temp2] = leaf_sim_sig(p,dtt,X1_temp,X2_temp,XX1_temp,g1,g2,sign,sigo,zLrow[iind],zRrow[iind],scp_prb);
  out_lfs[lfs_sim[i-1]]=out_lfs_temp;
  X2_temp=X2_temp2;
  if(ind_stop_temp==1){ind_stop=1;}
 }
  serial{XX2b=XX2b|X2_temp;}
}
if(n_lfs_temp<NL)
{
 serial decl ind_temp=dropr(range(0,NL-1)',lfs_sim);
 serial decl n_lfs_temp2=rows(ind_temp);
 parallel for(i=1;i<=n_lfs_temp2;++i)
 {
  decl iind2=ind_lfs_set2[ind_temp[i-1]];
  decl X2_temp=X2[iind2][];
  serial{XX2b=XX2b|X2_temp;}
 }
}
XX2b=XX2b[1:][];
XX2b=sortbyc(XX2b,<0,1,2>);
X2=XX2b;
}
}
}

// ---- accept / reject ----
if(I_atual==0){ sigret=sigo; cpret=cp; }
else
{
 sigret=sign; cpret=cpp;
 // refresh X1 Inti0 columns (8:11) at the accepted sigma
 for(ii=1;ii<=T;++ii)
 {
  rr=((ii-1)*nn)|(ii*nn-1);
  X1[rr[0]:rr[1]][8:]=Inti0(X1[rr[0]:rr[1]][4:7],dtt,g1,g2,sign,cpp,zpairs[ii-1][0]/sign,zpairs[ii-1][1]/sign);
 }
}

return {X1,X2,sigret,cpret,I_atual,ind_stop};
}
