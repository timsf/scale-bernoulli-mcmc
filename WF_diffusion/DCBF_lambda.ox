
Metropolis_lambda(const K, const k, const delta, lambdaAtual, const lambdaCand, const cAtual, const M_Y, const L, const U, const nL, const nU, const SiL, const SiU, const S, BetaS,
				  const par, const m, const Iind, const Cmu, const CCov, const s2S, const rc, const lims, const scp_prb, const LL, const NL, const NS,
				  const ind_lfs, const n_flp_lf, const n_std_lf, const Y_ind_sqr, const ind_flp_lf, const alphal, const betal, const cl)

{




decl I_atual;

decl muj=U*L/NS;

serial decl cont_lf=0;

serial decl ind_stop=0;


decl ind_lamb=0;
if(lambdaCand[k]>lambdaAtual[k]){ind_lamb=1;}

if(lambdaCand[k] < 0)
{
I_atual = 0;
}


else{

decl dif_lamb=fabs(lambdaAtual[k]-lambdaCand[k]);

decl Cand_par;
decl Atual_par;
decl i,j,l;
Cand_par=Atual_par=zeros(binomial(K,2),2);
j=0;
for(i=1;i<=(K-1);++i)
{
 for(l=(i+1);l<=K;++l)
{
 Cand_par[j][]=lambdaCand[i-1]~lambdaCand[l-1];
 Atual_par[j][]=lambdaAtual[i-1]~lambdaAtual[l-1];
 j+=1;
 }
}

decl prior_numer = sumc( (alphal-1)*log(lambdaCand)-betal*lambdaCand ) + sumc( log(1-exp(-cl*(fabs(Cand_par[][0]-Cand_par[][1])./sumr(sqrt(Cand_par))).^3)) );
decl prior_denom = sumc( (alphal-1)*log(lambdaAtual)-betal*lambdaAtual ) + sumc( log(1-exp(-cl*(fabs(Atual_par[][0]-Atual_par[][1])./sumr(sqrt(Atual_par))).^3)) );



decl C1=Y_ind_sqr.*log(lambdaAtual[k]);
decl C1c=Y_ind_sqr.*log(lambdaCand[k]);

decl C1n=zeros(NL,1);
decl C1d=zeros(NL,1);

decl nSK=zeros(NL,1);

parallel for (i=1;i<=NL;++i)
{
decl ak=vecindex(ind_lfs,i-1);
nSK[i-1]=rows(ak);
decl a = sumc(C1[ak]-C1c[ak]) + (1/NL)*prior_denom - (1/NL)*prior_numer; if(a==<>){a=0;}

C1n[i-1] = 1;
C1d[i-1] = 1 + exp(a);
}


decl C1p=C1n./C1d;



C1=muj*lambdaAtual[k] - Y_ind_sqr.*log(lambdaAtual[k]);
C1c=muj*lambdaCand[k] - Y_ind_sqr.*log(lambdaCand[k]);

C1n=zeros(NL,1);
C1d=zeros(NL,1);

parallel for (i=1;i<=NL;++i)
{
decl ak=vecindex(ind_lfs,i-1);
decl a = sumc(C1c[ak]-C1[ak]) - (-1/NL)*prior_denom + (-1/NL)*prior_numer;  if(a==<>){a=0;}

C1n[i-1] = 1;
C1d[i-1] = 1 + exp(a);
}


decl C1pf=C1n./C1d;

//decl temp=(1-ind_flp_lf).*C1p + ind_flp_lf.*C1pf;
//println(minc(temp)~maxc(temp));

serial decl cc=0;

serial decl out_lfs=zeros(NL,1);

parallel for(i=1;i<=NL;++i)
{
 if(ind_stop==0)
 {
  decl ind_stop_temp,cont_lf_temp;
  [ind_stop_temp,out_lfs[i-1],cont_lf_temp] = leaf_sim_lamb(i-1, k, ind_flp_lf,C1pf,C1p,nSK,delta,ind_lamb,dif_lamb,muj,ind_lfs,lims,
  													scp_prb,L,U,nL,nU,SiL,SiU,S,BetaS,par,m,s2S,cAtual);
  cont_lf+=cont_lf_temp;
  if(ind_stop_temp==1){ind_stop=1;}
  }
}

//println(out_lfs);

while(cc==0)
{

if(ind_stop==1){cc=1;I_atual=0;}

else
{

if(out_lfs==0){I_atual=0;cc=1;}

else if(out_lfs==1){I_atual=1;cc=1;}

else
{

decl lfs_sim=vecindex(DCBF(out_lfs)); //println(lfs_sim);
decl n_lfs_temp=rows(lfs_sim);

parallel for(i=1;i<=n_lfs_temp;++i)
{
 if(ind_stop==0)
 {
  decl ind_stop_temp,cont_lf_temp;
  [ind_stop_temp,out_lfs[lfs_sim[i-1]],cont_lf_temp] = leaf_sim_lamb(lfs_sim[i-1], k, ind_flp_lf,C1pf,C1p,nSK,delta,ind_lamb,dif_lamb,muj,ind_lfs,lims,
  													scp_prb,L,U,nL,nU,SiL,SiU,S,BetaS,par,m,s2S,cAtual);
  cont_lf+=cont_lf_temp;
  if(ind_stop_temp==1){ind_stop=1;}
  }
}

//println(out_lfs);

}  // else out_lfs not all 0 and not all 1
}  // else ind_stop=0


}  // while cc=0

}  // else lambdaCand[k] < 0



if(I_atual==1)
{

lambdaAtual=lambdaCand;

}






return {I_atual,lambdaAtual,cont_lf,ind_stop};
}