
DCBF_lambda(const K, const k, const delta, lambdaAtual, const lambdaCand, const cAtual, const M_Y, const L, const U, const nL, const nU, const SiL, const SiU, const S, BetaS,
				  const par, const m, const Iind, const Cmu, const CCov, const s2S, const lims, const scp_prb, const LL, const NL, const NS,
				  const ind_lfs, const n_flp_lf, const n_std_lf, const Y_ind_sqr, const ind_flp_lf)

{




decl I_atual;

decl muj=U*L/NS;

serial decl cont_lf=0;

serial decl ind_stop=0;

serial decl cont_w=0;
serial decl cont_phi=0;

decl ind_lamb=0;
if(lambdaCand[k]>lambdaAtual[k]){ind_lamb=1;}

if(lambdaCand[k] < 0)
{
I_atual = 0;
}


else{

decl dif_lamb=fabs(lambdaAtual[k]-lambdaCand[k]);

decl prior_numer = 0;
decl prior_denom = 0;


decl C1=Y_ind_sqr.*log(lambdaAtual[k]);
decl C1c=Y_ind_sqr.*log(lambdaCand[k]);

decl C1n=zeros(NL,1);
decl C1d=zeros(NL,1);

decl nSK=zeros(NL,1);

decl i;
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


serial decl cc=0;

serial decl out_lfs=zeros(NL,1);


parallel for(i=1;i<=NL;++i)
{
 if(ind_stop==0)
 {
  decl ind_stop_temp,cont_lf_temp;
  [ind_stop_temp,out_lfs[i-1],cont_lf_temp] = leaf_sim_lamb(i-1, k, ind_flp_lf,C1pf,C1p,nSK,delta,ind_lamb,dif_lamb,muj,ind_lfs,lims,
  													scp_prb,L,U,nL,nU,SiL,SiU,S,BetaS,par,m,s2S,cAtual);
  cont_phi=cont_phi+cont_lf_temp;
  cont_w=cont_w+1;
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
  cont_phi=cont_phi+cont_lf_temp;
  cont_w=cont_w+1;
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






return {I_atual,lambdaAtual,ind_stop,cont_phi,cont_w};
}