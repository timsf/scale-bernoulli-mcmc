
leaf_sim_lamb(const i, const k, const ind_flp_lf, const C1pf, const C1p, const nsk, const delta, const ind_lamb, const dif_lamb, const  muj, const ind_lfs, const lims, const scp_prb,
		 const  L, const U, const nL, const nU, const SiL, const SiU, const S, const BetaS, const par, const m, const s2S, const cAtual)
{



decl ind_stop,cont_lf,scp_ind,C1,NP1,NP2,NP3,a1,a2,a3,n1,n2,n3,U1,U2,U3,beta_est,ind_est1,ind_est2,ind_est3,ind_est,out_lfs,j;

decl ccc=0;

cont_lf=0;
ind_stop=0;


//if(ind_flp_lf[i] == 1){out_lfs=ranbinomial(1,1,1,C1pf[i]);}
//else{out_lfs=ranbinomial(1,1,1,C1p[i]);}


while(ccc==0)
{

decl cc=0;

cont_lf+=1;

scp_ind=ranbinomial(1,1,1,scp_prb);

if(scp_ind==1){ccc=1;ind_stop=1;out_lfs=0;}

else
{

if(ind_flp_lf[i] == 1)
{

C1=ranbinomial(1,1,1,C1pf[i]);

if(C1==1 && ind_lamb==0)
{

n1=nsk[i];

NP1=ranpoisson(n1,1,delta*dif_lamb*muj);

a1=vecindex(ind_lfs[][0],i);

//C2p = 1;

decl j=1;
while(j<=n1 && cc==0)
{
if(NP1[j-1]>0)
{
U1=runif(NP1[j-1],1,lims[a1[j-1]][0],lims[a1[j-1]][1]) ~ runif(NP1[j-1],1,lims[a1[j-1]][2],lims[a1[j-1]][3]);
beta_est=Sampling_SPG_cond2(U1, L, U, nL, nU, SiL, SiU, S, BetaS, par, m);
ind_est=RegiaoK4(k, cAtual, beta_est);
if(ind_est==0){cc=1;}
//C2p = C2p*prodc( (delta-(1-ind_est))/delta );
}j+=1;
}

//C2=ranbinomial(1,1,1,exp(C2p));
if(cc==0){ccc=1;out_lfs=1;}
}

else if(C1==1 && ind_lamb==1)
{
ccc=1;out_lfs=1; 
}

else if(C1==0 && ind_lamb==1)
{
n1=nsk[i];

NP1=ranpoisson(n1,1,delta*dif_lamb*muj);

a1=vecindex(ind_lfs[][0],i);

//C2p = 1;

decl j=1;
while(j<=n1 && cc==0)
{
if(NP1[j-1]>0)
{
U1=runif(NP1[j-1],1,lims[a1[j-1]][0],lims[a1[j-1]][1]) ~ runif(NP1[j-1],1,lims[a1[j-1]][2],lims[a1[j-1]][3]);
beta_est=Sampling_SPG_cond2(U1, L, U, nL, nU, SiL, SiU, S, BetaS, par, m);
ind_est=RegiaoK4(k, cAtual, beta_est);
if(ind_est==0){cc=1;}
//C2p = C2p*prodc( (delta-(1-ind_est))/delta );
}j+=1;
}
//C2=ranbinomial(1,1,1,C2p);
if(cc==0){ccc=1;out_lfs=0;}
}

else
{
ccc=1;out_lfs=0; 
}


} // if flp




else{

C1=ranbinomial(1,1,1,C1p[i]);

if(C1==1 && ind_lamb==1)
{
n1=nsk[i];

NP1=ranpoisson(n1,1,delta*dif_lamb*muj);

a1=vecindex(ind_lfs[][0],i);

//C2p = 1;
decl j=1;
while(j<=n1 && cc==0)
{
if(NP1[j-1]>0)
{
U1=runif(NP1[j-1],1,lims[a1[j-1]][0],lims[a1[j-1]][1]) ~ runif(NP1[j-1],1,lims[a1[j-1]][2],lims[a1[j-1]][3]);
beta_est=Sampling_SPG_cond2(U1, L, U, nL, nU, SiL, SiU, S, BetaS, par, m);
ind_est=RegiaoK4(k, cAtual, beta_est);
if(ind_est==1){cc=1;}
//C2p = C2p*prodc( (delta-ind_est)/delta );
}j+=1;
}

//C2=ranbinomial(1,1,1,C2p);
if(cc==0){ccc=1;out_lfs=1;}
}

else if(C1==1 && ind_lamb==0)
{
ccc=1;out_lfs=1; 
}

else if(C1==0 && ind_lamb==0)
{
n1=nsk[i];

NP1=ranpoisson(n1,1,delta*dif_lamb*muj);

a1=vecindex(ind_lfs[][0],i);

//C2p = 1;

decl j=1;
while(j<=n1 && cc==0)
{
if(NP1[j-1]>0)
{
U1=runif(NP1[j-1],1,lims[a1[j-1]][0],lims[a1[j-1]][1]) ~ runif(NP1[j-1],1,lims[a1[j-1]][2],lims[a1[j-1]][3]);
beta_est=Sampling_SPG_cond2(U1, L, U, nL, nU, SiL, SiU, S, BetaS, par, m);
ind_est=RegiaoK4(k, cAtual, beta_est);
if(ind_est==1){cc=1;}
//C2p = C2p*prodc( (delta-ind_est)/delta );
}j+=1;
}

//C2=ranbinomial(1,1,1,C2p);
if(cc==0){ccc=1;out_lfs=0;}
}

else
{
ccc=1;out_lfs=0; 
}

} 	// else flp




}	// else scp_ind=0

} // while ccc



return {ind_stop,out_lfs,cont_lf};

}