Metropolis_N(const K, const LL2, const npp, const lims, const delta,const lambda,  const c, const M_C, const L, const U, const nL, const nU, const SiL, const SiU, const S, const BetaS, const par, const m, const MuS)
{
							 

decl GGk, Lest, I, Alpha, u, i;


GGk=Gk(K,lambda,delta);		  
Lest= delta*max(lambda)-min(lambda);

serial decl M_C_novo=zeros(1,5);
decl npp_novo=npp;

decl npp_aux=zeros(LL2+1,1);
npp_aux[0]=0;

parallel for(i=1;i<=LL2;++i)
{
 npp_aux[i]=sumc(npp[:(i-1)]);
}

decl ar=L*U/LL2;

I=zeros(LL2,1);


decl M_C_at,N,M_C_cand,BetaNcand,Nk,Nc;


parallel for(i=1;i<=LL2;++i)
{

if(npp[i-1]>0){M_C_at=M_C[npp_aux[i-1]:(npp_aux[i]-1)][];}
else{M_C_at=<>;}


N=ranpoisson(1,1,ar*Lest);

if(N>0)
{

M_C_cand=runif(N,1,lims[i-1][0],lims[i-1][1])~runif(N,1,lims[i-1][2],lims[i-1][3])~runif(N,1,0,Lest);

BetaNcand=Sampling_SPG_cond3(M_C_cand[][0:1], L, U, nL, nU, SiL, SiU, S, BetaS, par, m);

if(npp[i-1]>0){Nk=RegiaoK(K,c, M_C_at[][3]);}else{Nk=zeros(1,K);}

Nc=RegiaoK(K,c, BetaNcand);


Alpha= sumr((Nc-Nk).*log(GGk));

u=log(ranu(1,1));


if(u <= Alpha)
{
M_C_novo|=M_C_cand~BetaNcand~(zeros(N,1)+i);
npp_novo[i-1]=N;
I[i-1]=1;
}

else
{
if(npp[i-1]>0){M_C_novo|=M_C_at~(zeros(npp[i-1],1)+i);}
npp_novo[i-1]=npp[i-1];
I[i-1]=0;
}
}


else
{
if(npp[i-1]>0){Nk=RegiaoK(K,c, M_C_at[][3]);}else{Nk=zeros(1,K);}	  

Nc=zeros(1,K);


Alpha= sumr((Nc-Nk).*log(GGk));

u=log(ranu(1,1));


if(u <= Alpha)
{
npp_novo[i-1]=N;
I[i-1]=1;
}

else
{
if(npp[i-1]>0){M_C_novo|=M_C_at~(zeros(npp[i-1],1)+i);}
npp_novo[i-1]=npp[i-1];
I[i-1]=0;
}
}


}


M_C_novo=M_C_novo[1:][];

M_C_novo=sortbyc(M_C_novo,4);

M_C_novo=M_C_novo[][:3];

return {M_C_novo,npp_novo,I};

}



