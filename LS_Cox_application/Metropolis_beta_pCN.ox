  																											
Metropolis_beta_pCN(const K, const delta, const lambda, const c, const M_C, const M_Y, const w, const L, const U, const nL, const nU, const SiL, const SiU, const S, const BetaS, const par, const m, const Iind, const Cmu, const CCov, const s2S)

{

decl Lest,eps,i,BetaS_cand,BetaYcand,BetaNcand,GGk,Yk,Nk,Yc,Nc,mu_cand,rtprior,Alpha,I_atual,M_C_novo,M_Y_novo,BetaS_novo,u,mu_novo,cnovo;


decl nn=rows(S)-1;

eps=zeros(nn+1,1);

eps[0]=rann(1,1);

for(i=1;i<=nn;++i)
{
eps[i]=Cmu[i]*eps[Iind[i]] + sqrt(CCov[i])*rann(1,1);
}


BetaS_cand=sqrt(1-s2S)*BetaS+sqrt(s2S)*eps;


decl epsYcand=Sampling_SPG_cond2(M_Y[][0:1], L, U, nL, nU, SiL, SiU, S, eps, par, m);
decl epsNcand=Sampling_SPG_cond2(M_C[][0:1], L, U, nL, nU, SiL, SiU, S, eps, par, m);


BetaYcand=sqrt(1-s2S)*M_Y[][2]+sqrt(s2S)*epsYcand;
BetaNcand=sqrt(1-s2S)*M_C[][3]+sqrt(s2S)*epsNcand;


//decl cCand = c;

GGk=Gk(K,lambda,delta);

Yk=RegiaoK(K, c, M_Y[][2]); 	 
Nk=RegiaoK(K, c, M_C[][3]); 	  

Yc=RegiaoK(K, c, BetaYcand);
Nc=RegiaoK(K, c, BetaNcand);


mu_cand=zeros(nn+1,1);
mu_cand[0]=0;


Alpha= sumr((Nc-Nk).*log(GGk)) + sumr((Yc-Yk).*log(lambda));


I_atual=0;
u=log(ranu(1,1));
if(u <= Alpha)
{
I_atual=1;
M_Y_novo= M_Y[][0:1]~ BetaYcand;
M_C_novo= M_C[][0:2] ~ BetaNcand;
BetaS_novo=BetaS_cand;
cnovo=c;


}

else{M_C_novo=M_C ; M_Y_novo=M_Y; BetaS_novo=BetaS; cnovo=c;}

decl temp=RegiaoK(K, c, BetaS_novo);
decl n=rows(BetaS_novo);
Yc=RegiaoK(K, c, M_Y_novo[][2]);
decl logVer =  sumr(-(temp/n).*lambda) + sumr(Yc.*log(lambda));


return {M_Y_novo,M_C_novo,BetaS_novo,I_atual,logVer};	   

}