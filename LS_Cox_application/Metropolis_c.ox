
Metropolis_c(const K, const cAtual, const delta,const lambda, const M_Y, const M_C, const y, const rc){
												   
decl  GGk, cCand, Ncand, Ycand, Alpha, I, cAtual_aceito, Nk, Yk, Lest, MC_at, i;


GGk=Gk(K,lambda,delta);

Nk = RegiaoK(K, cAtual, M_C[][3]); 	 
Yk = RegiaoK(K, cAtual, M_Y[][2]); 	 

cCand = runif(1,K-1, cAtual-rc, cAtual+rc);

if( (K>2 && minr(cCand[1:(K-2)]-cCand[0:(K-3)])>0) || K==2)
{

Ncand= RegiaoK(K, cCand, M_C[][3]); 	 
Ycand= RegiaoK(K, cCand, M_Y[][2]); 	 

Alpha = sumr((Ncand-Nk).*log(GGk)) + sumr((Ycand-Yk).*log(lambda));

if(log(ranu(1,1)) <= Alpha ){I=1;cAtual_aceito=cCand;}
else{I=0;cAtual_aceito=cAtual;}

}

else{I=0;cAtual_aceito=cAtual;}

return {cAtual_aceito, I};
}
