Sampling_SPG_cond3(const loc, const L, const U, const nL, const nU, const SiL, const SiU, const S, const BetaS, const par, const m)

{

decl n,Beta_loc,i,neig,D,SS,S12,S22,S22_inv,mu,var;


n=rows(loc);

Beta_loc=zeros(n,1);


for(i=1; i<=n; ++i)
{

 neig=locc(loc[i-1][], nL, nU, SiL, SiU, S, L, U, m);

 D=dmatrix(loc[i-1][]|S[neig][]);  

 SS=COV_GP(D, par);
 SS=(SS+SS')/2;

 S12=  SS[0][1:];	  
 S22=  SS[1:][1:];   

 S22_inv=invert(S22);

 SS=S12*(S22_inv);

 mu=SS*BetaS[neig]; 
 var= 1-SS*S12';

 Beta_loc[i-1]=mu+sqrt(var)*rann(1,1);
 
}

return  Beta_loc;	

}
