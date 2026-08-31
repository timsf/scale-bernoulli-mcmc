Rk(const K, const c, const S, const BetaS, const M_Y, const lambdaAtual)
{

decl m, a, aa, S1, S2, Yk, muk, D, lv, i;


m=rows(BetaS);

S1=zeros(m,1);
S2=zeros(m,K);

a=zeros(1,K);

decl cn=(-10^4)~c~(10^4);


for(i=1;i<=K;++i)
{
aa = vecindex(BetaS  .>= cn[i-1] .&& BetaS .<= cn[i]);
S1[aa] = lambdaAtual[i-1];
a[i-1]=rows(aa);
S2[aa][i-1] = 1;
}


muk=(a/m)*100;

Yk=RegiaoK(K, c, M_Y);

lv= -sumr(muk.*lambdaAtual) + sumr(Yk.*log(lambdaAtual));

D=-2*lv;


return{S1, S2, Yk, muk, D, lv};

}

