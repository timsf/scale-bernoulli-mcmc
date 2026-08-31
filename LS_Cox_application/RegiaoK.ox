RegiaoK(const K, const cAtual, const M){

decl Y, c, i;

Y=zeros(1,K);

c=(-10^4)~cAtual~(10^4);

for(i=1;i<=K;++i)
{
Y[i-1] = rows( vecindex(M  .>= c[i-1] .&& M .<= c[i]) );
}


return Y;

}