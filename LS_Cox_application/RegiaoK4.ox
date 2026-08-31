RegiaoK4(const k, const cAtual, const M){

decl Y, c, i;


c=(-10^4)~cAtual~(10^4);

Y=zeros(rows(M),1);

i = vecindex(M  .>= c[k] .&& M .< c[k+1]);

Y[i]=1;

return Y;

}