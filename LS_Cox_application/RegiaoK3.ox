RegiaoK3(const K, const cAtual, const M){

decl Y, c, i, Y_ind;

c=(-10^4)~cAtual~(10^4);

Y_ind=zeros(rows(M),1);

for(i=1;i<=K;++i)
{
Y= vecindex(M  .>= c[i-1] .&& M .< c[i]) ;
Y_ind[Y]=i-1;
}

return Y_ind;

}