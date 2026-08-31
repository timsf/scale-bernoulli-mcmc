
Gk(const K, const lambdaAtual,const delta){

decl LM, Ln, G, a;


LM = max(lambdaAtual);
Ln = min(lambdaAtual);
a=delta*LM -Ln;

G = (delta*LM - lambdaAtual)/a; 



return G;

}
