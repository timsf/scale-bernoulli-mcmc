COV_GP(const dist, const par)

{

decl t2,S;

t2=par;

S = exp( -( (dist).^1.95) / (2*t2) ); // specify the covariance function

return S;

}





