
dphi(const g1, const g2, const sig, const x)   // = phi'(x)  (exact first derivative)

{
  decl al,be,p,q,Ac,Bc,cp,cm,s,c,f;

al=g1*g2;
be=g1*(1-g2);
p=al-be;
q=al+be-sig^2;

Ac=(p^2+q^2)/(8*sig^2) - q/4;
Bc=p*q/(4*sig^2) - p/4;

cp=(Ac+Bc)/4;
cm=(Ac-Bc)/4;

s=sin(sig*x/2);
c=cos(sig*x/2);

f= sig*( cm*s/c^3 - cp*c/s^3 );

return f;


}
