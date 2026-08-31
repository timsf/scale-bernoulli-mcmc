
ddphi(const g1, const g2, const sig, const x)  // = phi''(x)  (exact second derivative)

{
  decl al,be,p,q,Ac,Bc,cp,cm,Y,f;

al=g1*g2;
be=g1*(1-g2);
p=al-be;
q=al+be-sig^2;

Ac=(p^2+q^2)/(8*sig^2) - q/4;
Bc=p*q/(4*sig^2) - p/4;

cp=(Ac+Bc)/4;
cm=(Ac-Bc)/4;

Y=(sin(sig*x/2))^2;

f= (2*cp/Y^3+2*cm/(1-Y)^3)*sig^2*Y*(1-Y) + (-cp/Y^2+cm/(1-Y)^2)*(sig^2/2)*(1-2*Y);

return f;


}
