
phi(const g1, const g2, const sig, const x)

{
  decl al,be,p,q,Ac,Bc,cp,cm,c0,Y,f;

// stationary law is Beta(al/sig^2 , be/sig^2);  p=al-be, q=al+be-sig^2
al=g1*g2;
be=g1*(1-g2);
p=al-be;
q=al+be-sig^2;

Ac=(p^2+q^2)/(8*sig^2) - q/4;
Bc=p*q/(4*sig^2) - p/4;

cp=(Ac+Bc)/4;                 // coeff of 1/Y
cm=(Ac-Bc)/4;                 // coeff of 1/(1-Y)
c0=-q^2/(8*sig^2);            // constant

Y=(sin(sig*x/2)).^2;          // = original WF value along the (sigma) Lamperti path

f= cp./Y + cm./(1-Y) + c0;

return f;


}
