
ddphi2(const g1n, const g2n, const g1o, const g2o, const sig, const x)  // d2/dx2 [phi(new)-phi(old)]

{
  decl al,be,p,q,Ac,Bc,cpn,cmn,cpo,cmo,dcp,dcm,Y,f;

al=g1n*g2n; be=g1n*(1-g2n); p=al-be; q=al+be-sig^2;
Ac=(p^2+q^2)/(8*sig^2)-q/4; Bc=p*q/(4*sig^2)-p/4;
cpn=(Ac+Bc)/4; cmn=(Ac-Bc)/4;

al=g1o*g2o; be=g1o*(1-g2o); p=al-be; q=al+be-sig^2;
Ac=(p^2+q^2)/(8*sig^2)-q/4; Bc=p*q/(4*sig^2)-p/4;
cpo=(Ac+Bc)/4; cmo=(Ac-Bc)/4;

dcp=cpn-cpo; dcm=cmn-cmo;

Y=(sin(sig*x/2))^2;

f= (2*dcp/Y^3+2*dcm/(1-Y)^3)*sig^2*Y*(1-Y) + (-dcp/Y^2+dcm/(1-Y)^2)*(sig^2/2)*(1-2*Y);

return f;


}
