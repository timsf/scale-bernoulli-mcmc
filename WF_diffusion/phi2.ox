
phi2(const g1n, const g2n, const g1o, const g2o, const sig, const x)   // phi(new;sig)-phi(old;sig)

{
  decl al,be,p,q,Ac,Bc,cpn,cmn,c0n,cpo,cmo,c0o,dcp,dcm,dc0,Y,f;

al=g1n*g2n; be=g1n*(1-g2n); p=al-be; q=al+be-sig^2;
Ac=(p^2+q^2)/(8*sig^2)-q/4; Bc=p*q/(4*sig^2)-p/4;
cpn=(Ac+Bc)/4; cmn=(Ac-Bc)/4; c0n=-q^2/(8*sig^2);

al=g1o*g2o; be=g1o*(1-g2o); p=al-be; q=al+be-sig^2;
Ac=(p^2+q^2)/(8*sig^2)-q/4; Bc=p*q/(4*sig^2)-p/4;
cpo=(Ac+Bc)/4; cmo=(Ac-Bc)/4; c0o=-q^2/(8*sig^2);

dcp=cpn-cpo; dcm=cmn-cmo; dc0=c0n-c0o;

Y=(sin(sig*x/2)).^2;

f= dcp./Y + dcm./(1-Y) + dc0;

return f;


}
