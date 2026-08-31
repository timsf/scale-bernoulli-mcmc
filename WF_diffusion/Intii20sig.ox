
// Tight bounds on Delta phi = phi(.;sigp) - phi(.;sig) over the standard-bridge layer bands,
// one interval. Bounds the DIFFERENCE directly so the width -> 0 as sigp -> sig
// (the old version bounded phi(sigp) and phi(sig) independently, giving a width equal to the
// sum of the two full phi-ranges, which does not shrink with the proposal -> very slow sigma step).
//
// Writing phi = cp/Y + cm/(1-Y) + c0 with Y = sin^2(sig*x/2), and along the path
// x = ell/sig + xdot so that sig*x = ell + sig*xdot (ell is sigma-free), the difference splits as
//   Delta phi = [ dcp/Y' + dcm/(1-Y') + dc0 ]            (P1: coefficient differences, ~O(sigp-sig))
//             + (Y - Y')*[ cp/(Y Y') - cm/((1-Y)(1-Y')) ] (P2: Y-mismatch, |Y-Y'| <= |sigp-sig|/2 * |xdot|)
// P1 is bounded by its range over Y' in [YminP,YmaxP] (endpoints + interior critical point);
// P2 by |Y-Y'|_max times the worst-case bracket. Both are O(sigp-sig).
Intii20sig(const xtt, const dtt, const g1, const g2, const sigp, const sig, const cpp, const cp, const zLe, const zRe)
{
  decl n,jvec,ta,tb,al,be,pp;
  decl cpS,cmS,c0S,cpP,cmP,c0P,dcp,dcm,dc0,AcS,BcS,qS,AcP,BcP,qP;
  al=g1*g2; be=g1*(1-g2); pp=al-be;
  qS=al+be-sig^2;  AcS=(pp^2+qS^2)/(8*sig^2)-qS/4;  BcS=pp*qS/(4*sig^2)-pp/4;
  cpS=(AcS+BcS)/4; cmS=(AcS-BcS)/4; c0S=-qS^2/(8*sig^2);
  qP=al+be-sigp^2; AcP=(pp^2+qP^2)/(8*sigp^2)-qP/4; BcP=pp*qP/(4*sigp^2)-pp/4;
  cpP=(AcP+BcP)/4; cmP=(AcP-BcP)/4; c0P=-qP^2/(8*sigp^2);
  dcp=cpP-cpS; dcm=cmP-cmS; dc0=c0P-c0S;

  n=rows(xtt);
  jvec=range(0,n-1)'; ta=jvec*dtt; tb=(jvec+1)*dtt;

  // standard-bridge band (sigma-free) and per-row max|xdot|
  decl stdlo=xtt[][0].*xtt[][2]+(1-xtt[][0]).*xtt[][3];
  decl stdhi=xtt[][0].*xtt[][3]+(1-xtt[][0]).*xtt[][2];
  decl umax=maxr(fabs(stdlo)~fabs(stdhi));

  // Lamperti bands at sig and sigp (affine reconstruction, as in Inti0)
  decl xLs=zLe/sig,  xRs=zRe/sig,  xLp=zLe/sigp, xRp=zRe/sigp;
  decl LaS=xLs+(xRs-xLs)*ta, LbS=xLs+(xRs-xLs)*tb;
  decl xmS=0.5*(LaS+LbS)-0.5*fabs(LaS-LbS)+stdlo, xMS=0.5*(LaS+LbS)+0.5*fabs(LaS-LbS)+stdhi;
  decl LaP=xLp+(xRp-xLp)*ta, LbP=xLp+(xRp-xLp)*tb;
  decl xmP=0.5*(LaP+LbP)-0.5*fabs(LaP-LbP)+stdlo, xMP=0.5*(LaP+LbP)+0.5*fabs(LaP-LbP)+stdhi;

  // Y = sin^2(sig*x/2) is increasing on the Lamperti domain, so band ends give Ymin/Ymax
  decl YaS=(sin(sig*xmS/2)).^2,  YbS=(sin(sig*xMS/2)).^2;
  decl YminS=YaS.*(YaS.<=YbS)+YbS.*(YbS.<YaS), YmaxS=YaS.*(YaS.>=YbS)+YbS.*(YbS.>YaS);
  decl YaP=(sin(sigp*xmP/2)).^2, YbP=(sin(sigp*xMP/2)).^2;
  decl YminP=YaP.*(YaP.<=YbP)+YbP.*(YbP.<YaP), YmaxP=YaP.*(YaP.>=YbP)+YbP.*(YbP.>YaP);

  // P1 range over Y' in [YminP,YmaxP]
  decl P1a=dcp./YminP+dcm./(1-YminP)+dc0;
  decl P1b=dcp./YmaxP+dcm./(1-YmaxP)+dc0;
  decl P1lo=P1a.*(P1a.<=P1b)+P1b.*(P1b.<P1a);
  decl P1hi=P1a.*(P1a.>=P1b)+P1b.*(P1b.>P1a);
  if(dcp*dcm>0){                                   // interior critical point of P1 in Y'
    decl r=sqrt(dcp/dcm), yc=r/(1+r);
    decl P1c=dcp/yc+dcm/(1-yc)+dc0;
    decl inb=(YminP.<yc).*(yc.<YmaxP);
    P1lo=(1-inb).*P1lo + inb.*(P1lo.*(P1lo.<=P1c)+P1c.*(P1c.<P1lo));
    P1hi=(1-inb).*P1hi + inb.*(P1hi.*(P1hi.>=P1c)+P1c.*(P1c.>P1hi));
  }

  // P2 : |Y-Y'| <= |sigp-sig|/2 * umax ; worst-case bracket uses smallest Y and 1-Y
  decl dY=(fabs(sigp-sig)/2).*umax;
  decl P2=dY.*( fabs(cpS)./(YminS.*YminP) + fabs(cmS)./((1-YmaxS).*(1-YmaxP)) );

  decl mm=P1lo-P2, MM_=P1hi+P2;
  return mm ~ MM_ ~ (MM_-mm) ~ (-dtt*mm);          // m ~ M ~ MM ~ I  (same layout as Inti0)
}
