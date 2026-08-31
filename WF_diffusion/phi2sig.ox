
// sigma-difference at a standard-bridge point: phi(L^sigp+xdot;sigp) - phi(L^sig+xdot;sig)
// ell = sigma-free linear interp of z at absolute time tabs;  Lamperti = ell/sigma + xdot
phi2sig(const g1, const g2, const sigp, const sig, const zLe, const zRe, const tabs, const xdot)
{
  decl ell;
  ell = zLe + (zRe-zLe).*tabs;
  return phi(g1,g2,sigp, ell/sigp + xdot) - phi(g1,g2,sig, ell/sig + xdot);
}
