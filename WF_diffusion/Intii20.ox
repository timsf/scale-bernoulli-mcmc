
Intii20(const X1, const dtt, const g1n, const g2n, const g1o, const g2o, const sig, const cp, const cphi, const xLrow, const xRrow)
{
  decl n,ind0,ind1,I,m,M,MM,al,be,p,q,Ac,Bc,kpn,kmn,kpo,kmo,dkp,dkm,xm,xM,rtn,i,
       ta,tb,La,Lb,Lmin,Lmax,stdlo,stdhi;
n=rows(X1);
al=g1n*g2n; be=g1n*(1-g2n); p=al-be; q=al+be-sig^2;
Ac=(p^2+q^2)/(8*sig^2)-q/4; Bc=p*q/(4*sig^2)-p/4; kpn=(Ac+Bc)/4; kmn=(Ac-Bc)/4;
al=g1o*g2o; be=g1o*(1-g2o); p=al-be; q=al+be-sig^2;
Ac=(p^2+q^2)/(8*sig^2)-q/4; Bc=p*q/(4*sig^2)-p/4; kpo=(Ac+Bc)/4; kmo=(Ac-Bc)/4;
dkp=kpn-kpo; dkm=kmn-kmo;
ind0=0; if( dkm>0 ){ind0=1;}
ind1=0; if( dkp>0 ){ind1=1;}
ta=X1[][1]*dtt;  tb=(X1[][1]+1)*dtt;
La=xLrow+(xRrow-xLrow).*ta;  Lb=xLrow+(xRrow-xLrow).*tb;
Lmin=0.5*(La+Lb)-0.5*fabs(La-Lb);  Lmax=0.5*(La+Lb)+0.5*fabs(La-Lb);
stdlo=X1[][4].*X1[][6]+(1-X1[][4]).*X1[][7];
stdhi=(1-X1[][4]).*X1[][6]+X1[][4].*X1[][7];
xm=Lmin+stdlo;  xM=Lmax+stdhi;
m=M=zeros(n,1);
if(ind0==1 && ind1==0){
 m=phi2(g1n,g2n,g1o,g2o,sig,xm); M=phi2(g1n,g2n,g1o,g2o,sig,xM); I=(-dtt*m); MM=M-m; }
else if(ind0==0 && ind1==1){
 m=phi2(g1n,g2n,g1o,g2o,sig,xM); M=phi2(g1n,g2n,g1o,g2o,sig,xm); I=(-dtt*m); MM=M-m; }
else if(ind0==1 && ind1==1){
 m=M=zeros(n,1);
 for(i=1;i<=n;++i){
  if(xm[i-1]>cp){ m[i-1]=phi2(g1n,g2n,g1o,g2o,sig,xm[i-1]); M[i-1]=phi2(g1n,g2n,g1o,g2o,sig,xM[i-1]); }
  else if(xM[i-1]<cp){ m[i-1]=phi2(g1n,g2n,g1o,g2o,sig,xM[i-1]); M[i-1]=phi2(g1n,g2n,g1o,g2o,sig,xm[i-1]); }
  else{ m[i-1]=cphi; M[i-1]=max(phi2(g1n,g2n,g1o,g2o,sig,xM[i-1]),phi2(g1n,g2n,g1o,g2o,sig,xm[i-1])); } }
 I=(-dtt*m); MM=M-m; }
else{
 m=M=zeros(n,1);
 for(i=1;i<=n;++i){
  if(xm[i-1]>cp){ m[i-1]=phi2(g1n,g2n,g1o,g2o,sig,xM[i-1]); M[i-1]=phi2(g1n,g2n,g1o,g2o,sig,xm[i-1]); }
  else if(xM[i-1]<cp){ m[i-1]=phi2(g1n,g2n,g1o,g2o,sig,xm[i-1]); M[i-1]=phi2(g1n,g2n,g1o,g2o,sig,xM[i-1]); }
  else{ m[i-1]=min(phi2(g1n,g2n,g1o,g2o,sig,xM[i-1]),phi2(g1n,g2n,g1o,g2o,sig,xm[i-1])); M[i-1]=cphi; } }
 I=(-dtt*m); MM=M-m; }
rtn=m~M~MM~I;
return rtn;
}
