
Inti0(const xtt, const dtt, const g1, const g2, const sig, const cp, const xL, const xR)
{
  decl ind0,ind1,cphi,m,M,xm,xM,MM,I,al,be,p,q,Ac,Bc,kp,km,rtn,n,i,
       jvec,ta,tb,La,Lb,Lmin,Lmax,stdlo,stdhi;
al=g1*g2; be=g1*(1-g2); p=al-be; q=al+be-sig^2;
Ac=(p^2+q^2)/(8*sig^2)-q/4; Bc=p*q/(4*sig^2)-p/4;
kp=(Ac+Bc)/4;  km=(Ac-Bc)/4;
ind0=0; if( km>0 ){ind0=1;}
ind1=0; if( kp>0 ){ind1=1;}
n=rows(xtt);
jvec=range(0,n-1)';
ta=jvec*dtt;  tb=(jvec+1)*dtt;
La=xL+(xR-xL)*ta;  Lb=xL+(xR-xL)*tb;
Lmin=0.5*(La+Lb)-0.5*fabs(La-Lb);  Lmax=0.5*(La+Lb)+0.5*fabs(La-Lb);
stdlo=xtt[][0].*xtt[][2]+(1-xtt[][0]).*xtt[][3];
stdhi=xtt[][0].*xtt[][3]+(1-xtt[][0]).*xtt[][2];
xm=Lmin+stdlo;  xM=Lmax+stdhi;
m=M=MM=I=zeros(n,1);
cphi=phi(g1,g2,sig,cp);
if(ind0==1 && ind1==0){
 m=phi(g1,g2,sig,xm); M=phi(g1,g2,sig,xM); I=(-dtt*m); MM=M-m; }
else if(ind0==0 && ind1==1){
 m=phi(g1,g2,sig,xM); M=phi(g1,g2,sig,xm); I=(-dtt*m); MM=M-m; }
else if(ind0==1 && ind1==1){
 m=M=zeros(n,1);
 for(i=1;i<=n;++i){
  if(xm[i-1]>cp){ m[i-1]=phi(g1,g2,sig,xm[i-1]); M[i-1]=phi(g1,g2,sig,xM[i-1]); }
  else if(xM[i-1]<cp){ m[i-1]=phi(g1,g2,sig,xM[i-1]); M[i-1]=phi(g1,g2,sig,xm[i-1]); }
  else{ m[i-1]=cphi; M[i-1]=max(phi(g1,g2,sig,xM[i-1]),phi(g1,g2,sig,xm[i-1])); } }
 I=(-dtt*m); MM=M-m; }
else{
 m=M=zeros(n,1);
 for(i=1;i<=n;++i){
  if(xm[i-1]>cp){ m[i-1]=phi(g1,g2,sig,xM[i-1]); M[i-1]=phi(g1,g2,sig,xm[i-1]); }
  else if(xM[i-1]<cp){ m[i-1]=phi(g1,g2,sig,xm[i-1]); M[i-1]=phi(g1,g2,sig,xM[i-1]); }
  else{ m[i-1]=min(phi(g1,g2,sig,xM[i-1]),phi(g1,g2,sig,xm[i-1])); M[i-1]=cphi; } }
 I=(-dtt*m); MM=M-m; }
rtn=m~M~MM~I;
return rtn;
}
