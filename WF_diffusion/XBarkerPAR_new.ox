

XBarkerPAR(const X0, XXX1, const X2, const T, const dt, const dtt, const g1, const g2, const sig, const cp, const MINI, const MAXI, const beta)

{
  
serial decl cont=0;	serial decl cont1=0;

serial decl X1=XXX1;

decl nn=dt/dtt;
decl t=range(dtt,dt-dtt,dtt)';

decl a11=sqrt(dtt)/4;
decl b11=sqrt(dtt)/4;
decl dl=sqrt(dtt)/8;

decl cphi=phi(g1,g2,sig,cp);

decl n=rows(X0);

serial decl X22=zeros(1,14);

serial decl acc=zeros(n,1);

decl i;
parallel for (i=1;i<=n;++i)
{

/////////////////Proposal (STANDARD bridge 0 -> 0)////////////////

decl x=X0[i-1][0];              // Lamperti endpoints (for the affine transform)
decl y=X0[i-1][1];

decl indX1=(((i-1)*nn)|(i*nn-1));
decl X1o=X1[indX1[0]:indX1[1]][];
decl indX2=vecindex(X2[][0],i-1);
indX2=(indX2[0]|indX2[rows(indX2)-1]);
decl X2o=X2[indX2[0]:indX2[1]][];

decl xt=zeros(nn,1);xt[0]=0;
decl tt=0|t;

xt[1:(nn-1)]=xlayer(0,0,dt,0,0-MINI,MAXI-0,t,MINI,MAXI);   // standard 0->0 layered bridge

xt=xt[1:];
decl xx=(0|xt)~(xt|0);
decl X1n=constant(i-1,nn,1)~range(0,nn-1)'~xx;

decl xtt=zeros(nn,5);
decl j,l;
for(j=1;j<=nn;++j)
{
 l=layerass(0,xx[j-1][0],dtt,xx[j-1][1],a11,b11,dl,MINI,MAXI);
 xtt[j-1][]=xlayerass0(0,xx[j-1][0],dtt,xx[j-1][1],l[0][0],l[0][1],l[1][0],l[1][1]);
}

decl X2na=constant(i-1,nn,1)~range(0,nn-1)'~zeros(nn,1)~xtt[][1]~xx[][0]~xtt[][2]~ones(nn,1)~zeros(nn,1)~zeros(nn,6);
decl X2nb=constant(i-1,nn,1)~range(0,nn-1)'~xtt[][1]~constant(dtt,nn,1)~xtt[][2]~xx[][1]~ones(nn,2)~zeros(nn,6);
decl X2n=sortbyc(X2na|X2nb,<1,2>);

X1n=X1n~xtt[][:3];

X1n=X1n~Inti0(xtt,dtt,g1,g2,sig,cp,x,y);

////////////////////////////////////////////

decl s0=sumc(X1o[][11]);
decl s1=sumc(X1n[][11]);
decl p=1/(1+exp(s0-s1));
decl MMMo=X1o[][10];
decl MMMn=X1n[][10];

decl c=0;
decl j1=0; decl j2=0; decl S,C1,np,inde,nnp,cc,indd,up,inddd,x2,phix,Lval;
while(c==0)
{cont+=1;

S = ranbinomial(1,1,1,beta);
if(S == 0)
{
c = 1;
X22|=X2o;	cont1+=1;

}

else
{

 C1=ranbinomial(1,1,1,p);

 if(C1==1)
 {
  np=ranpoisson(nn,1,dtt*MMMn);
  if(sumc(np)==0){c=1;}
  else
  {
   c=1;

   inde=round(vecindex(np));
   np=np[inde];
   nnp=rows(np);

   j=1;
   cc=0;
  while(j<=nnp && cc==0)
  {
   indd=vecindex(X1n[][1],inde[j-1]);
   up=runif(np[j-1],2,zeros(np[j-1],2),(zeros(np[j-1],1)+dtt)~(zeros(np[j-1],1)+MMMn[indd]));
   up=sortbyc(up,0);

    inddd=vecindex(X2n[][1],inde[j-1]);
	[x2,xt]=BBx(X1n[indd][],X2n[inddd][],up[][0],dtt);

    X2n=sortbyc(dropr(X2n,inddd)|x2,<1,2>);		//	renew X2
	Lval=x+(y-x)*(X1n[indd][1]*dtt+up[][0])+xt;
	phix=phi(g1,g2,sig,Lval)-X1n[indd][8];
	if(up[][1]>phix){j+=1;}
	else{cc=1;c=0;}
   
  }

  } //else np=0

  if(c==1)
  {
   X22|=X2n;
   X1[indX1[0]:indX1[1]][]=X1n;
   acc[i-1]=1;
  }

 } //if C1=1

 else //(C1=0)
 {

  np=ranpoisson(nn,1,dtt*MMMo);
  if(sumc(round(np))==0){c=1;}
  else
  {
   c=1;

   inde=vecindex(np);
   np=np[inde];
   nnp=rows(np);

   j=1;
   cc=0;
  while(j<=nnp && cc==0)
  {
   indd=vecindex(X1o[][1],inde[j-1]);
   up=runif(np[j-1],2,zeros(np[j-1],2),(zeros(np[j-1],1)+dtt)~(zeros(np[j-1],1)+MMMo[indd]));
   up=sortbyc(up,0);

    inddd=vecindex(X2o[][1],inde[j-1]);
	[x2,xt]=BBx(X1o[indd][],X2o[inddd][],up[][0],dtt);
		
    X2o=sortbyc(dropr(X2o,inddd)|x2,<1,2>);
	Lval=x+(y-x)*(X1o[indd][1]*dtt+up[][0])+xt;
	phix=phi(g1,g2,sig,Lval)-X1o[indd][8];
	if(up[][1]>phix){j+=1;}
	else{cc=1;c=0;}

  }

  } //else np=0

    if(c==1)
  {
   X22|=X2o;
  }
 
 } //else (C1=0)

} //else (S=0)
 
} //while c=0


}//for i

X22=X22[1:][];

X22=sortbyc(X22,<0,1,2>);

return {X1,X22,sumc(acc)/n,cont1/n};

}
