#include <oxstd.h>
#include <oxfloat.h>
#include <oxprob.h>
#include "rnorm.ox"
#include "dnorm.ox"
#include "rmvnorm.ox"
#include "runif.ox"
#include "BB2.ox"
#include "ivec.ox"
#include "sigmabar.ox"
#include "taubar.ox"
#include "deltaK1.ox"
#include "deltaK.ox"
#include "deltaL.ox"
#include "deltaKmin.ox"
#include "deltaKmax.ox"
#include "igaus.ox"
#include "condminim.ox"
#include "condmaxim.ox"
#include "xminim.ox"
#include "xmaxim.ox"
#include "xminim2.ox"
#include "xmaxim2.ox"
#include "layerass.ox"
#include "xlayer.ox"
#include "xlayerass0.ox"
#include "DI.ox"
#include "phi.ox"
#include "dphi.ox"
#include "ddphi.ox"
#include "nrphi.ox"
#include "phi2.ox"
#include "dphi2.ox"
#include "ddphi2.ox"
#include "nrphi2.ox"
#include "Inti0.ox"
#include "xseed.ox"
#include "deltaKK.ox"
#include "BBx.ox"
#include "XBarkerPAR_new.ox"
#include "leaf_sim_gamma.ox"
#include "DCBF.ox"
#include "rsample.ox"
#include "Intii20.ox"
#include "GBarker2_new.ox"
#include "phi2sig.ox"
#include "Intii20sig.ox"
#include "leaf_sim_sig.ox"
#include "SBarker.ox"



main()
{
  decl time,time2,Yt,z,zpairs,g1,g2,sig,T,dt,dtt,n,nn,M,X1,X2,l,i,s1,s2,s_sig,cp,
       MINI,MAXI,x0,xT,X0,xLrow,xRrow,gama,temp,temp1,acc_g,acc_s,
       ind_stop_temp,cont_phi_temp,cont_w_temp;

time=timer();
ranseed(1);

decl cont_g,cont1,cont_s,ind_stop_x,ind_stop_g,ind_stop_s; cont_g=cont1=cont_s=ind_stop_x=ind_stop_g=ind_stop_s=0;

// ---- data on the ORIGINAL WF scale (0,1); precompute sigma-free Lamperti z=2 asin(sqrt(Y)) ----
Yt=loadmat("cor3m_257.txt");
n=rows(Yt);
z=2*asin(sqrt(Yt));
x0=z[0]; xT=z[n-1];                 // A-endpoint (t=0 and t=T)
zpairs=z[:(n-2)]~z[1:];             // consecutive Lamperti pairs (sigma=1)


// ---- fixed standard-bridge band (derived from data + WF (0.0285,0.9675) at sig_ref=0.12) ----
MINI = -5;  MAXI = 7;

// ---- parameters / priors (SET AS DESIRED) ----

// ---- prior hyperparameters (global; seen by GBarker2 / SBarker) ----
decl PR_g1_a, PR_g1_b;   // theta1 = g1 ~ Gamma(shape,rate)
decl PR_g2_a, PR_g2_b;   // theta2 = g2 ~ Beta
decl PR_s_a,  PR_s_b;    // sigma^2 ~ InvGamma(shape,scale)

g1=0.2;  g2=0.40;  sig=0.1;

PR_g1_a=4.0;  PR_g1_b=20.0;   // theta1 ~ Gamma(shape,rate):      mean 0.2, sd 0.1
PR_g2_a=2.0;  PR_g2_b=3.0;    // theta2 ~ Beta:                   mean 0.4, sd 0.2
PR_s_a =1.2;  PR_s_b =0.005;  // sigma^2 ~ InvGamma(shape,scale): E[sigma]=0.1, sd[sigma]~0.12

T=n-1;  dt=1;  dtt=0.04;  nn=dt/dtt;

M=10000;

s1=0.03*3.5;  s2=0.022*3.5;  s_sig=0.01;   // RW radii (tune s_sig manually)

decl beta1=0.99;

decl NS=T/dtt;
decl ll1=4;  decl NL1=(2).^ll1;
decl ll2=4;  decl NL2=(2).^ll2;

decl scp_prb1=0.005/NL1;
decl scp_prb2=0.005/NL2;

decl gamaa=zeros(M+1,3);
gamaa[0][]=g1~g2~sig;

// ---- initial standard bridges ----
cp=nrphi(g1,g2,sig,M_PI/(2*sig),10^(-10));
[X1,X2]=xseed(zpairs/sig,dt,dtt,MINI,MAXI,cp,g1,g2,sig);

// ================= MCMC =================
for(l=1;l<=M;++l)
{

X0=zpairs/sig;                                  // Lamperti endpoint pairs at current sigma
xLrow=vec((zpairs[][0]/sig * ones(1,nn))');     // per-row left endpoint  (NS x 1)
xRrow=vec((zpairs[][1]/sig * ones(1,nn))');     // per-row right endpoint
cp=nrphi(g1,g2,sig,M_PI/(2*sig),10^(-10));

// X (bridges)
time2=timer();
[X1,X2,temp,temp1]=XBarkerPAR(X0,X1,X2,T,dt,dtt,g1,g2,sig,cp,MINI,MAXI,beta1);
cont1+=temp;
ind_stop_x+=temp1;



time2=timer();
[X2,gama,cp,acc_g,ind_stop_temp,cont_phi_temp,cont_w_temp]=GBarker2(X1,X2,dtt,g1,g2,sig,xLrow,xRrow,cp,s1,s2,MINI,MAXI,x0,xT,NS,NL1,ll1,scp_prb1,PR_g1_a,PR_g1_b,PR_g2_a,PR_g2_b);
cont_g+=acc_g;
ind_stop_g+=ind_stop_temp;
g1=gama[0]; g2=gama[1];



// sigma
time2=timer();
[X1,X2,sig,cp,acc_s,ind_stop_temp]=SBarker(X1,X2,sig,g1,g2,s_sig,zpairs,x0,xT,dt,dtt,T,NS,NL2,ll2,scp_prb2,PR_s_a,PR_s_b);
cont_s+=acc_s;
ind_stop_s+=ind_stop_temp;


gamaa[l][]=g1~g2~sig;
if(fmod(l,10)==0)
{
println(l~g1~g2~sig);
println("acc - x, g , s = ",(cont1/l)~(cont_g/l)~(cont_s/l));
println("Portkey - x, g , s = ",(ind_stop_x/l)~(ind_stop_g/l)~(ind_stop_s/l));
}
}

println("Time = ",timespan(time));
println("acc: X=",cont1/M," th1=",cont_g/M," sigma=",cont_s/M);
savemat("gama_data2.mat",gamaa,1);
savemat("X1_data2.mat",X1);
savemat("X2_data2.mat",X2);
}
