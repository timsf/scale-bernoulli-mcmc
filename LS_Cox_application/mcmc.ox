#include <oxstd.h>
#include <oxfloat.h>
#include <oxprob.h>

#include "dmatrix.ox"
#include "dmatrix2.ox"
#include "rmvnorm.ox"
#include "runif.ox"
#include "rnorm.ox"
#include "cov_func.ox"
#include "GK.ox"
#include "locc.ox"
#include "Sampling_SPG_cond3.ox"
#include "Sampling_SPG_cond2.ox"
#include "RegiaoK.ox"
#include "RK.ox"
#include "Metropolis_beta_pCN.ox"
#include "Metropolis_c.ox"
#include "Metropolis_N.ox"

#include "RegiaoK3.ox"
#include "RegiaoK4.ox"
#include "asng_vl_mt.ox"
#include "mtxfind.ox"
#include "DCBF.ox"

#include "leaf_sim_lamb.ox"
#include "DCBF_lambda.ox"


main()
{

ranseed(555);  // setting seed

decl Y, y, MuS, L, U, LU, delta, lambda0, Lest0, lambdaCand, Lest_cand, w, np, Xn, Yn, Zn, N0; 
decl N0_Lest_cand, Yk0, Nk0, BetaY0, BetaN0, M_C, M_Y, par, L1, L2,L3;
decl rc, cov_prop, i, dist, ind, i0, SS, S12, S22, S22_inv,j,BetaSa,K,aux;




Y= loadmat("alseis_blackiana.txt");			// reading data file
y=rows(Y);	 println(y);


///////////// Setting some parameter values //////////////

K=3;			// value of K

L=82.4; U=41.2;    // L = horizontal length  ,  U = vertical length

MuS= L*U;
LU=min(L,U);

delta=5;		// delta parameter from the pseudo-marginal estimator
decl deltab=1;		// delta parameter from the DCBF estimator
par=2;			// range parameter tau2 of the covariance function of the GP



cov_prop=(<0.08,0.14,0.25>);

///////////// Initial value and prior hyperparameters for parameter lambda //////////////


lambda0=<1.3 , 3, 6>;

//lambda0=zeros(1,K)+(y/MuS);


decl lambdaAtual=lambda0;

Lest0 = delta*max(lambda0)-min(lambda0);


////////////// Initial value for c //////////////

decl cAtual;
cAtual= <-0.2 , 1>;



/////////// Definition of the partition to update N //////////////

w=Lest0;

decl LL=25;			// S is split into LL*(LL*U/L) parts
decl d=L/LL;
LL=LL~(ceil(U/d));


decl LL2=LL[0]*LL[1];

decl limss1=range(0,L,L/LL[0])';
decl limss2=range(0,U,U/LL[1])';

decl lims=zeros(1,4);

for(i=1;i<=LL[0];++i)
{
  lims|=(zeros(LL[1],2)+(limss1[i-1][0]~limss1[i][0]))~(limss2[:(LL[1]-1)][0]~limss2[1:][0]);
}

lims=lims[1:][];

decl npp=ranpoisson(LL2,1,w*MuS/LL2);

N0=zeros(1,3);

for(i=1;i<=LL2;++i)
{
N0|= runif(npp[i-1],1,lims[i-1][0],lims[i-1][1]) ~ runif(npp[i-1],1,lims[i-1][2],lims[i-1][3]) ~ runif(npp[i-1],1,0,w);
}

N0=N0[1:][];



////////////// NNGP specifications /////////////////

decl n1=75;		 	// number of values in the horizontal coordinate - if L=U, the total number of points is nxn
decl n2=floor(n1*(U/L));		// number of values in the vertical coordinate - if L=U, the total number of points is nxn

decl m=16;		// number of nearest neighbors used in the NNGP prior

decl s2S=((0.012/4)^2);     // proposal variance for the RW proposal for the GP at S

decl S=zeros(1,2);

decl SiL=range(L/(2*n1),L-L/(2*n1),L/n1)';
decl SiU=range(U/(2*n2),U-U/(2*n2),U/n2)';

decl nL=rows(SiL);
decl nU=rows(SiU);


for(i=1;i<=nL;++i)
{
 S=S|(SiL[i-1]~SiU);
}
S=S[1:][];				   // S is the reference set of locations from the NNGP prior

SiL=0|SiL|L;
SiU=0|SiU|U;



decl Cmu,CCov,Iind;

Cmu=Iind= new array[nL*nU];
Cmu[0]=0;
CCov=zeros(nL*nU,1);
CCov[0]=1;
Iind[0]=0;



for(i=1;i<=(nL*nU-1);++i)
{

 dist=S[i][]|S[:(i-1)][];
 dist=dmatrix2(dist);
 ind=sortcindex(dist);
 
 i0=min(m,i);
 Iind[i]= ind[:(i0-1)];			   // sets of the m closest neighbors for each of the n^2 locations
 
 dist=S[i][]|S[Iind[i]][];
 dist=dmatrix(dist);

 SS=COV_GP(dist, par);
 SS=(SS+SS')/2;
 S12=  SS[0][1:];	  
 S22=  SS[1:][1:];
 S22_inv=invert(S22);
 SS=S12*S22_inv;

 Cmu[i]=SS;						// the scale factor of the conditional mean for each of the n^2 locations
 CCov[i]=1-SS*S12';				// conditional variance for each of the n^2 locations

}



 /////////// Initial value of the GP at locations Y, N and S ////////
 

decl BetaS;

BetaS=zeros(rows(S),1);

BetaS[0]=0;

for(i=1;i<=(rows(S)-1);++i)
{
BetaS[i]=Cmu[i]*BetaS[Iind[i]] + sqrt(CCov[i])*rann(1,1);
}


BetaS= loadmat("BetaS0.mat");

BetaY0=Sampling_SPG_cond2(Y, L, U, nL, nU, SiL, SiU, S, BetaS, par, m);
BetaN0=Sampling_SPG_cond2(N0[][:1], L, U, nL, nU, SiL, SiU, S, BetaS, par, m);

M_C = N0~BetaN0;
M_Y = Y~BetaY0;



decl nnn=nU*nL;


decl temp=RegiaoK(K, cAtual, BetaS);
decl n=rows(BetaS);
decl Yc=RegiaoK(K, cAtual, BetaY0);

decl logVer =  sumr(-(temp/n).*lambdaAtual) + sumr(Yc.*log(lambdaAtual));



///////////////////////////////////////////////
///////////////////////////////////////////////
///////////////////////////////////////////////


///////////// Parameters for the DCBF of lambdas //////////////


decl LLk=<6;6;6>;	 // number of levels
decl NLk=2 .^LLk;	 // number of leaves
decl NSk=4 .^LLk;	 // number of squares

println(LLk,NLk,NSk);



decl n_temp=sqrt(NSk[0]);
limss1=range(0,L,L/n_temp)';
limss2=range(0,U,U/n_temp)';

decl lims1=zeros(1,4);

for(i=1;i<=n_temp;++i)
{
  lims1|=(zeros(n_temp,2)+(limss1[i-1][0]~limss1[i][0]))~(limss2[:(n_temp-1)][0]~limss2[1:][0]);
}

lims1=lims1[1:][];

//////////////

n_temp=sqrt(NSk[1]);
limss1=range(0,L,L/n_temp)';
limss2=range(0,U,U/n_temp)';

decl lims2=zeros(1,4);

for(i=1;i<=n_temp;++i)
{
  lims2|=(zeros(n_temp,2)+(limss1[i-1][0]~limss1[i][0]))~(limss2[:(n_temp-1)][0]~limss2[1:][0]);
}

lims2=lims2[1:][];


//////////////

n_temp=sqrt(NSk[2]);
limss1=range(0,L,L/n_temp)';
limss2=range(0,U,U/n_temp)';

decl lims3=zeros(1,4);

for(i=1;i<=n_temp;++i)
{
  lims3|=(zeros(n_temp,2)+(limss1[i-1][0]~limss1[i][0]))~(limss2[:(n_temp-1)][0]~limss2[1:][0]);
}

lims3=lims3[1:][];


//////////////



decl scp_prbk=<1/(4^3*1000),1/(4^3*1000),1/(4^3*1000)>;		// scape probability for the 2-coin



decl Yk=RegiaoK3(K, cAtual, M_Y[][2]);

decl Yk1=Y[vecindex(Yk,0)][];
decl Yk2=Y[vecindex(Yk,1)][];
decl Yk3=Y[vecindex(Yk,2)][];


n_temp=NSk[0];
decl Y_ind_sqr1=zeros(n_temp,1);

decl ntem=rows(Yk1);
parallel for(i=1;i<=ntem;++i)
{
decl limstemp=lims1;
decl Yk1temp=Yk1[i-1][];
decl ind1=vecindex( limstemp[][0].<= Yk1temp[][0] );
decl ind2=vecindex( limstemp[][1].>  Yk1temp[][0] );
decl ind3=vecindex( limstemp[][2].<= Yk1temp[][1] );
decl ind4=vecindex( limstemp[][3].>  Yk1temp[][1] );

decl ind_temp=intersection(ind1,ind2)';
ind_temp=intersection(ind_temp,ind3)';
ind_temp=intersection(ind_temp,ind4)';

Y_ind_sqr1[ind_temp]+=1;

}

n_temp=NSk[1];
decl Y_ind_sqr2=zeros(n_temp,1);

ntem=rows(Yk2);
parallel for(i=1;i<=ntem;++i)
{
decl limstemp=lims2;
decl Yk2temp=Yk2[i-1][];
decl ind1=vecindex( limstemp[][0].<= Yk2temp[][0] );
decl ind2=vecindex( limstemp[][1].>  Yk2temp[][0] );
decl ind3=vecindex( limstemp[][2].<= Yk2temp[][1] );
decl ind4=vecindex( limstemp[][3].>  Yk2temp[][1] );

decl ind_temp=intersection(ind1,ind2)';
ind_temp=intersection(ind_temp,ind3)';
ind_temp=intersection(ind_temp,ind4)';

Y_ind_sqr2[ind_temp]+=1;

}

n_temp=NSk[2];
decl Y_ind_sqr3=zeros(n_temp,1);

ntem=rows(Yk3);
parallel for(i=1;i<=ntem;++i)
{
decl limstemp=lims3;
decl Yk3temp=Yk3[i-1][];
decl ind1=vecindex( limstemp[][0].<= Yk3temp[][0] );
decl ind2=vecindex( limstemp[][1].>  Yk3temp[][0] );
decl ind3=vecindex( limstemp[][2].<= Yk3temp[][1] );
decl ind4=vecindex( limstemp[][3].>  Yk3temp[][1] );

decl ind_temp=intersection(ind1,ind2)';
ind_temp=intersection(ind_temp,ind3)';
ind_temp=intersection(ind_temp,ind4)';

Y_ind_sqr3[ind_temp]+=1;

}


////////////// Leaf allocation and flipped indicator for lambdas ////////////////////

decl u1,u2,beta_sqr;

decl ind_flp1=zeros(NSk[0],1);

ntem=NSk[0];

parallel for (i=1;i<=ntem;++i)
{

u1=((lims1[i-1][0]+lims1[i-1][1])/2);
u2=((lims1[i-1][2]+lims1[i-1][3])/2);

beta_sqr= Sampling_SPG_cond2(u1~u2, L, U, nL, nU, SiL, SiU, S, BetaS, par, m);

decl rgn=RegiaoK(K, cAtual, beta_sqr);
if(rgn[0]==1){ind_flp1[i-1]=1;}


}  // for i

decl ind_std1=1-ind_flp1;


decl iind_flp=mtxfind(ind_flp1,1);
decl n_flp=rows(iind_flp);

decl iind_std=mtxfind(ind_std1,1);
decl n_std=rows(iind_std);

decl n_flp_lf1=ceil((n_flp/(n_flp+n_std))*NLk[0]);
decl n_std_lf1=NLk[0]-n_flp_lf1;

temp=(ranshuffle(NLk[0],range(0,NLk[0]-1)'))';
decl ind_flp_lf1=zeros(NLk[0],1);
if(n_flp_lf1>0){ind_flp_lf1[temp[0:(n_flp_lf1-1)]]=1;}

decl index_flp_lf=vecindex(ind_flp_lf1);
decl index_std_lf=vecindex(1-ind_flp_lf1);

decl ind_lfs_flp=reshape(index_flp_lf,n_flp,1);
decl ind_lfs_std=reshape(index_std_lf,n_std,1);
decl temp_flp=fmod(n_flp,n_flp_lf1);
decl temp_std=fmod(n_std,n_std_lf1);

decl ind_lfs1=zeros(NSk[0],1);
			 
if(n_flp>0)
{
ind_lfs_flp=ranshuffle(n_flp,ind_lfs_flp)';
ind_lfs1=asng_vl_mt(ind_lfs1, ind_lfs_flp, iind_flp);
}

if(n_std>0)
{
ind_lfs_std=ranshuffle(n_std,ind_lfs_std)';
ind_lfs1=asng_vl_mt(ind_lfs1, ind_lfs_std, iind_std);
}



//////////////

decl ind_flp2=zeros(NSk[1],1);

ntem=NSk[1];

parallel for (i=1;i<=ntem;++i)
{

u1=((lims2[i-1][0]+lims2[i-1][1])/2);
u2=((lims2[i-1][2]+lims2[i-1][3])/2);

beta_sqr= Sampling_SPG_cond2(u1~u2, L, U, nL, nU, SiL, SiU, S, BetaS, par, m);

decl rgn=RegiaoK(K, cAtual, beta_sqr);
if(rgn[1]==1){ind_flp2[i-1]=1;}


}  // for i

decl ind_std2=1-ind_flp2;


iind_flp=mtxfind(ind_flp2,1);
n_flp=rows(iind_flp);

iind_std=mtxfind(ind_std2,1);
n_std=rows(iind_std);

decl n_flp_lf2=ceil((n_flp/(n_flp+n_std))*NLk[1]);
decl n_std_lf2=NLk[1]-n_flp_lf2;

temp=(ranshuffle(NLk[1],range(0,NLk[1]-1)'))';
decl ind_flp_lf2=zeros(NLk[1],1);
if(n_flp_lf2>0){ind_flp_lf2[temp[0:(n_flp_lf2-1)]]=1;}

index_flp_lf=vecindex(ind_flp_lf2);
index_std_lf=vecindex(1-ind_flp_lf2);

ind_lfs_flp=reshape(index_flp_lf,n_flp,1);
ind_lfs_std=reshape(index_std_lf,n_std,1);
temp_flp=fmod(n_flp,n_flp_lf2);
temp_std=fmod(n_std,n_std_lf2);

decl ind_lfs2=zeros(NSk[1],1);
			 
if(n_flp>0)
{
ind_lfs_flp=ranshuffle(n_flp,ind_lfs_flp)';
ind_lfs2=asng_vl_mt(ind_lfs2, ind_lfs_flp, iind_flp);
}

if(n_std>0)
{
ind_lfs_std=ranshuffle(n_std,ind_lfs_std)';
ind_lfs2=asng_vl_mt(ind_lfs2, ind_lfs_std, iind_std);
}


//////////////

decl ind_flp3=zeros(NSk[2],1);

ntem=NSk[2];

parallel for (i=1;i<=ntem;++i)
{

u1=((lims3[i-1][0]+lims3[i-1][1])/2);
u2=((lims3[i-1][2]+lims3[i-1][3])/2);

beta_sqr= Sampling_SPG_cond2(u1~u2, L, U, nL, nU, SiL, SiU, S, BetaS, par, m);

decl rgn=RegiaoK(K, cAtual, beta_sqr);
if(rgn[2]==1){ind_flp3[i-1]=1;}


}  // for i

decl ind_std3=1-ind_flp3;


iind_flp=mtxfind(ind_flp3,1);
n_flp=rows(iind_flp);

iind_std=mtxfind(ind_std3,1);
n_std=rows(iind_std);

decl n_flp_lf3=ceil((n_flp/(n_flp+n_std))*NLk[2]);
decl n_std_lf3=NLk[2]-n_flp_lf3;

temp=(ranshuffle(NLk[2],range(0,NLk[2]-1)'))';
decl ind_flp_lf3=zeros(NLk[2],1);
if(n_flp_lf3>0){ind_flp_lf3[temp[0:(n_flp_lf3-1)]]=1;}

index_flp_lf=vecindex(ind_flp_lf3);
index_std_lf=vecindex(1-ind_flp_lf3);

ind_lfs_flp=reshape(index_flp_lf,n_flp,1);
ind_lfs_std=reshape(index_std_lf,n_std,1);
temp_flp=fmod(n_flp,n_flp_lf3);
temp_std=fmod(n_std,n_std_lf3);

decl ind_lfs3=zeros(NSk[2],1);
			 
if(n_flp>0)
{
ind_lfs_flp=ranshuffle(n_flp,ind_lfs_flp)';
ind_lfs3=asng_vl_mt(ind_lfs3, ind_lfs_flp, iind_flp);
}

if(n_std>0)
{
ind_lfs_std=ranshuffle(n_std,ind_lfs_std)';
ind_lfs3=asng_vl_mt(ind_lfs3, ind_lfs_std, iind_std);
}



////////////// Parameters for the MCMC chain ////////////////////


decl NumSim=11000;  	// total number of iterations

decl burn=1000;	   // number of iterations discarded for burn-in

decl burn1=100;

//////////////////////////////////////////////////////////////////////////////

decl c, NkAtual, YkAtual;
decl  sS1, sS2, slambda, lambda, Ibeta, Ilambda1, Ilambda2, Ilambda3, Ilambda4, S1, S2, Ick, Inn, Il1, Il2, Il3, Il4;
decl In,Ib,Il,Ic;
decl smuk,sYk,sD,muk,D;

sS1=zeros(rows(S),1);
sS2=zeros(rows(S),K);
smuk=zeros(1,K);
sYk=zeros(1,K);
sD=zeros(1,1);
slambda=zeros(1,K);
Il1=Il2=Il3=Il4=Ib=0;
In=zeros(1,LL2);

decl Ykk, IntIF, loglik, loglikPM, mmuk;

loglik=zeros(NumSim+1,1);
loglik[0]=logVer;

lambda = zeros(NumSim+1,K);
lambda[0][]=lambdaAtual;
c = zeros(NumSim+1,K-1);
c[0][] = cAtual;
IntIF  = zeros(NumSim,1);
mmuk = zeros(NumSim,K);

decl cont_lf1=zeros(NumSim,1);
decl cont_lf2=zeros(NumSim,1);
decl cont_lf3=zeros(NumSim,1);
decl cont_lf_temp,ind_stop_temp,ind_stop1,ind_stop2,ind_stop3,ind_stop4,ind_stop1_temp,ind_stop2_temp,ind_stop3_temp,ind_stop4_temp;

ind_stop1=ind_stop2=ind_stop3=ind_stop4=0;
ind_stop1_temp=ind_stop2_temp=ind_stop3_temp=ind_stop4_temp=0;


decl M_C_cand, w_cand, npp_cand, cont_phi_temp, cont_w_temp;

decl cont_w1=zeros(NumSim,1);
decl cont_phi1=zeros(NumSim,1);
decl cont_w2=zeros(NumSim,1);
decl cont_phi2=zeros(NumSim,1);
decl cont_w3=zeros(NumSim,1);
decl cont_phi3=zeros(NumSim,1);



decl time=timer();
decl time2=0;
decl time3=0;

decl time_temp,s,ss,lv;



for(i=1;i<=NumSim;++i)	  	  // begin MCMC iterations
{




lambdaCand = lambdaAtual + runif(1,K,-cov_prop,cov_prop);


time_temp=timer();
[Ilambda1,lambdaAtual,ind_stop_temp,cont_phi_temp,cont_w_temp] = DCBF_lambda(K, 0, deltab, lambdaAtual, lambdaCand[0]~lambdaAtual[1:2], cAtual, M_Y, L, U, nL, nU, SiL, SiU, S, BetaS, par, m, Iind, Cmu, CCov, s2S,
														lims1, scp_prbk[0], LLk[0], NLk[0], NSk[0], ind_lfs1, n_flp_lf1, n_std_lf1, Y_ind_sqr1, ind_flp_lf1);

Il1+=Ilambda1;
cont_w1[i-1]=cont_w_temp;
cont_phi1[i-1]=cont_phi_temp;
ind_stop1+=ind_stop_temp;
ind_stop1_temp+=(1-ind_stop_temp);



[Ilambda2,lambdaAtual,ind_stop_temp,cont_phi_temp,cont_w_temp] = DCBF_lambda(K, 1, deltab, lambdaAtual, lambdaAtual[0]~lambdaCand[1]~lambdaAtual[2], cAtual, M_Y, L, U, nL, nU, SiL, SiU, S, BetaS, par, m, Iind, Cmu, CCov, s2S,
														lims2, scp_prbk[1], LLk[1], NLk[1], NSk[1], ind_lfs2, n_flp_lf2, n_std_lf2, Y_ind_sqr2, ind_flp_lf2);


Il2+=Ilambda2;
cont_w2[i-1]=cont_w_temp;
cont_phi2[i-1]=cont_phi_temp;
ind_stop2+=ind_stop_temp;
ind_stop2_temp+=(1-ind_stop_temp);





[Ilambda3,lambdaAtual,ind_stop_temp,cont_phi_temp,cont_w_temp] = DCBF_lambda(K, 2, deltab, lambdaAtual, lambdaAtual[0:1]~lambdaCand[2], cAtual, M_Y, L, U, nL, nU, SiL, SiU, S, BetaS, par, m, Iind, Cmu, CCov, s2S,
														lims3, scp_prbk[2], LLk[2], NLk[2], NSk[2], ind_lfs3, n_flp_lf3, n_std_lf3, Y_ind_sqr3, ind_flp_lf3);


Il3+=Ilambda3;
cont_w3[i-1]=cont_w_temp;
cont_phi3[i-1]=cont_phi_temp;
ind_stop3+=ind_stop_temp;
ind_stop3_temp+=(1-ind_stop_temp);


														time_temp=timespan(time_temp);
sscan(time_temp, "%d.%d", &s, &ss);  // parse hours, minutes, seconds
time2 +=  s + ss;

if(i==100){println(time2);}


lambda[i][] = lambdaAtual;



/////////// update N ///////////

time_temp=timer();

[M_C,npp,Inn] = Metropolis_N(K, LL2, npp, lims, delta, lambdaAtual, cAtual, M_C, L, U, nL, nU, SiL, SiU, S, BetaS, par, m, MuS);
In+=Inn';


/////////// update beta ///////////


[M_Y,M_C, BetaS, Ibeta, logVer] = Metropolis_beta_pCN(K, delta, lambdaAtual, cAtual, M_C, M_Y, w, L, U, nL, nU, SiL, SiU, S, BetaS, par, m, Iind, Cmu, CCov, s2S);

Ib+=Ibeta;
loglik[i]=logVer;


time_temp=timespan(time_temp);
sscan(time_temp, "%d.%d", &s, &ss);  // parse hours, minutes, seconds
time3 +=  s + ss;


//////////////////////////////////////////////
//////////////////////////////////////////////
//////////////////////////////////////////////

if(Ibeta==1)
{

decl Yk=RegiaoK3(K, cAtual, M_Y[][2]);

decl Yk1=Y[vecindex(Yk,0)][];
decl Yk2=Y[vecindex(Yk,1)][];
decl Yk3=Y[vecindex(Yk,2)][];


ind_flp1=zeros(NSk[0],1);

ntem=NSk[0];

parallel for (j=1;j<=ntem;++j)
{

u1=((lims1[j-1][0]+lims1[j-1][1])/2);
u2=((lims1[j-1][2]+lims1[j-1][3])/2);

beta_sqr= Sampling_SPG_cond2(u1~u2, L, U, nL, nU, SiL, SiU, S, BetaS, par, m);

decl rgn=RegiaoK(K, cAtual, beta_sqr);
if(rgn[0]==1){ind_flp1[j-1]=1;}


}  // for i

ind_std1=1-ind_flp1;


iind_flp=mtxfind(ind_flp1,1);
n_flp=rows(iind_flp);

iind_std=mtxfind(ind_std1,1);
n_std=rows(iind_std);

n_flp_lf1=ceil((n_flp/(n_flp+n_std))*NLk[0]);
n_std_lf1=NLk[0]-n_flp_lf1;

temp=(ranshuffle(NLk[0],range(0,NLk[0]-1)'))';
ind_flp_lf1=zeros(NLk[0],1);
if(n_flp_lf1>0){ind_flp_lf1[temp[0:(n_flp_lf1-1)]]=1;}

index_flp_lf=vecindex(ind_flp_lf1);
index_std_lf=vecindex(1-ind_flp_lf1);

ind_lfs_flp=reshape(index_flp_lf,n_flp,1);
ind_lfs_std=reshape(index_std_lf,n_std,1);
temp_flp=fmod(n_flp,n_flp_lf1);
temp_std=fmod(n_std,n_std_lf1);

ind_lfs1=zeros(NSk[0],1);
		 
if(n_flp>0)
{
ind_lfs_flp=ranshuffle(n_flp,ind_lfs_flp)';
ind_lfs1=asng_vl_mt(ind_lfs1, ind_lfs_flp, iind_flp);
}

if(n_std>0)
{
ind_lfs_std=ranshuffle(n_std,ind_lfs_std)';
ind_lfs1=asng_vl_mt(ind_lfs1, ind_lfs_std, iind_std);
}


//////////////

ind_flp2=zeros(NSk[1],1);

ntem=NSk[1];

parallel for (j=1;j<=ntem;++j)
{

u1=((lims2[j-1][0]+lims2[j-1][1])/2);
u2=((lims2[j-1][2]+lims2[j-1][3])/2);

beta_sqr= Sampling_SPG_cond2(u1~u2, L, U, nL, nU, SiL, SiU, S, BetaS, par, m);

decl rgn=RegiaoK(K, cAtual, beta_sqr);
if(rgn[1]==1){ind_flp2[j-1]=1;}


}  // for i

ind_std2=1-ind_flp2;


iind_flp=mtxfind(ind_flp2,1);
n_flp=rows(iind_flp);

iind_std=mtxfind(ind_std2,1);
n_std=rows(iind_std);

n_flp_lf2=ceil((n_flp/(n_flp+n_std))*NLk[1]);
n_std_lf2=NLk[1]-n_flp_lf2;

temp=(ranshuffle(NLk[1],range(0,NLk[1]-1)'))';
ind_flp_lf2=zeros(NLk[1],1);
if(n_flp_lf2>0){ind_flp_lf2[temp[0:(n_flp_lf2-1)]]=1;}

index_flp_lf=vecindex(ind_flp_lf2);
index_std_lf=vecindex(1-ind_flp_lf2);

ind_lfs_flp=reshape(index_flp_lf,n_flp,1);
ind_lfs_std=reshape(index_std_lf,n_std,1);
temp_flp=fmod(n_flp,n_flp_lf2);
temp_std=fmod(n_std,n_std_lf2);

ind_lfs2=zeros(NSk[1],1);
			 
if(n_flp>0)
{
ind_lfs_flp=ranshuffle(n_flp,ind_lfs_flp)';
ind_lfs2=asng_vl_mt(ind_lfs2, ind_lfs_flp, iind_flp);
}

if(n_std>0)
{
ind_lfs_std=ranshuffle(n_std,ind_lfs_std)';
ind_lfs2=asng_vl_mt(ind_lfs2, ind_lfs_std, iind_std);
}



//////////////

ind_flp3=zeros(NSk[2],1);

ntem=NSk[2];

parallel for (j=1;j<=ntem;++j)
{

u1=((lims3[j-1][0]+lims3[j-1][1])/2);
u2=((lims3[j-1][2]+lims3[j-1][3])/2);

beta_sqr= Sampling_SPG_cond2(u1~u2, L, U, nL, nU, SiL, SiU, S, BetaS, par, m);

decl rgn=RegiaoK(K, cAtual, beta_sqr);
if(rgn[2]==1){ind_flp3[j-1]=1;}


}  // for i

ind_std3=1-ind_flp3;


iind_flp=mtxfind(ind_flp3,1);
n_flp=rows(iind_flp);

iind_std=mtxfind(ind_std3,1);
n_std=rows(iind_std);

n_flp_lf3=ceil((n_flp/(n_flp+n_std))*NLk[2]);
n_std_lf3=NLk[2]-n_flp_lf3;

temp=(ranshuffle(NLk[2],range(0,NLk[2]-1)'))';
ind_flp_lf3=zeros(NLk[2],1);
if(n_flp_lf3>0){ind_flp_lf3[temp[0:(n_flp_lf3-1)]]=1;}

index_flp_lf=vecindex(ind_flp_lf3);
index_std_lf=vecindex(1-ind_flp_lf3);

ind_lfs_flp=reshape(index_flp_lf,n_flp,1);
ind_lfs_std=reshape(index_std_lf,n_std,1);
temp_flp=fmod(n_flp,n_flp_lf3);
temp_std=fmod(n_std,n_std_lf3);

ind_lfs3=zeros(NSk[2],1);
			 
if(n_flp>0)
{
ind_lfs_flp=ranshuffle(n_flp,ind_lfs_flp)';
ind_lfs3=asng_vl_mt(ind_lfs3, ind_lfs_flp, iind_flp);
}

if(n_std>0)
{
ind_lfs_std=ranshuffle(n_std,ind_lfs_std)';
ind_lfs3=asng_vl_mt(ind_lfs3, ind_lfs_std, iind_std);
}


}

//////////////////////////////////////////////
//////////////////////////////////////////////
//////////////////////////////////////////////


[S1, S2, Ykk, muk, D, lv]=Rk(K, cAtual, S, BetaS, M_Y[][2], lambdaAtual);

IntIF[i-1]=(MuS/nnn)*sumc(S1);
mmuk[i-1][]=muk;


if(i>burn)
{
sD+=D; 
sYk+=Yk;	
smuk+=muk;		  
sS1+=S1;
sS2+=S2;
slambda+=lambdaAtual;
}  

/////////// printings ///////////

if(fmod(i,burn1)==0)
{
println("i=",i);
println(Ib/(i));
println(Il1/(i));
println(Il2/(i));
println(Il3/(i));


println(quantiler(In/(i),range(0,1,0.2)));
println(lambdaAtual|muk,cAtual);


println((ind_stop1_temp/(i))~(ind_stop2_temp/(i))~(ind_stop3_temp/(i)));


}


if(fmod(i,10)==0){println("Time up to iteration ",i," = ",timespan(time));}

if(fmod(i,1)==0){println(i);}


if(i>0 && fmod(i,1000)==0)
{


decl mmm = (i-burn);
decl S1_A = sS1./mmm;

savemat("S1_A.mat",S1_A,1) ;
savemat("lambda.mat",lambda[:i][],1) ;
savemat("BetaS0.mat",BetaS);

}


}	  // end MCMC iterations



decl mmm,S1_A,S2_a, aa, S2_A,bb,b1,b2,b3,b4,b5,media_lambda,D_hat,Yk_hat,muk_hat;

mmm = NumSim - burn;
media_lambda = slambda./mmm;
D_hat = sD./mmm;	  
Yk_hat= sYk./mmm;
muk_hat=  smuk./mmm;


S1_A = sS1./mmm;

S2_a = sS2./mmm; aa=maxr(S2_a);bb=zeros(nL*nU,1);

for(j=0;j<=(nL*nU-1);++j){bb[j][]=vecrindex(S2_a[j][] .== aa[j]);} 

b1=vecrindex(bb .== 0);
b2=vecrindex(bb .== 1);
if(K>2){b3=vecrindex(bb .== 2);}
if(K>3){b4=vecrindex(bb .== 3);}
if(K>4){b5=vecrindex(bb .== 4);}

S2_A = zeros(nL*nU,1);
S2_A[b1] = media_lambda[0];
S2_A[b2] = media_lambda[1];
if(K>2){S2_A[b3] = media_lambda[2];}
if(K>3){S2_A[b4] = media_lambda[3];}
if(K>4){S2_A[b5] = media_lambda[4];}



decl time1 = timespan(time);

println("Total time = ",time1);

println("Lambda time = ",time2);
println("GP + N time = ",time3);

savemat("TotalTime.mat",time1,1) ;
savemat("Times.mat",time2|time3,1) ;


savemat("lambda.mat",lambda,1) ;			// MCMC chain for lambda

savemat("S1_A.mat",S1_A,1) ;				// output to produce IF estimate plot
savemat("S2_A.mat",S2_A,1) ;				// output to produce IF estimate plot
savemat("media_lambda.mat",media_lambda) ;	// posterior mean of lambda

println("Il1=",Il1/NumSim);					    // acceptance rate of lambda
println("Il2=",Il2/NumSim);					    // acceptance rate of lambda
println("Ib=",Ib/NumSim);					// acceptance rate 1 of beta
println("In_avg=",meanr(In/NumSim)~sqrt(varr(In/NumSim)));					// acceptance rate of N

savemat("Il1.mat",Il1/NumSim,1) ;				// acceptance rate of lambda
savemat("Il2.mat",Il2/NumSim,1) ;				// acceptance rate of lambda
savemat("Il3.mat",Il3/NumSim,1) ;				// acceptance rate of lambda
savemat("Ib.mat",Ib/NumSim,1) ;			// acceptance rate of beta
savemat("In.mat",In/NumSim,1) ;			// acceptance rate of N
savemat("loglik.mat",loglik,1) ;			// MCMC chain for log-likelihood

savemat("BetaS0.mat",BetaS);					 // last iteration of Beta at S
savemat("M_Y.mat",M_Y);						 // last iteration of Beta at Y
savemat("M_C.mat",M_C);						 // last iteration of N and Beta at N

savemat("cont_data.mat",cont_w1~cont_w2~cont_phi1~cont_phi2,1);

} 
