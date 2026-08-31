

DCBF(const out_lfs)
{


decl n=rows(out_lfs);

decl LL=round(log(n)/log(2));

decl ind_lfs=zeros(n,1);

decl temp,ind0_temp,i;

decl tree=out_lfs;

decl ind_dif=fabs(tree[range(0,(n-2),2)']-tree[range(1,(n-1),2)']);

ind_lfs[range(0,(n-2),2)']=ind_dif;
ind_lfs[range(1,(n-1),2)']=ind_dif;

decl ind2=vecindex(ind_dif,0);

decl tree2=zeros(n/2,1)+2;
tree2[ind2]=tree[2*ind2];

tree=tree2;

decl l=1;

decl cc=0;

decl ind0a,ind0b,ind0,m,ind,ind1;

while(cc==0 && l<=(LL-1))
{
n=rows(tree);

ind0a=vecindex(tree[range(0,(n-2),2)'],<0;1>);
ind0b=vecindex(tree[range(1,(n-1),2)'],<0;1>);
ind0=intersection(ind0a,ind0b)';

m=rows(ind0);
if(m==0){cc=1;}

else
{

ind_dif=fabs(tree[2*ind0]-tree[2*ind0+1]);

ind1=vecindex(ind_dif);
ind2=vecindex(ind_dif,0);

tree2=zeros(n/2,1)+2;

tree2[ind0[ind2]]=tree[2*ind0[ind2]];

tree=tree2;

temp=rows(ind1);
if(temp>0)
{
for(i=1;i<=temp;++i)
{
ind0_temp=ind0[ind1[i-1]];
ind=((2^(l+1))*(ind0_temp+1))-1;
ind_lfs[(ind-(2^(l+1)-1)):ind]=1;
}
}

l+=1;


} //else rows(ind0)=0


}  // while cc=0


return(ind_lfs);

}




