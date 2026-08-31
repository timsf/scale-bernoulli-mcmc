
dmatrix(const s)

{
  decl n,tt1,dt1,tt2,dt2,WW;

 n=rows(s);

 tt1=s[][0]*ones(1,n);
 tt2=s[][1]*ones(1,n);

 WW=(((tt1-tt1').^2)+((tt2-tt2').^2)).^(1/2);

 return WW;
}