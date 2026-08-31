
dmatrix2(const s)

{
  decl n,tt1,dt1,tt2,dt2,WW;

 n=rows(s);

 WW=(((s[0][0]-s[1:][0]).^2)+((s[0][1]-s[1:][1]).^2)).^(1/2);

 return WW;
}