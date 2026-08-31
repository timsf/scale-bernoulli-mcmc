
locc(const s, const nL, const nU, const SiL, const SiU, const S, const L, const U, const mm)

{
  decl m0,m1,ind0,ind1,ind,dist;
  
  ind0=(s[0]|SiL)~(ones(1,1)|zeros(nL+2,1));
  ind1=(s[1]|SiU)~(ones(1,1)|zeros(nU+2,1));

  ind0=sortbyc(ind0,0);
  ind1=sortbyc(ind1,0);

  ind0=vecindex(ind0[][1]);
  ind1=vecindex(ind1[][1]);

  ind0=range(max(0,ind0-5),min(nL-1,ind0+2))';
  ind1=range(max(0,ind1-5),min(nU-1,ind1+2))';


  m0=rows(ind0);
  m1=rows(ind1);
  
  ind0=sortc(reshape(nU*ind0,m0*m1,1));
  ind1=reshape(ind1,m0*m1,1);

  ind=ind0+ind1;

  dist=dmatrix2(s|S[ind][]);

  ind=ind[(sortcindex(dist))[:(mm-1)]];

  
 return ind;
}