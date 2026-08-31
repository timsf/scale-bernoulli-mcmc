
mtxfind(const A,const x)
 {
  decl r,c,MM,ind,rr,iind,cc;

  r=rows(A);
  c=columns(A);

  MM=vecr(A);

  ind=vecindex(MM,x);

  rr=floor(ind/c);
  cc=fmod(ind,c);

  iind=rr~cc;

  
  return iind;
 }
