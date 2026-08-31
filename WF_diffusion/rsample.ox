
rsample(const r, const c, const m)
 {
  decl u,ind;

  u = ranu(r,c);

  ind=floor(u*m);
  
  return ind;
 }
