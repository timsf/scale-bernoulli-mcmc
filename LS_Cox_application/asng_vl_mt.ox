
asng_vl_mt(A, const x, const pos)
 {

 decl i;

for (i = 0; i < rows(x); ++i)
{
 A[pos[i][0]][pos[i][1]] = x[i];
}
  
  return A;
 }
