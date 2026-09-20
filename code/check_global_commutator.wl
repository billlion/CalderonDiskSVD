(* Direct symbolic check of the global kernel identity.
   Run in Mathematica and expect the result 0. *)
ClearAll[r, rho, psi, den, ker, Pz, Pw];
den = 1 - 2 r rho Cos[psi] + r^2 rho^2;
ker = 1/(Pi den);
Pz[u_] := 1/(4 r) D[r (1 - r^2)^2 D[u, r], r]
  + 1/4 (r^2 + 1/r^2) D[u, {psi, 2}] + r^2 u;
Pw[u_] := 1/(4 rho) D[rho (1 - rho^2)^2 D[u, rho], rho]
  + 1/4 (rho^2 + 1/rho^2) D[u, {psi, 2}] + rho^2 u;
FullSimplify[Pz[ker] - Pw[ker],
 Assumptions -> 0 < r < 1 && 0 < rho < 1 && Element[psi, Reals]]
