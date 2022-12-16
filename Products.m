// Given a matrix g, return the dimension of the 1-eigenspace of g
Dim1 := function(g)
  G := Parent(g);
  ell := #BaseRing(G);
  return Integers() ! Log(ell,#Kernel(Matrix(g)-Matrix(One(G))));
end function;

// Given a subgroup D <= GL(2,ell) x GL(2,ell), return the set of pairs <dim_1 h_1, dim_1 h_2> as <h_1,h_2>
// varies over all elements of D.
Dim1Sub := function(D)
  return {<Dim1(Submatrix(h,1,1,2,2)),Dim1(Submatrix(h,3,3,2,2))> : h in D};
end function;

// Given a subgroup D <= GL(2,ell) x GL(2,ell) and i = 1 or 2, return the i-th projection map D -> GL(2,ell)
Pr := function(D,i)
  G := GL(2,BaseRing(D));
  return hom<D->G | [<d,Submatrix(d,2*i-1,2*i-1,2,2)> : d in Generators(D)]>;
end function;

// Given a subgroup G of GL(2,ell), return the group G x_det G, realized in GL(4,ell)
ConstructDelta := function(G)
  ell := #BaseRing(G);
  GG := DirectProduct(G,G);
  pr1 := Pr(GG,1);
  pr2 := Pr(GG,2);
  L := LowIndexSubgroups(GG,ell-1);
  for H in L do
    if #H * (ell - 1) eq #GG and #pr1(H) eq #G and #pr2(H) eq #G then
      return H;
    end if;
  end for;
end function;

// This function verifies Lemma 6.1.
// > DistinguishedDims(2);
// { <0, 1>, <1, 2>, <1, 0>, <2, 1> }
// > DistinguishedDims(3);
// { <1, 2>, <2, 1> }
DistinguishedDims := function(ell)
  G := GL(2,ell);
  Delta := ConstructDelta(G);
  pr1 := Pr(Delta,1);
  pr2 := Pr(Delta,2);
  L := Subgroups(Delta);
  S := {};
  for H in L do
    if #H`subgroup ne #Delta and #pr1(H`subgroup) eq #G and #pr2(H`subgroup) eq #G then
      S := S join Dim1Sub(H`subgroup);
    end if;
  end for;
  return Dim1Sub(Delta) diff S;
end function;

// The code below is a companion to Example 6.2. 
E1 := EllipticCurve([0, 0, 0, 1, 10]);
E2 := EllipticCurve([0, 0, 0, -362249, 165197113]);

// Verify that 17 is not a bad prime for E1 or E2.
BadPrimes(E1); 
// [ 2, 13 ]
BadPrimes(E2);
// [ 2, 13, 19 ]

// Compute the Frobenius trace for E1 and E2 at 17
TraceOfFrobenius(E1,17);
// 6
TraceOfFrobenius(E2,17);
// -7

// Verify that a_p(E_1) = LegendreSymbol(-1,p) * a_p(E_2) (mod 13) for all 
// primes p of good reduction for E_1 and E_2 up to p <= 1000.
for p in PrimesUpTo(1000) do
  if not p in [2,13,19] then
    p, TraceOfFrobenius(E1,p) mod 13 eq (LegendreSymbol(-1,p) * TraceOfFrobenius(E2,p)) mod 13;
  end if;
end for;

// Verify surjectivity using Lemma 6.1 by calling Sutherland's EpSig function, and
// looking at the third entry, which is the dimension of the 1-eigenspace of the image
// of Frob_p; we consider p = 73 below.
EpSigs(ChangeRing(E1,GF(73)),3)[1][3];
// 1
EpSigs(ChangeRing(E2,GF(73)),3)[1][3];
// 2