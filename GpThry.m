/////////////////////////////////////////////////
// All of our code was run on a machine with the
// following specifications.
// CPU: 2.9 GHz 6-Core Intel Core i9
// Memory: 32 GB 2400 MHz DDR4
// OS: macOS Montery Version 12.5.1
// Magma Version: V2.27-3
/////////////////////////////////////////////////

/////////////////////////////////////////////////
// Given i = 1 or 2, positive integers d_1, d_2
// and a subgroup D of GL(d_1,R) x GL(d_2,R), for
// some ring R, output the natural projection map 
// D -> GL(d_i,R) onto the i-th component.
//////////////////// EXAMPLE ////////////////////
// > D := DirectProduct(GL(2,2),GL(2,2));
// > Pr(D,2,2,1);
// Homomorphism of MatrixGroup(4, GF(2)) of order 2^2 * 3^2 into GL(2, GF(2)) induced by
//    [1 1 0 0]
//    [0 1 0 0]
//    [0 0 1 0]
//    [0 0 0 1] |--> [1 1]
//    [0 1]
//    [0 1 0 0]
//    [1 0 0 0]
//    [0 0 1 0]
//    [0 0 0 1] |--> [0 1]
//    [1 0]
//    [1 0 0 0]
//    [0 1 0 0]
//    [0 0 1 1]
//    [0 0 0 1] |--> [1 0]
//    [0 1]
//    [1 0 0 0]
//    [0 1 0 0]
//    [0 0 0 1]
//    [0 0 1 0] |--> [1 0]
//    [0 1]
/////////////////////////////////////////////////
Pr := function(D,d1,d2,i)
  if i eq 1 then
    G := GL(d1,BaseRing(D));
    return hom<D->G | [<d,Submatrix(d,1,1,d1,d1)> : d in Generators(D)]>;
  end if;
  if i eq 2 then
    G := GL(d2,BaseRing(D));
    return hom<D->G | [<d,Submatrix(d,d1+1,d1+1,d2,d2)> : d in Generators(D)]>;
  end if;
end function;

/////////////////////////////////////////////////
// Given subgroups G1 and G2 of GSp(2g_1,R) and
// GSp(2g_2, R) (respectively), for some positive 
// integers g_1, g_2 and some ring R, output a
// group that is conjugate to G1 x_mult G2 inside
// G1 x G2, realized as a subgroup of 
// GL(2g_1+2g_2, R), together with the natural
// projection maps to G1 and G2. This function 
// throws an error if it is unable to uniquely
// determine G1 x_mult G2 (up to conjugacy). While 
// this error does not occur in any of the cases 
// where we call the function, a more general 
// version should someday be written.
//////////////////// EXAMPLE ////////////////////
// > G1 := CSp(4,2);
// > G2 := GL(2,2);
// > Delta, pr1, pr2 := ConstructDelta(G1,G2);
// > Delta;
// MatrixGroup(6, GF(2)) of order 2^5 * 3^3 * 5
// Generators:
//    [0 0 1 1 0 0]
//    [1 1 1 0 0 0]
//    [1 1 0 0 0 0]
//    [1 0 0 0 0 0]
//    [0 0 0 0 1 0]
//    [0 0 0 0 0 1]
//
//    [1 0 0 0 0 0]
//    [1 1 0 0 0 0]
//    [0 0 1 0 0 0]
//    [0 0 1 1 0 0]
//    [0 0 0 0 1 0]
//    [0 0 0 0 0 1]
//
//    [1 0 0 0 0 0]
//    [0 1 0 0 0 0]
//    [0 0 1 0 0 0]
//    [0 0 0 1 0 0]
//    [0 0 0 0 1 1]
//    [0 0 0 0 0 1]
//
//    [1 0 0 0 0 0]
//    [0 1 0 0 0 0]
//    [0 0 1 0 0 0]
//    [0 0 0 1 0 0]
//    [0 0 0 0 0 1]
//    [0 0 0 0 1 1]
/////////////////////////////////////////////////
ConstructDelta := function(G1,G2)
  R := BaseRing(G1);
  NumUnits := #UnitGroup(R);
  GG := DirectProduct(G1,G2);
  pr1 := Pr(GG,Dimension(G1),Dimension(G2),1);
  pr2 := Pr(GG,Dimension(G1),Dimension(G2),2);
  BigSubs := LowIndexSubgroups(GG,NumUnits);
  GoodSubs := [H : H in BigSubs | Index(GG,H) eq NumUnits and pr1(H) eq G1 and pr2(H) eq G2];
  assert #GoodSubs eq 1;
  return GoodSubs[1], pr1, pr2;
end function;

/////////////////////////////////////////////////
// Given a matrix h in GL(d,ell), for some
// positive integer d and some prime ell, output 
// dim_1(h), i.e., the dimension of the 
// 1-eigenspace of h.
//////////////////// EXAMPLE ////////////////////
// > h := Matrix(GF(3),[[1,1],[0,1]]);
// > Dim1(h);
// 1
// > h := Matrix(GF(5),[[1,0,0,0],[0,1,2,0],[0,0,1,0],[0,0,0,1]]);
// > Dim1(h);
// 3
/////////////////////////////////////////////////
Dim1 := function(h)
  G := Parent(h);
  ell := #BaseRing(G);
  return Integers() ! Log(ell,#Kernel(Matrix(h)-Matrix(One(G))));
end function;

/////////////////////////////////////////////////
// Given a subgroup D of GL(d,ell) x GL(d,ell),
// for some positive integer d and some prime ell,
// output the set of pairs <dim_1(h_1), dim_1(h_2)> 
// as <h_1,h_2> ranges over all elements of D.
//////////////////// EXAMPLE ////////////////////
// > G := GL(2,3);
// > Delta, pr1, pr2 := ConstructDelta(G,G);
// > Dim1Sub(Delta);
// { <2, 2>, <1, 2>, <0, 1>, <0, 2>, <0, 0>, <1, 1>, <2, 0>, <2, 1>, <1, 0> }
/////////////////////////////////////////////////
Dim1Sub := function(D)
  d := Integers() ! (Dimension(D)/2);
  return {<Dim1(Submatrix(h,1,1,d,d)),Dim1(Submatrix(h,d+1,d+1,d,d))> : h in D};
end function;

/////////////////////////////////////////////////
// Given subgroups G1 and G2 of GSp(2g_1,R) and
// GSp(2g_2, R) (respectively), for some positive 
// integers g_1, g_2 and some ring R, output the set
// of ordered pairs (d_1,d_2) such that there 
// exists (h_1,h_2) in G1 x_mult G_2 for which
// (dim_1(h_1), dim_1(h_2)) = (d_1, d_2) yet for
// all proper subgroups K of G1 x_mult G_2
// that surject onto G1 and G2 via the respective
// projection maps, there does not exist a pair
// (k_1,k_2) in K such that 
// (dim_1(k_1), dim_1(k_2)) = (d_1, d_2).
////////////// Used in Lemma 6.2 ////////////////
// > DistinguishedDims(GL(2,2),GL(2,2));
// { <0, 1>, <1, 2>, <1, 0>, <2, 1> }
// > DistinguishedDims(GL(2,3),GL(2,3));
// { <1, 2>, <2, 1> }
/////////////////////////////////////////////////
DistinguishedDims := function(G1,G2)
  Delta, pr1, pr2 := ConstructDelta(G1,G2);
  MaxSubs := [H`subgroup : H in MaximalSubgroups(Delta)];
  SubDims := {};
  for H in MaxSubs do
    if H ne Delta and pr1(H) eq G1 and pr2(H) eq G2 then
      SubDims := SubDims join Dim1Sub(H);
    end if;
  end for;
  return Dim1Sub(Delta) diff SubDims;
end function;

/////////////////////////////////////////////////
// Give a positive integer g and a prime ell,
// output `true' if Corollary 3.8 holds for the 
// given pair (g,ell) and `false' otherwise. 
////////////// Used in Example 6.4 //////////////
// > time Cor38Holds(2,3);
// true
// Time: 16.980
/////////////////////////////////////////////////
Cor38Holds := function(g,ell)
  G := CSp(2*g,ell);
  Delta, pr1, pr2 := ConstructDelta(G,G);
  MaxSubs := [H`subgroup : H in MaximalSubgroups(Delta)];
  for H in MaxSubs do
    if pr1(H) eq G and pr2(H) eq G then
      for h in H do
        if not Trace(pr1(h)) in {Trace(pr2(h)), -Trace(pr2(h))} then
          return false;
        end if;
      end for;
    end if;
  end for;
  return true;
end function;

/////////////////////////////////////////////////
// Given distinct positive integers g_1, g_2 and
// a prime ell, output `true' if Lemma 3.9 also 
// holds for (g_1,g_2,ell) and `false' otherwise. 
////////////// Used in Example 6.5 //////////////
// > Lemma39Holds(2,1,3);
// true
/////////////////////////////////////////////////
Lemma39Holds := function(g1,g2,ell)
  G1 := CSp(2*g1,ell);
  G2 := CSp(2*g2,ell);
  Delta, pr1, pr2 := ConstructDelta(G1,G2);
  MaxSubs := [H`subgroup : H in MaximalSubgroups(Delta)];
  for H in MaxSubs do
    if pr1(H) eq G1 and pr2(H) eq G2 then
      return false;
    end if;
  end for;
  return true;
end function;

/////////////////////////////////////////////////
// The following code is used in Example 6.5 to
// verify the surjectivity of the mod 2 Galois
// representation of A x E (as in the example). 
// It does so by computing the order of the 
// Galois group of the 2-division field of A x E.
/////////////////////////////////////////////////
// > R<x> := PolynomialRing(Rationals());
// > f1 := x^6-2*x^4+2*x^3+5*x^2+2*x+1;
// > f2 := x^3-17977*x-927735;
// > #GaloisGroup(f1);
// 720
// > #GaloisGroup(f2);
// 6
// > #GaloisGroup(f1*f2);
// 2160
/////////////////////////////////////////////////

/////////////////////////////////////////////////
// The following code uses Andrew V. Sutherland's
// 'galrep' (https://math.mit.edu/~drew/galrep/).
// The code is used in Example 6.2 to verify the
// surjectivity of the mod 3 Galois 
// representation of E_1 x E_2, where E_1 and E_2
// are as in the example. To do so, it calls
// Sutherland's EpSig function, and looks at the
// third entry, which is the dimension of the 
// 1-eigenspace of the image of Frob_p. We use
// p = 73; surjectivity follows from Lemma 6.2.
/////////////////////////////////////////////////
// > ChangeDirectory("/Users/.../galrep"); // User should modify this path to the path of 'galrep' on their machine.
// > load "algorithms.m";
// Loading "algorithms.m"
// Loading "subgroups.m"
// Loading "distinguish.m"
// > E1 := EllipticCurve([0, 0, 0, 1, 10]);
// > E2 := EllipticCurve([0, 0, 0, -362249, 165197113]);
// > EpSigs(ChangeRing(E1,GF(73)),3)[1][3];
// 1
// > EpSigs(ChangeRing(E2,GF(73)),3)[1][3];
// 2
/////////////////////////////////////////////////
