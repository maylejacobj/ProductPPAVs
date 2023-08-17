#################################################
## All of our code was run on a machine with the
## following specifications.
## CPU: 2.9 GHz 6-Core Intel Core i9
## Memory: 32 GB 2400 MHz DDR4
## OS: macOS Montery Version 12.5.1
## SageMath Version: 10.0
#################################################


def FindLambda(C1,C2,B=1000,Verbose=False):
    r"""
    Return the nonnegative integer Lambda (as in Algorithm 6.1).

    INPUT:

    - ``C1`` -- hyperelliptic curve

    - ``C2`` -- hyperelliptic curve

    - ``bound`` -- integer (default: `1000`)

    - ``verbose`` -- boolean (default: `False`)

    OUTPUT: the nonnegative integer Lambda (as in Algorithm 6.1)

    EXAMPLES:

    Example 6.3::

        sage: E1 = EllipticCurve([0, 0, 0, 1, 10])
        sage: E2 = EllipticCurve([0, 0, 0, -362249, 165197113])
        sage: FindLambda(E1,E2,Verbose=True)
        (13, [[17, 6, -7]])

    Example 6.4::

        sage: R.<x> = PolynomialRing(QQ)
        sage: C1 = HyperellipticCurve(R([-1, 1, -3, 2, -3, 1, -1]), R([1, 1, 0, 1]))
        sage: C2 = HyperellipticCurve(R([0, 0, -1, 2, 0, -2, 1]), R([1, 1]))
        sage: FindLambda(C1,C2,Verbose=True)
        (1, [[5, 2, -1], [17, 3, 1]])

    Example 6.5::

        sage: R.<x> = PolynomialRing(QQ)
        sage: C = HyperellipticCurve(R([0, 0, 1, 0, -1]), R([1, 1, 0, 1]))
        sage: E = EllipticCurve([0, 0, 0, -17977, -927735])
        sage: FindLambda(C,E,Verbose=True)
        (1, 'The genera of the curves are different.')

    """
    if C1.genus() != C2.genus(): 
        if Verbose: return 1, "The genera of the curves are different."
        return 1
    Lambda = 0
    VerboseOut = []
    for p in primes(B):
        try:
            C1p = C1.change_ring(GF(p))
            C2p = C2.change_ring(GF(p))
            frob1 = C1p.frobenius_polynomial()
            frob2 = C2p.frobenius_polynomial()
            a1p = -frob1.coefficients(sparse=False)[frob1.degree()-1]
            a2p = -frob2.coefficients(sparse=False)[frob2.degree()-1]
            M = (a1p-a2p)*(a1p+a2p)
            G = gcd(Lambda,M)
            if (Lambda == 0 and G != 0) or G < Lambda:
                VerboseOut.append([p,a1p,a2p])
                Lambda = G
                if Lambda == 1: 
                    if Verbose: return Lambda,VerboseOut
                    return Lambda
        except:
            continue
    if Verbose: return Lambda,VerboseOut
    return Lambda
