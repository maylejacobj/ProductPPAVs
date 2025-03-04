# ProductPPAVs

This repository is associated with the article *An effective open image theorem for products of principally polarized abelian varieties* [(preprint)](https://arxiv.org/abs/2212.11472) [(journal)](https://www.sciencedirect.com/science/article/abs/pii/S0022314X25000472) by Jacob Mayle and Tian Wang.

The function `FindLambda` in `FindLambda.sage` is a Sage implementation of Algorithm 6.1. It is used for bounding the largest nonsurjective prime associated with the product of two Jacobians of hyperelliptic curves, as described in the statement of Algorithm 6.1. The file `GpThry.m` contains Magma code that is refered to in Section 6 for checking various claims. 

Installation instructions for `FindLambda.sage`:
1. Install the latest version of SageMath from [https://doc.sagemath.org/html/en/installation/index.html](https://doc.sagemath.org/html/en/installation/index.html).
2. Download and run `FindLambda.sage` from this repository.


Installation instructions for `GpThry.m`:
1. Install the latest version of Magma from [http://magma.maths.usyd.edu.au/magma/](http://magma.maths.usyd.edu.au/magma/).
2. Download and run `GpThry.m` from this repository.
3. For the last code snippet only, Andrew V. Sutherland's `galrep` is used. For this, download the files from [https://math.mit.edu/~drew/galrep/](https://math.mit.edu/~drew/galrep/) and move them into a folder called `galrep`. Then modify the path appearing in Line 248 of `GpThry.m` to match the path of the `galrep` folder on your machine.

We welcome any questions, comments, or suggestions.
