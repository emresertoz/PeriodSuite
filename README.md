# PeriodSuite

Software to compute periods of hypersurfaces, with additional functions to make use of computed periods (e.g. for Picard/Hodge group computation).

This is an implementation of the algorithm introduced in [Computing periods of hypersurfaces](https://arxiv.org/abs/1803.08068).

Picard rank, or more generally Hodge rank, computations are now available with the caveat that they lack certification. A detailed analysis of the reliability of the Hodge rank computations appears in [A numerical transcendental method in algebraic geometry](https://arxiv.org/abs/1811.10634).

## Requirements

You need a Unix based system, [SageMath](http://www.sagemath.org/) version 9+, and [Magma](https://magma.maths.usyd.edu.au/). If you must use SageMath 8+ then checkout and use the branch "2019" instead, however I will not keep that branch up to date. 

## Installation

1) Clone this repository to a directory. Replace `/path/to/PeriodSuite` below with the path you want to install this package.

    ```
    git clone https://github.com/emresertoz/PeriodSuite.git /path/to/PeriodSuite
    ```

2) Change the current directory to the cloned repository and type `./make`.

    ```
    cd /path/to/PeriodSuite && ./make
    ```

The installation is now complete. To use PeriodSuite, open Magma and attach the file `suite.mag` before each session, see below. You may want to add this to your `startup` file. To check if everything is working type `X:=test();`.

    Magma> Attach("suite.mag");
    Magma> X:=test();

#### Remarks:

- The newest, most stable versions of all software used are recommended and often mandatory: SageMath \(v9 or higher\) and Magma \(v2.23.8 or higher\). 

- If anything goes wrong please contact me at `emresertoz` at `gmail.com` and I will try to fix the issue.

- If you must use SageMath v8 then type `git checkout 2019` after the first step (but after you `cd` into the PeriodSuite folder). Then type `./make`. This is a commit frozen in time, I do not plan on keeping it up to date.

******

# Usage

Open a Magma session and attach `suite.mag`:

    Magma> Attach("/path/to/suite/suite.mag");

Write the polynomial of the hypersurface of interest:

    Magma> P<x,y,z,w>:=PolynomialRing(Rationals(),4);
    Magma> f:=-4*x^4 + 7*x*y^3 + 2*z^3*w + 6*z^2*w^2 + 6*w^4;

Now you may call the function `PeriodHomotopy` on this polynomial as follows:

    Magma> X:=PeriodHomotopy(f);

To increase precision, set the `precision` parameter to the desired number of digits. The default value is 100.

    Magma> X:=PeriodHomotopy(f : precision:=500);

To access the periods of the hypersurface and the number of digits to which the periods are correct, type:

    Magma> X`periods;
    Magma> X`precision;

 The integral basis for homology is only implicit. However, the basis for cohomology used for the period matrix is explicit and can be accessed by typing:

    Magma> X`cohomBasis;

### Hodge rank

The output `X` of `PeriodHomotopy` can be fed in to `HodgeLattice`, which will return the record `X` after computing the (virtual) Hodge lattice and its polarization. We recommend the following usage:

    Magma> X:=HodgeLattice(X);

Here *virtual* refers to the fact that the computed lattice may in principle differ from the correct Hodge lattice. This will happen if the precision of the period matrix is too low. 

The main attributes of the (virtual) Hodge lattice of `X` can be accessed as follows:

    Magma> Lambda:=X`hodgeLattice;
    Magma> Lambda`lattice;
    Magma> Lambda`polarization;

We offer the following shortcut:

    Magma> X`lat;
    Magma> X`pol;

## Important tip

By default, the deformation path to the target hypersurface is created randomly according to simple heuristics. If a computation takes more than a few seconds, kill the process and run `PeriodHomotopy` again. Repeating this many times increases your chances of discovering a deformation path for which the function terminates quickly.

This is very important! The result of a computation may take hours or a fraction of a second, depending on how good the deformation path is. For harder computations, it is worth spending time on the discovery of a good path.

******

# Advanced usage

There are many options and additional functions. We list the most important ones here. Recall the following definition.

A polynomial of the form `c0*x0^d + ... + cn*xn^d` where `c0,..,cn` are non-zero, is called of *Fermat type*.

### Output

The output `X` of `PeriodHomotopy` is a record of type `PeriodBundle`. The following attributes of `X` should be of interest:

```
X`deformation : contains the deformation path used in period homotopy. 
X`cohomBasis : contains the polynomial part of the forms whose periods are to be computed, one list per family.
X`primitivePeriods : the primitive periods of X.
X`periods : the periods of X. When X is odd dimensional this coincides with the primitive periods of X.
```

Auxiliary data are written in to the directory `ode_storage/incinerator`, where it is read by `integrator.sage`. As the name suggests, files in this folder are not safe as the folder will be cleared of content at the beginning of each execution of `PeriodHomotopy`. 

### Path finder

The integrator needs a path in the complex plane for every one of the systems it integrates. By default it will try the straight path 0 to 1, it works most of the time. However, if there are singularities along the way the computations might not terminate. In this case we recommend setting the option `pathfinder:=1`:

    Magma> PeriodHomotopy(f: pathfinder:=1);

### Straight deformation

The following will turn off finding a randomized, monomial deformation path to target hypersurface and will go straight to `f` from the best possible Fermat type hypersurface.

    Magma> PeriodHomotopy(f:straight:=true);

If you want to go straight from a hand picked Fermat type hypersurface `g` to `f` then type: 

    Magma> X := PeriodHomotopy([g,f]);

### Custom deformation path

If you have a sequence of smooth, degree `d` polynomials `f0,f1,..,fr` where `f0` is of Fermat type, then you can type in this sequence into `PeriodHomotopy`. Here is an example with a sequence of 5 elements:

    Magma> PeriodHomotopy([f0,f1,f2,f3,f4,f5]);

The period homotopy will be performed on this sequence with `f5` being the target hypersurface and `fi`'s the only intermediate hypersurfaces.

### Monodromy / Picard-Lefschetz

If you would like to compute the monodromy operator of going counterclockwise around a singular hypersurface, then you can do so using PeriodHomotopy which uses the periods of nearby fibers. If you give a singular target hypersurface, this intention is understood and the necessary computations are performed. For example:

    Magma> P3<x,y,z,w>:=PolynomialRing(Rationals(),4);
    Magma> def:=[ x^4 - y^4 - z^4 - w^4, x^4 - y^4 - y^2*z^2 - z^4 - w^4, x^4 - y^2*z^2 ];
    Magma> X:=PeriodHomotopy(def : bound_pole_order:=false, precision:=100);
    Magma> X:=PicardLefschetz(X);
    Magma> X`monodromy;
    [ 0  0  0  0  0  0 -1  0 -1  0  0  1  0  1  1  1  1  0  0 -1  0]
    [ 0  0  0 -1  0  0 -1  0 -1  0  0  1  0  1  1  2  0  0  1  0  1]
    [ 0  1  1  0 -1 -1 -2 -1 -2 -1  0  2  0  3  0  1  1  1  0  0  0]
    [ 0  0  0  1  0  0 -1  0 -1 -1 -1  1 -1  0  1  1  1  1  0 -1  0]
    [ 0  0  0  0  0  0 -1  0 -1 -1  0  1  0  1  1  1  1  0  0  0  1]
    [ 0  0  0 -1  0  0 -1  0 -1  0  0  1  0  2  0  1  0  0  1  0  1]
    [ 1  0  1  0 -1 -1  0 -1  0 -1 -1  0 -1 -1  0  1 -1  1  0  1  1]
    [ 0  0  0  0  0  0 -1  0 -1 -1  0  1  0  1  1  1  1  0  1 -1  0]
    [ 1  2  1  1 -2 -2 -1 -2 -1 -2 -1  1 -1  1 -1  0  0  2  0  1 -1]
    [ 0  1  1  0 -1 -1 -1 -1 -1  0  0  1  0  1  0  1  0  1 -1  0 -1]
    [ 0  0  0  1  0  0 -1  0 -1 -1  0  1  0  1  1  0  1  0  0 -1  0]
    [ 1  1  1  1 -2 -1 -1 -1 -1 -2 -1  1 -1  1  0  0  0  2 -1  1  0]
    [ 0 -1 -1  0  1  1  0  1  0  0  0  0  0 -1  1  1  1 -1  1 -1  1]
    [ 1  1  1  0 -1 -2  0 -1  0 -1 -1  0 -1  0 -1  0 -1  1  0  1  0]
    [ 0  1  1  0 -1 -1 -1 -1  0  0  0  0  0  1  0  0 -1  1 -1  1  0]
    [ 1  0  0  0  0  0  0 -1  0 -1 -1  0 -1 -1  0  1  0  1  1  0  0]
    [ 0  1  1  0 -1 -1  0 -1 -1  0  0  0  0  1 -1  0  0  1  0  1 -1]
    [ 0  0  0  1  0  0  0  0  0 -1  0  0  0 -1  1  0  1  0  0 -1  0]
    [ 0  1  1  1 -1 -1 -1 -1 -1 -1 -1  1  0  1  0  0  0  1  0  0 -1]
    [ 0  1  1  1 -1 -1 -1 -1 -1 -1  0  0  0  2  0  0  0  1 -1  1 -1]
    [ 0  0  0  0  0  0  0  0  0  0  0  0 -1 -1  0  1  0  0  0  0  1]

If you are reading this, then you are on an experimental branch. Please let me know of any bugs you encounter.
