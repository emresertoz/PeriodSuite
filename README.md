# period-suite

Software to compute periods of hypersurfaces, with additional functions to make use of computed periods.

This is an implementation of the algorithm introduced in [Computing periods of hypersurfaces](https://arxiv.org/abs/1803.08068).

Picard rank, or more generally Hodge rank, computations are now available with the caveat that they lack certification. A detailed analysis of the reliability of the Hodge rank computations will appear in an upcoming paper.

## Installation

Two parts are needed to compute periods. One part, provided here, has everything except for the numerical integrator. The second part, the numerical integrator, requires the installation of the developmental branch of the `ore_algebra` package, which we explain below.

#### Remarks:

- This software will only work on Unix based systems, e.g., Linux and Mac. 
- If you do not have the text processor `awk` already, please install it first.

### Part 1 - Picard-Fuchs equations

1) Clone this repository to a directory.

    git clone https://github.com/period-suite/period-suite.git /path/to/period-suite

2) From the directory into which you have cloned `period-suite`, open `suite.mag` with a text editor and set `pathToSuite` to the *complete* path of this directory. Make sure that the path name ends with a forward slash `/`. 

3) Open `integrator.sage` with a text editor and set the value of `pathToSuite` as above.

(Optional) Open `integrator.ipynb` and `find_path.ipynb` to set the value of `pathToSuite` as before. For instance:

    sage --notebook=jupyter integrator.ipynb

    sage> pathToSuite="/path/to/period-suite";

### Part 2 - Integration

We will need Mezzarrobba's `analytic` extension of the `ore_algebra` package. You may download and install the current version by typing:

    sage -pip install --user git+https://github.com/mkauers/ore_algebra.git@analytic

In case this doesn't work, check the newest instructions on how to install `ore_algebra--analytic` [here](http://marc.mezzarobba.net/code/ore_algebra-analytic/).

### Recommendations

The newest, most stable versions of all software used are recommended and often mandatory.

- If you have a version of `ore_algebra-analytic` older than March 2017, please renew it for much better performance.
- Make sure you have the latest version of SageMath \(v8.1 or higher\) and Magma \(v2.23.8 or higher\). 
- If you do not use `PathFinder2` then Magma v.2.23.1 or higher will be sufficient.

### Optional evaluation-interpolation scheme

For advanced usage, clone Lairez' evaluation-interpolation scheme for computing Picard-Fuchs equations. The `linhomotopy` branch of his repository is compatible with `period-suite`. We recommend cloning the `linhomotopy` branch directly into a folder within `period-suite`:

    git clone -b linhomotopy https://github.com/lairez/periods.git lairez-periods 

Now edit the variables `lairezAvailable` and `lairezSpec` in `suite.mag` accordingly.

******

# Usage

Navigate to the directory in which you have cloned this repository. Open a Magma session and load the file `suite.mag`.

    Magma> load "suite.mag";

Now you are set and can start using the code. We automatically define a few rings to start using the function but also to manipulate the output of the main function `PeriodHomotopy`. In particular `P` is already defined to be the polynomial ring with four variables `x,y,z,w`. These definitions are made in `suite.mag` --- feel free to remove them if you wish.

    Magma> f:=-4*x^4 + 7*x*y^3 + 2*z^3*w + 6*z^2*w^2 + 6*w^4;
    Magma> ph:=PeriodHomotopy([f]);

Note the extra square brackets around `f`! See advanced usage below on how to give sequences to `PeriodHomotopy`.  

The function `PeriodHomotopy` will produce the necessary ODEs and the initial conditions for the numerical integrator. Then, it will call SageMath and start integration.

Once computations are complete, type the following to read in the results of the integration manually:

    Magma> load "lastPeriods";
    Magma> X`primitivePeriods:=periods;
    Magma> X`precision:=precision;

At this point, one could view `X` as representing the target hypersurface, but with the additional data of its periods as well as the means used to compute the periods.

### Hodge rank

Let us suppose you have completed a period homotopy and read in the results of the integration as indicated above. Then `X` can be fed in to `HodgeLattice` which will return the record `X` after completing the primitive periods and computing the (virtual) Hodge lattice.

    Magma> X:=HodgeLattice(X);

Here *virtual* refers to the fact that the computed lattice may in principle differ from the Hodge rank. This will happen if the precision of the period matrix has been computed to low precision. We will offer an analysis of the reliability of this method in an upcoming paper.

The attributes of the (virtual) Hodge lattice of `X` can be accessed as follows:

    Magma> X`hodgeLattice`lattice;
    Magma> X`hodgeLattice`polarization;

## Examples

See `examples.mag` and `demo.mag` for many examples.

## Important tip

By default, the deformation path to the target hypersurface is created randomly according to simple heuristics. If a computation takes more than a few seconds, kill the process and run `PeriodHomotopy` again. Repeating this many times increases your chances of discovering a deformation path for which the function terminates quickly.

This is very important! The result of a computation may take hours or a fraction of a second, depending on how good the deformation path is. For harder computations, it is worth spending time on the discovery of a good path.

******

# Advanced usage

There are many options and additional functions. We list the most important ones here. Recall the following definition.

A polynomial of the form `c0*x0^d + ... + cn*xn^d` where `c0,..,cn` are non-zero, is called of *Fermat type*.

### Output

Let `X` be the output of `PeriodHomotopy` as above. It is a record of type `PeriodBundle`. Once `PeriodHomotopy` has terminated, the following attributes of `X` are the most relevant.

```
X`deformation : contains the deformation path used in period homotopy. 
X`cohomBasis : contains the polynomial part of the forms whose periods are to be computed, one list per family.
X`odes : contains the ODEs corresponding to the periods of the forms in X`cohomBasis.
```

The main bulk of the output, including the initial conditions and transition matrices between families, is written in the file `current.sage`, where it is read by `integrator.sage`.

### Path finder

The integrator needs a path in the complex plane for every one of the systems it integrates. By default it will try the straight path 0 to 1, it works most of the time. However, if there are singularities along the way the computations might not terminate. In this case we recommend setting one of the following options on `PeriodHomotopy`:

    Magma> PeriodHomotopy([f]: pathfinder1:=true);

or 

    Magma> PeriodHomotopy([f]: pathfinder2:=true);

Both have advantages and disadvantages. But we recommend that you use the second one only when the first one doesn't terminate.

If either of these options are set, you can open the notebook `find_path.ipynb` and run the *first* block to see where the singularities of your families lie \(blue points\) and the path of integration \(red lines\).

If you want to edit the integration paths manually, open the file `current.sage` and change the end points of your rectilinear paths by setting `paths`.

### Straight deformation

The following will turn off finding a randomized, monomial deformation path to target hypersurface and will go straight to `f` from the best possible Fermat type hypersurface.

    Magma> PeriodHomotopy([f]:straight:=true);

If you want to go straight from a hand picked Fermat type hypersurface `g` to `f` then type: 

    Magma> ph := PeriodHomotopy([g,f]:straight:=true);

This option is not recommended for typical usage.

### Custom deformation path

If you have a sequence of smooth degree `d` polynomials `f0,f1,..,fr` where `f0` is of Fermat type, then you can type in this sequence into `PeriodHomotopy`. For example, with a sequence of 5 elements:

    Magma> PeriodHomotopy([f0,f1,f2,f3,f4,f5]);

Then the period homotopy will be performed on this sequence with `f5` being the target hypersurface.

### Non-random deformation path

Set the option `randomizePath` to `false`.

    Magma> PeriodHomotopy([f]:randomizePath:=false);

### Precision

If you want to increase or decrease the precision to which the integration is to be performed, set the `precision` option to the desired *number of digits*. By default we try to compute 100 digits of accuracy.

    Magma> PeriodHomotopy([f]:precision:=500);

### Classical periods

If for some reason you want the entire period matrix of the target hypersurface and not just the rows necessary to identify the Hodge decomposition, you can turn off the option `classicalPeriods`.

    Magma> PeriodHomotopy([f]:classicalPeriods:=false);
