# period-suite

Software to compute periods of hypersurfaces, with additional functions to make use of computed periods.

This is an implementation of the algorithm introduced in [Computing periods of hypersurfaces](https://sertozemre.files.wordpress.com/2018/03/computing_periods.pdf).

## Help

If you need help using this software, feel free to contact me at emresertoz \[at\] gmail \[dot\] com.  

Bug reports are most welcome.

## Installation

Two parts are needed to compute periods. One part, provided here, has everything except for the numerical integrator. The second part, the inegrator, has to be obtained from Marc Mezarrabbo's webpage.

#### Remarks:

- This software will only work on Unix based systems, e.g., Linux and Mac. 
- If you do not have the text processor `awk` already, please install it first.

### Part 1 - Picard-Fuchs equations

Clone this repository to a directory.

    git clone https://github.com/period-suite/period-suite.git /path/to/period-suite

From within the directory, open `suite.mag` and set `pathToSuite` to the *complete* path of this directory. 

Now open `integrator.ipynb` and `find_path.ipynb` to set the value of `pathToSuite` as before. For instance:

    sage --notebook=jupyter integrator.ipynb

    sage> pathToSuite="/path/to/period-suite";

### Part 2 - Integration

Navigate to [Mezzarobba's homepage](http://marc.mezzarobba.net/code/ore_algebra-analytic/). Read the instructions to install the *latest developmental* branch of `ore_algebra-analytic`.  His code is available on his private repository:

> http://marc.mezzarobba.net/code/ore_algebra-analytic.git/

### Recommendations

The newest, most stable versions of all software used are recommended and often mandatory.

- If you have a version of `ore_algebra-analytic` older than March 2017, please renew it for much better performance.
- Make sure you have the latest version of SageMath \(v8.1 or higher\) and Magma \(v2.23.8 or higher\). 
- If you do not use `PathFinder2` then Magma v.2.23.1 or higher will be sufficient.

******

# Usage

Navigate to the directory in which you have cloned this repository. Open a Magma session and load the file `suite.mag`.

    Magma> load "suite.mag";

Now you are set and can start using the code. We automatically define a few rings to start using the function but also to manipulate the output of the main function `PeriodHomotopy`. In particular `P` is already defined to be the polynomial ring with four variables `x,y,z,w`. These definitions are made in `suite.mag`; feel free to remove them if you wish.

    Magma> f:=-4*x^4 + 7*x*y^3 + 2*z^3*w + 6*z^2*w^2 + 6*w^4;
    Magma> ph:=PeriodHomotopy([f]);

Note the extra square brackets around `f`! See advanced usage below on how to give sequences to `PeriodHomotopy`.  The function `PeriodHomotopy` will produce the necessary ODEs and the initial conditions for the numerical integrator. 

Now, from another terminal, open the integrator in a SageMath jupyter notebook session:

    sage --notebook=jupyter integrator.ipynb

When you hit enter on the first code block, integration will start automatically. You can follow the progress of integration from the printed messages. Once integration is complete, the periods are printed on the screen and are written into the file `lastPeriods` in a way readable by Magma. 

From within SageMath, you can use the value `periods` or from within Magma type `load "lastPeriods"`. The latter command sets `M` to be the period matrix.

## Examples

See `examples.mag` for some examples.

## Important tip

By default, the deformation path to the target hypersurface is created randomly according to simple heuristics. If a computation takes more than a few seconds, run `PeriodHomotopy` several times and see if you can discover a deformation path for which the function terminates quickly.

This is very important! The result of a computation may take hours or a fraction of a second, depending on how good the deformation path is. For harder computations, it is worth spending time on the discovery of a good path.

******

# Advanced usage

There are many options and additional functions. We list the most important ones here. Recall the following definition.

A polynomial of the form `c0\*x0^d + ... + cn\*xn^d` where `c0,..,cn` are non-zero, is called of *Fermat type*.

### Output

Let `ph` be the output of `PeriodHomotopy` as above. It is a list with 3 elements:

```
ph[1] contains the deformation path used in period homotopy. 
ph[2] contains the polynomial part of the forms whose periods are to be computed, one list per family.
ph[3] contains the ODEs for the periods of the forms in ph[2].
```

But the main bulk of the output, including the initial conditions and transition matricies between families, is written in the file `current.sage`, where it is read by `integrator.ipynb`.

### Path finder

The integrator needs a path in the complex plane for every one of the systems it integrates. By default it will try the straigh path 0 to 1, it works most of the time. However, if there are singularities along the way the computations might not terminate. In this case we recommend setting one of the following options on `PeriodHomotopy`:

    Magma> ph := PeriodHomotopy([f]: pathfinder1:=true);

or 

    Magma> ph := PeriodHomotopy([f]: pathfinder2:=true);

Both have advantages and disadvanteges. But we recommend that you use the second one only when the first one doesn't terminate.

If either of these options are set, you can open the notebook `find_path.ipynb` and run the *first* block to see where the singularities of your families lie \(blue points\) and the path of integration \(red lines\).

If you want to edit the integration paths manually, open the file `current.sage` and change the end points of your rectilinear paths.

### Straight deformation

The following will turn off finding a randomized, monomial deformation path to target hypersurface and will go straight to `f` from the best possible Fermat type hypersurface.

    Magma> ph := PeriodHomotopy([f]:straight:=true);

If you want to go straight from a hand picked Fermat type hypersurface `g` to `f` then type: 

    Magma> ph := PeriodHomotopy([g,f]:straight:=true);

This option is not recommended for typical usage.

### Custom deformation path

If you have a sequence of smooth degree `d` polynomials `f0,f1,..,fr` where `f0` is of Fermat type, then you can type in this sequence into `PeriodHomotopy`. For example, with a sequence of 5 elements:

    Magma> ph:=PeriodHomotopy([f0,f1,f2,f3,f4,f5]);

Then the period homotopy will be performed on this sequence with `f5` the target hypersurface.

### Non-random deformation path

Set the option `randomizePath` to `false`.

    Magma> ph:=PeriodHomotopy([f]:randomizePath:=false);

### Precision

If you want to increase or decrease the precision to which initial conditions are computed, set the `precision` option to the desired *number of digits*. By default we compute 200 digits.

    Magma> ph:=PeriodHomotopy([f]:precision:=500);

### Classical periods

If for some reason you want the entire period matrix of the target hypersurface and not just the rows necessary to identify the Hodge decomposition, you can turn off the option `classicalPeriods`.

    Magma> ph:=PeriodHomotopy([f]:classicalPeriods:=false);
