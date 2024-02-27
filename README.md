# Attempt to do smoothing of calibration solutions that retains solution integrity

## Folder contents

- `aocal.py`: Stolen from [https://github.com/johnsmorgan/aocalpy](https://github.com/johnsmorgan/aocalpy) (git commit 7cc2d0e)
- `SB38969.B1934-638.beam0.aocalibrate.bin`: A test data set

## Regularisation

Regularisation by smoothing operates by weighting a "goodness-of-fit" cost function against a model "smoothness" cost function,
```math
Q(\hat{y}) = Q_\text{fit}(y, \hat{y}) + \lambda Q_\text{smooth}(\hat{y}),
```
where $`\lambda`$ is the so-called "regularisation parameter" that lets you decide whether the goodness-of-fit is more or less important relative to model smoothness, $`y`$ is the data, and $`\hat{y}`$ is the smoothed model we are trying to find.

For "simple" datasets (e.g. 1D functions), both cost functions can be defined in such a way that the problem of finding the model $`\hat{y}(x)`$ that minimises $`Q(\hat{y})`$ ls a linear problem and can be solved by standard least-squares techniques (e.g. see [Stickel, 2010](https://www.sciencedirect.com/science/article/pii/S0098135409002567)).
However, like any optimisation problem, it can also be done numerically, which may be necessary for more complicated cost functions.

## Calibration data sets

In this case, the data to which we are trying to fit smoothed models are calibration solutions, and the dimension we are smoothing over is frequency, represented hereafter by the subscript $`f`$.
Each data point, "$`y`$", is therefore a *set* of complex-valued Jones matrices of the form
```math
{\bf J} = \begin{bmatrix} j_{xx} & j_{xy} \\ j_{yx} & j_{yy} \end{bmatrix},
```
where a "set" consists of one Jones matrix per antenna element.
That is, for the $`f`$th frequency, we can write
```math
y_f = \left\{{\bf J}_{f,i} \,\middle|\, i \text{ is an antenna index}\right\}
```
Naturally, our model point, $`\hat{y}_f`$, will have the same form.

We must now define **goodness-of-fit** and **smoothness** cost functions.

## Goodness-of-fit

The data, $`y_f`$, are themselves the result of an optimisation problem whose cost function includes information (sky model, visibilities) that we will assume are no longer available to us.
Any goodness-of-fit cost function that is defined *without* that information is always at risk of producing a model solution, $`\hat{y}_f`$, that ultimately performs worse than the original solution.
It *may* perform better, but only if the "noise" component of the fits is due to some incorrect assumption (i.e. an inaccurate sky model or a badly modelled primary beam) which are otherwise expected to vary relatively smoothly over frequency.

The Jones matrices are "coupled" in the sense that each pair of antennas, i.e. a "baseline", form a sky measurement (for each frequency, so the $`f`$ subscript is suppressed here)
```math
{\bf E} = {\bf B}^{-1} {\bf J}_i^{-1} {\bf V} ({\bf J}_j^H)^{-1} ({\bf B}^H)^{-1},
```
where $`i`$ and $`j`$ are both antenna indices, $`\cdot^H`$ is the Hermitian operation, $`{\bf V}`$ is the masured visibility (a Hermitian matrix), $`{\bf B}`$ is the primary beam response, assumed here to be identical for all antennas, and 
```math
{\bf E} = \frac{1}{2} \begin{bmatrix} I + Q & U - Vi \\ U + Vi & I - Q \end{bmatrix}
```
is the measured sky represented as a coherency matrix.

We want to know the degree to which the set of $`{\bf J}_i`$'s can be altered without affecting the final measured sky.
Although there are some degrees of freedom hat exist under certain conditions (an unpolarised sky, and a "well-behaved" beam model), here we will make the least number of assumptions possible.
In this case, the only freedom is that, taken altogether, the set $`y_i`$ is equivalent up to multiplication by a unit phasor, $`e^{i\theta}`$, for arbitrary $`\theta`$.
The reason is that any scalar factor $`e^{i\theta}`$ that enters $`{\bf E}`$ via its presence in $`{\bf J}_i`$ cancels with its conjugate, $`e^{-i\theta}`$ present in $`{\bf J}_j^H`$.
(This fact is what makes possible the "trick" of dividing the Jones matrix elements by a reference antenna in order to better evaluate the frequency structure of the solutions.)
