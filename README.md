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

## Cost functions

In this case, the models we are trying to fit are for calibration solutions which are themselves the converged solutions to a previous numerical optimisation problem.
Each data point, "$`y`$", is in fact a complex-valued Jones matrix,
```math
J = \begin{bmatrix} j_{xx} & j_{xy} \\ j_{yx} & j_{yy} \end{bmatrix},
```

