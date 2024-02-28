# Attempt to do smoothing of calibration solutions that retains solution integrity

## Folder contents

- `aocal.py`: Stolen from [https://github.com/johnsmorgan/aocalpy](https://github.com/johnsmorgan/aocalpy) (git commit 7cc2d0e)
- `SB38969.B1934-638.beam0.aocalibrate.bin`: A test data set

## Regularisation

Regularisation by smoothing operates by weighting a "goodness-of-fit" cost function against a model "smoothness" cost function,
```math
C(\hat{y}) = C_\text{fit}(y, \hat{y}) + \lambda C_\text{smooth}(\hat{y}),
```
where $`\lambda`$ is the so-called "regularisation parameter" that lets you decide whether the goodness-of-fit is more or less important relative to model smoothness, $`y`$ is the data, and $`\hat{y}`$ is the smoothed model we are trying to find.

For "simple" datasets (e.g. 1D functions), both cost functions can be defined in such a way that the problem of finding the model $`\hat{y}(x)`$ that minimises $`C(\hat{y})`$ ls a linear problem and can be solved by standard least-squares techniques (e.g. see [Stickel, 2010](https://www.sciencedirect.com/science/article/pii/S0098135409002567)).
However, like any optimisation problem, it can also be done numerically, which may be necessary for more complicated cost functions.

## Calibration data sets

In this case, the data to which we are trying to fit smoothed models are calibration solutions, and the dimension we are smoothing over is frequency, represented hereafter by the subscript $`f`$.
Each data point, "$`y`$", is therefore a *set* of complex-valued Jones matrices of the form
```math
{\bf J} = \begin{bmatrix} j_{XX} & j_{XY} \\ j_{YX} & j_{YY} \end{bmatrix},
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
{\bf E} \equiv \frac{1}{2} \begin{bmatrix} I + Q & U - Vi \\ U + Vi & I - Q \end{bmatrix}
```
is the measured sky represented as a coherency matrix.

We want to know the degree to which the set of $`{\bf J}_i`$'s can be altered without affecting the final measured sky.
Although there are some degrees of freedom hat exist under certain conditions (an unpolarised sky, and a "well-behaved" beam model), here we will make the least number of assumptions possible.
In this case, the only freedom is that, taken altogether, the set $`y_i`$ is equivalent up to multiplication by a unit phasor, $`e^{i\theta}`$, for arbitrary $`\theta`$.
The reason is that any scalar factor $`e^{i\theta}`$ that enters $`{\bf E}`$ via its presence in $`{\bf J}_i`$ cancels with its conjugate, $`e^{-i\theta}`$ present in $`{\bf J}_j^H`$.
This fact is what makes possible the "trick" of dividing the Jones matrix elements by a reference antenna in order to better evaluate the frequency structure of the solutions.

Rather than relying on the use of an arbitrary reference antenna, we would like a goodness-of-fit that can "absorb" this factor in its very definition.
Suppose we define a residual Jones matrix that depends on a choice of $`\theta`$ (still suppressing the subscript $`f`$ for now):
```math
{\bf R}_{i,\theta} \equiv y_i - \hat{y}_i = \left\{{\bf J}_i - e^{i\theta} \hat{\bf J}_i \middle| i \text{ is an antenna index}\right\},
```
where care must be taken to distinguish the antenna index $`i`$ in the subscripts from the imaginary number $`i`$ in the exponent.
We use the Frobenius norm to convert the residual Jones into a real-valued objective function that can be minimised:
```math
C_{\text{fit},\theta} = \sum_i \lVert {\bf R}_{i,\theta} \rVert_F^2
```
The "best", or "global" goodness-of-fit is the minimum $`C_{\text{fit},\theta}`$ in the set:
```math
C_\text{fit} \equiv \min \left\{C_{\text{fit},\theta} \middle| \theta \in [0,2\pi)\right\}.
```

### Finding the optimal $`\theta`$

$`\theta`$ is *not* a free parameter to be fitted; the optimal value of $`\theta`$ which minimises the cost functionmust be found at each step of the numerical solution.
Fortunately, it can be solved for analytically, as shown below.

First, note that the use of the Frobenius norm to construct the objective function treats each element of the Jones matrices on the same footing.
That is, for the purposes of this smoothing exercise, the arrangement of the elements $`j_{XX}`$ into a $`2 \times 2`$ matrix is irrelevant; they could have equivalently been arranged in a single vector containing the elements from all antennas, for example.
The theta-dependent objective function could therefore be equivalently expressed as
```math
C_{\text{fit},\theta} = \sum_i \sum_{p\in\{X,Y\}} \sum_{q\in\{X,Y\}} |j_{pq,i} - e^{i\theta} \hat{j}_{pq,i}|^2.
```
Each individual term has an expansion
```math
|j_{pq,i} - e^{i\theta} \hat{j}_{pq,i}|^2
    = (j_{pq,i} - e^{i\theta} \hat{j}_{pq,i})(j_{pq,i}^\ast - e^{-i\theta} \hat{j}_{pq,i}^\ast)
    = |j_{pq,i}|^2 + |\hat{j}_{pq,i}|^2 - e^{i\theta} j_{pq,i}^\ast \hat{j}_{pq,i} - e^{-i\theta} j_{pq,i} \hat{j}_{pq,i}^\ast.
```
Therefore, minimising $`C_{\text{fit},\theta}`$ is equivalent to maximising
```math
\sum_{p,q,i} \left( e^{i\theta} j_{pq,i}^\ast \hat{j}_{pq,i} + e^{-i\theta} j_{pq,i} \hat{j}_{pq,i}^\ast \right)
    = \sum_{p,q,i} 2 \Re \left[ e^{i\theta} j_{pq,i}^\ast \hat{j}_{pq,i} \right]
    = 2 \Re \left[ e^{i\theta} \sum_{p,q,i} j_{pq,i}^\ast \hat{j}_{pq,i} \right],
```
since the sum of the real parts equals the real part of the sum, and since the common factor of $`e^{i\theta}`$ can be factored out.
This will be maximised when the sum itself is already purely real, which occures when $`\theta`$ is set to
```math
\theta_\text{max}
    = -\arg \left[ \sum_{p,q,i} j_{pq,i}^\ast \hat{j}_{pq,i} \right]
    = \arg \left[ \sum_{p,q,i} j_{pq,i} \hat{j}_{pq,i}^\ast \right].
```

The objective function is then simply
```math
C_\text{fit} = \sum_i \lVert {\bf R}_{i,\theta_\text{max}} \rVert_F^2.
```

## Smoothness

Smoothness refers to how sets Jones matrices compare across nearby frequency bins.
