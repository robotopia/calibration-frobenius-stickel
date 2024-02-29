# Equivalence classes of Jones matrices

Let:
- $`{\bf V}`$ be a set of measured visibilities (a Hermitian $`2 \times 2`$ matrix) for baseline $`ij`$;
- $`{\bf J}_i`$ be the Jones matrix representing the $`i`$th antenna element's response;
- $`{\bf B}`$ be the primary beam response (assumed identical for all antenna elements); and
- $`{\bf E}`$ be the coherency matrix of the sky signal:
```math
{\bf E} \equiv \frac{1}{2} \begin{bmatrix} I + Q & U - Vi \\ U + Vi & I - Q \end{bmatrix}
```

The relationship between the sky, $`{\bf E}`$, and the measured visibilities, $`{\bf V}`$, is:
```math
{\bf E} = {\bf B}^{-1} {\bf J}_i^{-1} {\bf V_{ij}} ({\bf J}_j^H)^{-1} ({\bf B}^H)^{-1},
```

An interesting question is: under what conditions do equivalence classes of instrumental Jones matrices, $`{\bf J}_i`$, exist, in the sense that swapping out a particular Jones matrix with another member of its class doesn't change the predicted sky signal?

## Stokes I skies

When $`{\bf E}`$ is purely unpolarised, it has the form
```math
{\bf E} = \frac{1}{2} \begin{bmatrix} I & 0 \\ 0 & I \end{bmatrix}
        = \frac{I}{2} \begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix}
```
Since this just a scaled identity matrix, it will commute with any other matrix.
Thus, if we multiply $`{\bf E}`$ on the left by any arbitrary matrix, $`{\bf M}`$, we can ensure that the LHS remains unchanged by multiplying it on the right by $`{\bf M}^{-1}`$.

For an individual baseline, it will always be possible (in the manner shown below) to adjust $`{\bf J}_i`$ and $`{\bf J}_j`$ in such a way to absorb (respectively) the effects of $`{\bf M}`$ and $`{\bf M}^{-1}`$.
However, we are interested in finding just the subset of matrices, $`{\bf M}`$, for which the same transformation can be applied to all $`{\bf J}`$'s.
