# cookflip_code

This repository contains Matlab code for the paper "[The Mathematics of
Burger Flipping][3]" by [Jean-Luc Thiffeault][1].

## Functions

The following MATLAB functions are part of this repository:

* heat        - Temperature profile for the heat equation.
* heatsteady  - Steady temperature profile for heat equation.
* heateigval  - Eigenvalues (spectrum) of the heat operator.
* heateigfun  - Eigenfunctions of the heat operator.
* tcookthru   - The time to cook food through without flipping.
* flipheatfix - The fixed point of the flip-heat operation for fixed time.
* tcooksym    - Time to cook for equal flips and symmetric boundary conditions.
* cooktime    - Cooking time on a hot plate for several flips of the food.
* mincooktime - Minimize cooking time for given number of flips of the food.
* flipheatbl  - Boundary layer for the flip-heat operator for rapid flip time.
* flipheatop  - The flip-heat operator.
* flipop      - The "flipping" operator acting on eigenfunctions.

## Citing this work

Thiffeault, J.-L. "The Mathematics of burger flipping," _Physica D:
Nonlinear Phenomena_ **439**, 133410 (2022).
DOI: [10.1016/j.physd.2022.133410][4]

BibTeX entry:
```latex
@Article{Thiffeault2022,
  title =    {The mathematics of burger flipping},
  journal =  {Physica D: Nonlinear Phenomena},
  volume =   439,
  pages =    133410,
  year =     2022,
  issn =     {0167-2789},
  doi =      {10.1016/j.physd.2022.133410},
  url =      {https://doi.org/10.1016/j.physd.2022.133410},
  author =   {Jean-Luc Thiffeault},
  keywords = {Heat equation, Optimization, Cooking},
  abstract = {What is the most effective way to grill food? Timing
              is everything, since only one surface is exposed to
              heat at a given time. Should we flip only once, or
              many times? We present a simple model of cooking by
              flipping, and some interesting observations
              emerge. The rate of cooking depends on the spectrum
              of a linear operator, and on the fixed point of a
              map. If the system has symmetric thermal properties,
              the rate of cooking becomes independent of the
              sequence of flips, as long as the last point to be
              cooked is the midpoint. After numerical
              optimization, the flipping intervals become roughly
              equal in duration as their number is increased,
              though the final interval is significantly
              longer. We find that the optimal improvement in
              cooking time, given an arbitrary number of flips, is
              about 29\% over a single flip. This toy problem has
              some characteristics reminiscent of turbulent
              thermal convection, such as a uniform average
              interior temperature with boundary layers.}
}
```

## License

This code is released under the MIT License.  See the file
[LICENSE][2] for copying permission.

[1]: http://www.math.wisc.edu/~jeanluc/
[2]: http://github.com/jeanluct/cookflip_code/raw/main/LICENSE
[3]: https://arxiv.org/abs/2206.13900
[4]: https://doi.org/10.1016/j.physd.2022.133410
