# cookflip_code

This repository contains Matlab code for the paper "The Mathematics of
Burger Flipping" by [Jean-Luc Thiffeault][1].

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

## License

This code is released under the MIT License.  See the file
[LICENSE][2] for copying permission.

[1]: http://www.math.wisc.edu/~jeanluc/
[2]: http://github.com/jeanluct/cookflip_code/raw/main/LICENSE
