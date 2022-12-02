# Linear mixed-effects models workshop

This is an introductory 2 hours course on linear multilevel models using R and the `lme4` library.

See the slides online here [mlisi.xyz/RHUL-stats/workshops.html](https://mlisi.xyz/RHUL-stats/workshops.html).

### Outline

+ **Part 1**: introduction to LMM and `lme4` package
+ **Part 2**: 
  1. Estimation methods, ML and REML
  2. Parametric bootstrapping
  3. Power analyses via fake-data simulation
  4. Convergence warnings

The folder `exercises` contains some worked examples and exercises.

---

The slides are compiled using the `xaringan` library. If you want to compile them locally, you may have first to run some code bits (in particular `code_bits_part2.R`) for the second part of the course, as these generate comutationally demanding outpurs such as bootstrap distributions of simulations. (you may have to change the path to that of your machine)
