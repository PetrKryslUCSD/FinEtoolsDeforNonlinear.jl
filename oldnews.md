# Past news

- 12/19/2019: Fixed flawed scaling of the threaded calculation: Julia threads do not cooperate with BLAS threads. All the BLAS calls needed to be eliminated from the explicit code in order to obtain good parallel efficiency.
- 12/13/2019: Instrumented an example of transient (explicit) dynamics so that runs in parallel on multiple threads.
- 12/09/2019: Added an example of transient (explicit) dynamics.
- 10/12/2019: Corrected a design flaw in the matrix utilities module.
- 07/28/2019: Implemented automatic differentiation in material models.
- 06/11/2019: Applications are now separated  out from the `FinEtools` package.
