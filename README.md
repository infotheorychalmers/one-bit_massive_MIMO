# One-bit Massive MIMO uplink  
=============================

Numerical routines for the computation of the achievable uplink throughput with PSK/QAM for a one-bit quantized massive MIMO system with a MRC/ZF receiver. The system operates over a Rayleigh fading channel with no a priori CSI available at the transmitters and the receiver. Hence, the channel fading coefficients has to be estimated based on coarsely quantized data.

Please cite the following paper when using the code:
  * S. Jacobsson, G. Durisi, M. Coldrey, U. Gustavsson, and C. Studer, “One-bit massive MIMO: channel estimation and high-order modulations,” in Proc. IEEE Int. Conf. Commun. (ICC), London, U.K., Jun. 2015, to appear.

Dependencies
------------
The code is written in MATLAB. The following additional MATLAB routines (available under the BSD license) are required:

 * mtimesx (http://www.mathworks.com/matlabcentral/fileexchange/25977-mtimesx-fast-matrix-multiply-with-multi-dimensional-support),
 * multinv (http://www.mathworks.com/matlabcentral/fileexchange/31222-inversion-every-2d-slice-for-arbitrary-multi-dimension-array),
 * histcn (http://www.mathworks.com/matlabcentral/fileexchange/23897-n-dimensional-histogram).

Running
-------

In order to run the simulations. The above mentioned routines (provided in this release) have to be placed in the same directory. Note that mtimesx routine requires an installed C-compiler on your system, and has to be compiled before being used.

The following MATLAB routines are provided in the repository

  * simRate.m,
  * pilotOptSoft.m,
  * plotRateVsN.m,
  * plotRateVsSNR.m,
  * plotRateVsT.m.

Details on how to use the routines can be found in each file.
