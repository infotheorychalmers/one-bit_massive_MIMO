# One-bit Massive MIMO  
======================
Numerical routines for the computation of the achievable uplink throughput with PSK/QAM for a one-bit quantized massive MIMO system with a MRC/ZF receiver. The system operates over a Rayleigh fading channel with no a priori CSI available at the transmitters and the receiver. Hence, the channel fading coefficients has to be estimated based on coarsely quantized data.

Please cite the following paper when using the code:
  * S. Jacobsson, G. Durisi, M. Coldrey, U. Gustavsson, and C. Studer, “One-bit massive MIMO: channel estimation and high-order modulations,” in Proc. IEEE Int. Conf. Commun. (ICC), London, U.K., Jun. 2015, to appear.

Dependencies
------------
The code is written in MATLAB. The following additional MATLAB routines (available under the BSD license) are required:

 * MTIMESX (http://www.mathworks.com/matlabcentral/fileexchange/25977-mtimesx-fast-matrix-multiply-with-multi-dimensional-support),
 * MULTINV (http://www.mathworks.com/matlabcentral/fileexchange/31222-inversion-every-2d-slice-for-arbitrary-multi-dimension-array),
 * HISTCN (http://www.mathworks.com/matlabcentral/fileexchange/23897-n-dimensional-histogram).

Running
-------

The following MATLAB routines are provided in the repository

  * simRate.m,
  * pilotOptSoft.m,
  * plotRateVsN.m,
  * plotRateVsSNR.m,
  * plotRateVsT.m.

Details on how to use the routines can be found in each file.

Note that the mtimesx routines requires an installed C-compiler on your system, and has to be compiled before being used. Instructions on how to do some can be found in the MTIMESX manual (found in the mtimesx folder). However, two precompiled configurations that should be able to run on 64-bit Windows and OSX systems are provided in the release.


