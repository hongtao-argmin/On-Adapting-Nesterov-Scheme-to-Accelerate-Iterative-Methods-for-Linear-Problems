# On Adapting Nesterov's Scheme to Accelerate Iterative Methods for Linear Problems

Tao Hong and Irad Yavneh, ``On Adapting Nesterov’s Scheme to Accelerate Iterative Methods for Linear Problem'', Accepted to Numerical Linear Algebra with Applications, 2021.

In this paper, we propose a closed-form solution to decide the optimal parameter inside Nesterov's scheme when the eigenvalues of the iteration matrix are real and the smallest and largest eigenvalues are given. Moreover, we show a sufficient condition which explicitly depicts a complex domain that the optimal parameter obtained through the smallest and largest real eigenvalues are still optimal when the iteration matrix has complex eigenvlaues. 

Paper Link: https://arxiv.org/abs/2102.09239.

``Demo_Laplacian.m'': the demo of solving the Laplacian problem.
``Fig6_DampingFactorJacobi.m'': reproduce Fig. 6 of our paper.
``Figs4and5_ChebyVSNesterACFComplexEigenvalues'': reproduce Fig.s 4 and 5 of our paper.

Please add ``utilities'' folder in your own matlab path before running the codes. 

If you find any bugs in this toolbox, feel free to contact the author Tao Hong through email: hongtao@cs.technion.ac.il.   
