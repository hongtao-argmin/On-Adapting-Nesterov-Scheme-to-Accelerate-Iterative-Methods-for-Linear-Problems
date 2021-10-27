# On Adapting Nesterov's Scheme to Accelerate Iterative Methods for Linear Problems

[Tao Hong](https://hongtao.cswp.cs.technion.ac.il) and [Irad Yavneh](https://irad.cs.technion.ac.il), ``[On Adapting Nesterov’s Scheme to Accelerate Iterative Methods for Linear Problem](https://onlinelibrary.wiley.com/doi/full/10.1002/nla.2417)'', Accepted to Numerical Linear Algebra with Applications, 2021.

In this paper, we propose a closed-form solution to decide the optimal parameter inside Nesterov's scheme when the eigenvalues of the iteration matrix are real and the smallest and largest eigenvalues are given. Moreover, we show a sufficient condition which explicitly depicts a complex domain that the optimal parameter obtained through the smallest and largest real eigenvalues are still optimal when the iteration matrix has complex eigenvlaues. 

ArXiv Paper Link: https://arxiv.org/abs/2102.09239.

``Demo_Laplacian.m'': the demo of solving the poisson problem.

``Fig6_DampingFactorJacobi.m'': reproduce Fig. 6 of our paper.

``Figs4and5_ChebyVSNesterACFComplexEigenvalues'': reproduce Fig.s 4 and 5 of our paper.

Please add ``utilities'' folder in your own matlab path before running the codes. 

Feel free to shoot me (Tao Hong) an email: <hongtao@cs.technion.ac.il> if you find any bugs in this software or have any question regarding our paper.   
