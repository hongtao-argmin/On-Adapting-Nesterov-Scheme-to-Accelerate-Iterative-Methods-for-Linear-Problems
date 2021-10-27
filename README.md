# On Adapting Nesterov's Scheme to Accelerate Iterative Methods for Linear Problems -- To appear in Numerical Linear Algebra with Applications

In this paper, we propose a closed-form solution to decide the optimal parameter inside Nesterov's scheme when the eigenvalues of the corresponding iteration matrix are real and the smallest and largest eigenvalues are given. Moreover, we show a sufficient condition which explicitly depicts a complex domain that the optimal parameter obtained through the smallest and largest real eigenvalues are still valid when the iteration matrix has complex eigenvlaues.


Paper Link: https://arxiv.org/abs/2102.09239

To run these codes, please download all of the codes and the corresponding necessary external toolboxes. Note that do not forget to add their path in your own matlab path. 

Two toolboxes named omp and ksvd are needed to run this toolbox which can be download from http://www.cs.technion.ac.il/~ronrubin/software.html.

The optimization toolbox named minFunc can be found in https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html.

The online algorithm shown in 
J. Mairal, F. Bach, J. Ponce, and G. Sapiro, ``Online dictionary learning for sparse coding,’’  Proceedings of the 26th annual international conference on machine learning. ACM, 2009 can be download from http://spams-devel.gforge.inria.fr/.  
 
In our case, we also redo this algorithm which is named ` Online_DIC_MBPS09.m’.  

The file named ‘Experiment A.m’, ‘Experiment  B.m’and ‘Experiment C.m’ are used to redo the corresponding experiments mentioned in the paper. 
The main algorithm in this paper is named ‘Robust_project_matrix.m’.

Remark:
The traindata and testdata should be stored as X and X_test.

To obtain how to utilize the related functions in this toolbox, just use `help function_name’.


If you find any bugs in this toolbox, feel free to contact the author Tao Hong through email: hongtao@cs.technion.ac.il.   

Implementation of  ``On Adapting Nesterov’s Scheme to Accelerate Iterative Methods for Linear Problem''
