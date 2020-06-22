# Distributed Branch and Bound  

Branch and Bound algorithm with parallelized NLP subproblems. An ADMM algorithm is provided as the default method for distributed optimiztion of the NLP subproblems. Standard branch and bound method also included. This code uses the open-source tool CasADi:  
https://web.casadi.org/  

Input problem must be partitioned:  
**obj_funs:** cell array of CasADi functions  
**cons_funs:** cell array of CasADi functions  
**lbx:** cell array of vectors  
**ubx:** cell array of vectors  
**lbz:** cell array of vectors    
**ubz:** cell array of vectors  
**A:** cell array of matrices  
**B:** cell array of matrices  
**c:** vector  
**x0:** cell array of vectors  
**z0:** cell array of vectors  
**params:** struct  
**opts:** struct  

Params will assume default values if left empty. Possible params to set:  
**U:** double. User provided upper bound. Default: inf;   
**L:** double. User provided lower bound. Default: -inf;  
**eq_relax:** 1,0. Solve input problem with relaxed consensus constraints at root node. Default: 0   
**cont_relax:** 1,0. Solve continuous relaxtion of input problem at root node. Default: 1  
**nav_strat_U:** 'depth-first', 'breadth-first', 'best-first', 'best bounds', 'improve L'. Node navigation strategy until an upper bound has been established. Default: 'best-first'  
**nav_strat_L:** 'depth-first', 'breadth-first', 'best-first', 'best bounds', 'improve L'. Node navigation strategy after an upper bound has been established. Default: 'improve L'  
**eps:** double. Termination tolerance. Default: 10^-3;  
**timeout:** double. Time limit for algorithm. Default: 300  
**maxiter:** int. Iteration limit. Default: 1000;   

Opts must be set. Algorithm options are:  
**NLP_solver:** 'ipopt', 'ALADIN', 'ADMM', 'ADMM_exp_lasso'  
**NLPsolver_opts:** used to pass options to the NLP solver   
**custom_weights:** optional node weighting.   

The syntax for using the algorithm is shown in the EXAMPLE_MIP. The option 'ADMM_exp_lasso' calls on a specialized ADMM algorithm for a problem of the form shown in MI_lasso_problem.

To cite use:  
@Article{Murray_2020e,
  author  = {Murray, A. and Hagenmeyer, V.},
  title   = {A New Distributed Optimization Algorithm for MINLPs with Affinely Coupled Decision Variables},
  journal = {4th International Conference on Algorithms, Computing and Systems},
  year    = {2020},
}
