function [Sig,lam0,dist_alg_opts] = DBB_subalg_setup(rvars,ivars,NC,opts)
    import casadi.*
    n = length(rvars);
    if strcmp(opts.NLP_solver,'ALADIN')
        for i = 1:n
            Sig{i} = diag(ones(rvars(i)+ivars(i),1));
        end
        lam0 = zeros(NC,1);
        dist_alg_opts.solveQP = 'backslash';
        dist_alg_opts.MINLPsolver_opts.print_time = 0;
        dist_alg_opts.rho0 = 100;
        dist_alg_opts.mu0 = 1000;
        dist_alg_opts.cent_sol = [];
        dist_alg_opts.maxiter = 200;
        dist_alg_opts.eps = 10^-3;
        dist_alg_opts.min_iter = 1;
        dist_alg_opts.rhoMax = 5*10^4;
        dist_alg_opts.muMax = 10^8;
        dist_alg_opts.rhoUpdate = 1.1;
        dist_alg_opts.muUpdate = 1.2;
    elseif strcmp(opts.NLP_solver,'ADMM')
        for i = 1:n
            Sig{i} = sparse(diag(ones(NC,1)));
            lam0{i} = 10^-5*ones(NC,1); %for ADMM sub-alg
        end
        dist_alg_opts.rho0 = 100;
        dist_alg_opts.maxiter = 100;
        dist_alg_opts.eps = 10^-3;
        dist_alg_opts.rhoUpdate = 1.1;
    elseif strcmp(opts.NLP_solver,'ADMM_exp_lasso')
        dist_alg_opts.QUIET    = opts.exp_lasso_opts.QUIET;
        dist_alg_opts.MAX_ITER = opts.exp_lasso_opts.MAX_ITER;
        dist_alg_opts.ABSTOL   = opts.exp_lasso_opts.ABSTOL;
        dist_alg_opts.RELTOL   = opts.exp_lasso_opts.RELTOL;
        Sig = [];
        lam0 = [];
    elseif strcmp(opts.NLP_solver,'IPOPT')
    else
        error('Invalid solver selected. Options are: ALADIN, ADMM, ADMM_exp_lasso, IPOPT')
    end
end