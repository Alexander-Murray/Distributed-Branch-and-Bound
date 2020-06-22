function [timeout,eq_relax,cont_relax,nav_strat_U,nav_strat_L,U,L,eps,maxiter] = BB_param_check(params)
    if isfield(params,'timeout')
        timeout = params.timeout;
    else
        timeout = 300;
    end
    if isfield(params,'eq_relax')
        eq_relax = params.eq_relax;
    else eq_relax = 0;
    end
    if isfield(params,'cont_relax')
        cont_relax = params.cont_relax;
    else cont_relax = 1;
    end
    if isfield(params,'nav_strat_U')
        nav_strat_U = params.nav_strat_U;
    else nav_strat_U = 'best-bounds';
    end
    if isfield(params,'nav_strat_L')
        nav_strat_L = params.nav_strat_L;
    else nav_strat_L = 'improve L';
    end
    if isfield(params,'U')
        U = params.U;
    else U = inf;
    end
    if isfield(params,'L')
        L = params.L;
    else L = -inf;
    end
    if isfield(params,'eps')
        eps = params.eps;
    else eps = 10^-3;
    end
    if isfield(params,'maxiter')
        maxiter = params.maxiter;
    else maxiter = 1000;
    end
end
