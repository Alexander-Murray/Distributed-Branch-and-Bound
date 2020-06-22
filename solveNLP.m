function [ sol, log] = solveNLP( f,ineq,lbx,ubx,x0,opts)
    import casadi.*
    nx= length(x0);
    z = SX.sym('z',nx,1);

    fcas    = f(z);
    
    ineqcas = ineq(z); 
    
    gdim = length(ineqcas);

    nlp = struct('x',z,'f',fcas,'g',ineqcas);
    
    tic
    if strcmp(opts.NLP_solver,'ipopt') 
        cas = nlpsol('solver','ipopt',nlp, opts.NLPsolver_opts);
        sol = cas('x0' , x0,...
                'lbx', lbx,...
                'ubx', ubx,...
                'lbg', -inf(gdim,1), ...
                'ubg', zeros(gdim,1));
    else
        error('No other NLP solvers currently implemented. Use IPOPT!')
    end
    log.time = toc;
    log.status = cas.stats.return_status;
    log.obj = full(sol.f);
    log.lambda = full(sol.lam_g);
    
end
