function [ sol, stats] = solveMINLP( f,ineq,lbx,ubx,x0, discrete )
    import casadi.*
    nx= length(x0);
    z = SX.sym('z',nx,1);

    fcas    = f(z);
    
    ineqcas = ineq(z); 
    
    gdim = length(ineq);
    
    opts.print_time = 0;
    opts.verbose = 0; 
    opts.discrete = discrete;
    
    nlp = struct('x',z,'f',fcas,'g',ineqcas);
    
    try
        cas = qpsol('solver','gurobi',nlp, opts);
    catch
        opts.verbose_init = 0;
        opts.bonmin.print_level = 0;
        cas = nlpsol('solver','bonmin',nlp, opts);
    end
    sol = cas('x0' , x0,...
            'lbx', lbx,...
            'ubx', ubx,...
            'lbg', -inf*ones(gdim,1), ...
            'ubg', zeros(gdim,1));
   
    stats = cas.stats.return_status;
    
end
