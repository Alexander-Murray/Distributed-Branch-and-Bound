% Takes as input a cell-array of CasADi functions, g_fun, 
% and a vector, vars, where entry i is the number of input variables for
% the corresponding g_fun{i}

function g_glob_fun = concatCasadiFuns(g_fun,rvars,ivars)
    import casadi.*
    g_glob = [];
    for i = 1:length(g_fun)
        rvar{i} = SX.sym(strcat('x_',num2str(i)),[rvars(i),1]);
        ivar{i} = SX.sym(strcat('z_',num2str(i)),[ivars(i),1]);
        g_glob = [g_glob;g_fun{i}([rvar{i};ivar{i}])];
    end
    g_glob_fun = Function('f',{[vertcat(rvar{:});vertcat(ivar{:})]},{g_glob});
end