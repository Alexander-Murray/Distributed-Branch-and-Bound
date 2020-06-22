%check if the current node is at the bottom of the branch and bound tree.
%If so, then return a value for a possible new upper bound

function U = get_U(dec_cur,x_sol,z_sol,f_fun)
    if isempty(dec_cur == 0)
        U = 0;
        for i = 1:length(f_fun)
           U = U + full(f_fun{i}([x_sol;z_sol])); 
        end
    else
        U = inf;
    end
end