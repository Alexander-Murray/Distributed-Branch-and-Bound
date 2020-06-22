% Given a mixed integer program with constraints
% g(x,z)<0
% Ax + Bz = c
% lbx< x <ubx
% lbz< z <ubz
% Then return a boolean value feas = 1 if the point (x0,z0) satisfies the
% above constraints and feas = 0 if not

function feas = check_feas_MIP(eps,g,A,B,c,lbx,lbz,ubx,ubz,x0,z0)
    
    feas = 0;
    if sum(abs(A*x0 + B*z0 - c))<eps
        ineqi = length(full(g([x0;z0])));
        if sum(full(g([x0;z0]))<=eps*ones(ineqi,1)) == ineqi
            if and(sum(x0<=ubx)==length(x0), sum(x0>=lbx)==length(x0)) %only for one var per subproblem
                if and(sum(z0<=ubz)==length(z0), sum(z0>= lbz)==length(z0))
                    feas = 1;               
                end
            end
        end
    end
end