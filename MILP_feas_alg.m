% Finds a feasible solution to Ax+Bz=C 
% s.t. lbx <= x <= ubx and lbz <= z <= ubz
%
% To cite use:
% @Article{Murray_2020d,
%   author  = {Murray, A. and Hagenmeyer, V.},
%   title   = {On Convergence of Mixed-Integer ALADIN},
%   journal = {4th International Conference on Algorithms, Computing and Systems},
%   year    = {2020},
% }

function [x_sol,z_sol, flag] = MILP_feas_alg(rA,A,C,lbx,ubx,lbz,ubz)
flag = 0;

rvars = length(lbx);
ivars = length(lbz);

    if or(sum(lbz>ubz)>0, sum(lbx>ubx)>0) % infeasible bounds
        finish = 1;
    elseif rA*lbx+A*lbz==C % lower bounds satisfy consensus constraints
        finish = 1;
        flag = 1;
        x_sol = lbx;
        z_sol = lbz;
    elseif rA*ubx+A*ubz==C % upper bounds satisfy consensus constraints
        finish = 1;
        flag = 1;
        x_sol = ubx;
        z_sol = ubz;
    else % the problem isn't so easy...
        % initial guess
        x_sol = lbx; 
        z_sol = lbz;
        
        % full consensus matrix
        B = [rA,A]; 
        
        % use gaussian elimination to reduce the size of the system
        n = size(B,2);
        reduced_sys = rref([B,C]);
        reduced_sys = remove_rows_of_zeros(reduced_sys);
        B = reduced_sys(:,1:n);
        C = reduced_sys(:,n+1);
        
        %eliminate empty columns
        empty_indices = [];
        for i = 1:n
           if sum(abs(B(:,i)))==0 
              empty_indices = [empty_indices,i]; 
           end
        end
        B(:,empty_indices)=[];
        lb = [lbx;lbz];
        ub = [ubx;ubz];
        lb(empty_indices) = [];
        ub(empty_indices) = [];
        
        %put identity submatrix first
        [B,new_indices] = rrefFixer(B);

        %parameterize!
        r = rank(B);
        N = size(B,1);
        n = size(B,2);
        v0 = [C;zeros(n-N,1)];
        Vs = [-B(1:r,r+1:n);diag(ones(n-r,1))];

        % init
        z_star = v0; % some point in the hyperplane
        dist = distance(z_star,lb,ub); % distance to box constraints
        if sum(z_star>ub)+sum(z_star<lb)==0 % check if inside of box constraints
           flag = 1;
           finish = 1;
        else
            finish = 0;
        end
    end
    
    if finish ~=1
        % main loop
        cycle_breaker = []; % this will be used to prevent iterations from repeating
        dim = length(z_star); % length of solution vector
        fail = 0;
        while and(finish == 0, fail>=dim)
            val1 = zeros(n-r,1); val2 = zeros(n-r,1);
            for i = 1:n-r
               test1 = z_star + Vs(:,i); % candidate vector
               test2 = z_star - Vs(:,i); % candidate vector
 
               val1(i) = distance(test1,lb,ub); %how much closer do we get by applying Vs(:,i)?
               val2(i) = distance(test2,lb,ub); %how much closer do we get by applying -Vs(:,i)?
            end
            update_candidates = setdiff(find(min([val1;val2])==[val1;val2]),cycle_breaker); % prevent cyclically choosing bases
            if isempty(update_candidates)
               finish = 1; 
            else
                update_ind = update_candidates(1); % just in case multiple solutions are available
                if update_ind <=n-r
                    z_star = z_star + Vs(:,update_ind);
                else
                    z_star = z_star - Vs(:,update_ind-n+r);
                end
            end
            % check whether we are actually closer. if not, then the feasible set does not include the hyperplane
            new_dist = distance(z_star,lb,ub);
            if new_dist == 0 % current iteration satisfies constraints
                flag = 1;
                finish = 1;
                fail = 0;
            elseif new_dist>dist % we have moved further from the box constraints
               finish = 1;
               fail = fail + 1;
            elseif new_dist == dist % we haven't moved closer to the box constraints
                cycle_breaker = [cycle_breaker,update_ind];
                fail = fail + 1;
            else
                dist = new_dist;
                cycle_breaker = [];
                fail = 0;
            end
        end

        % output solution
        reindexed_zstar = sortrows([new_indices',z_star]);
        sol = reindexed_zstar(:,2);
        nonempty_indices = setdiff(1:rvars+ivars,empty_indices);
        x_sol(setdiff(1:rvars,empty_indices)) = sol(find(nonempty_indices<=rvars)); 
        z_sol(setdiff(1:ivars,empty_indices(find(empty_indices>rvars))-rvars)) = sol(find(nonempty_indices>rvars));
    end
end