%%% NOT COMPATIBLE WITH R2017b OR OLDER! %%% (issue with conversion between strings and integer arrays)
function [XZsol,log] = Dist_BnB(f_fun,g_fun,lbx,lbz,ubx,ubz,A,B,c,x0,z0,ints,int_parts,params,opts)
 
good_statuses = {'Maximum_Iterations_Exceeded','Restoration_Failed','Solve_Succeeded','SUCCESS'};  

import casadi.*

log.root_time = 0;
log.milp_time = [];
log.nlp_time = [];
log.branch_time = [];

[timeout,eq_relax,cont_relax,nav_strat_U,nav_strat_L,U,L,eps,maxiter] = BB_param_check(params);

% set up vars
N = length(f_fun);
real_vars = length(vertcat(x0{:}));
int_vars = length(vertcat(z0{:}));
x = SX.sym('x', [real_vars, 1]);
z = SX.sym('z', [int_vars, 1]);
XZ = [x;z];

%initialize branch and bound
for i = 1:N
    rvars(i) = length(x0{i}); %number of real vars in partition i
    ivars(i) = length(z0{i}); %number of integer vars in partition i
    discrete{i} = [zeros(1,rvars(i)), ones(1,ivars(i))];
    lbz_temp{i} = lbz{i};
    ubz_temp{i} = ubz{i}; 
    
    init{i} = [x0{i};z0{i}];
    AB{i} = [A{i},B{i}];
end
[Sig,lam0,dist_alg_opts] = DBB_subalg_setup(rvars,ivars,length(c),opts); % distributed algorithm params
for i = 1:length(ints)
    d(i) = length(ints{i}); %number of integer decisions in partition i
end
x_sol = x0;
z_sol = z0;
XZsol = [vertcat(x0{:});vertcat(z0{:})];
lbz_temp=lbz;
ubz_temp=ubz;
finish = 0;
iter = 1;
nlp_time = 0;
milp_time = 0;
cur_dec = string(num2str(zeros(1,length(d))));
node_stack = [];

% some helpful functions
f_glob = 0;
for i = 1:N
   f_glob = f_glob + f_fun{i}([x(sum(rvars(1:i-1))+1:sum(rvars(1:i)));...
                               z(sum(ivars(1:i-1))+1:sum(ivars(1:i)))]);
end
f_glob_fun = Function('f',{XZ},{f_glob});

g_glob = concatCasadiFuns(g_fun,rvars,ivars);
g_glob_with_EQcons = [g_glob([x;z]);[horzcat(A{:}), horzcat(B{:})]*[x;z]-c;c-[horzcat(A{:}), horzcat(B{:})]*[x;z]];
g_glob_with_EQcons = Function('f',{XZ},{g_glob_with_EQcons});

%get a first upper bound if (x0,z0) is already feasible 
feas_init = check_feas_MIP(eps,g_glob,horzcat(A{:}),horzcat(B{:}),c,vertcat(lbx{:}),vertcat(lbz{:}),vertcat(ubx{:}),vertcat(ubz{:}),vertcat(x0{:}),vertcat(z0{:}));
if feas_init == 1 
    init_obj = full(f_glob_fun([vertcat(x0{:});vertcat(z0{:})]));
    if init_obj<U
        U = init_obj;
        disp(["U updated by root node: ",U]);
    end
end

% find a point satisfying box and linear equality constraints
[xsol, zsol, eq_feas] = MILP_feas_alg(horzcat(A{:}),horzcat(B{:}),c,vertcat(lbx{:}),vertcat(ubx{:}),vertcat(lbz{:}),vertcat(ubz{:}));
if eq_feas == 0
   disp("Coupling constraints cannot be satisfied")
   finish = 1;
   log.status = "INFEASIBLE";
else
    %get a first upper bound if (xsol, zsol) is feasible w.r.t. inequalities 
    feas2 = check_feas_MIP(eps,g_glob,horzcat(A{:}),horzcat(B{:}),c,vertcat(lbx{:}),vertcat(lbz{:}),vertcat(ubx{:}),vertcat(ubz{:}),xsol,zsol);
    if feas2 == 1
        if full(f_glob_fun([xsol;zsol]))<U
            U = full(f_glob_fun([xsol;zsol]));
            XZsol = [xsol;zsol];
            disp(["U updated by MILP_feas_alg: ",U]);
        end
    end
end

if finish ~= 1
    f1 = 0;
    subtime = zeros(N+1,1);
    %get some info by relaxing equality constraints
    if eq_relax ~=0
        for i = 1:N
            tic
            [relax_sol{i},minlp_status] = solveMINLP(f_fun{i},g_fun{i},[lbx{i};lbz{i}],[ubx{i};ubz{i}],[x_sol{i};z_sol{i}], discrete{i});
            if ~or(strcmp(minlp_status,'OPTIMAL'),strcmp(minlp_status,'SUCCESS'))
               keyboard
               finish = 1;
               log.status = "INFEASIBLE";
               disp('Problem infeasible! Decoupled subproblem not solvable.')
            end
            subtime(i) = toc;
            f1 = f1 + full(relax_sol{i}.f);
        end
    else
        f1 = -inf;
    end
    if finish~=1
        % get some info by relaxing integrality conditions
        if cont_relax~=0
            NLPsolver_opts.print_time = 0;
            [relax_sol2, nlp_log] = solveNLP( f_glob_fun,g_glob_with_EQcons,[vertcat(lbx{:});vertcat(lbz{:})],[vertcat(ubx{:});vertcat(ubz{:})],[vertcat(x_sol{:});vertcat(z_sol{:})],'ipopt',NLPsolver_opts);
            subtime(N+1) = nlp_log.time;
            if ~any(strcmp(good_statuses,nlp_log.status))
                   finish = 1;
                   log.status = "INFEASIBLE";
                   disp('Problem infeasible! Relaxed subproblem not solvable.')
            end
            f2 = full(relax_sol2.f);
        else
            f2 = -inf;
        end
        L = max(f1,f2);
        log.root_time = max(subtime);
        disp(["L updated by root node: ",L]);
        tree = cur_dec; % list of explored nodes
        tree_vals = L; % values of each explored node
        node_stack = [node_stack; getChildren(d,cur_dec)];
    end
else
    log.root_time = 0;
end

log.U = U;
log.L = L;

% main loop
while finish ~= 1
disp(['Starting iteration ',num2str(iter)])
    
tic    
    eq_feas = 0;
    while eq_feas ~= 1
        % get weights associated with each candidate node
        if U==inf
              nodal_weights = [node_stack,getNodeWeights(d,tree,tree_vals,node_stack,nav_strat_U,opts.custom_weights)];
        else
              nodal_weights = [node_stack,getNodeWeights(d,tree,tree_vals,node_stack,nav_strat_L,opts.custom_weights)];        
        end
        % navigate decision tree
        [finish,cur_dec,lbz_temp,ubz_temp,log.status] = navigateBnBTree(nodal_weights,U,L,d,lbz,ubz,ints,int_parts,eps);
        if finish == 1
           break 
        end
        
        %check if local inequality constraints can be satisfied
        [xsol, zsol, eq_feas] = MILP_feas_alg(horzcat(A{:}),horzcat(B{:}),c,vertcat(lbx{:}),vertcat(ubx{:}),vertcat(lbz_temp{:}),vertcat(ubz_temp{:}));
   
        for n = 1:N
           loc_x_cur_decex = 1+sum(rvars(1:n-1)):sum(rvars(1:n-1))+rvars(n);
           x_sol{n} = xsol(loc_x_cur_decex);
           loc_z_cur_decex = 1+sum(ivars(1:n-1)):sum(ivars(1:n-1))+ivars(n);
           z_sol{n} = zsol(loc_z_cur_decex);
        end
          
        if eq_feas == 0
            node_stack(find_tree_ind(node_stack,cur_dec),:) = [];
            [finish,log.status] = terminationCheckBnB(node_stack,U,L,eps);
            disp("Node cannot satisfy equality constraints")
            if finish == 1
               break 
            end
        end

    end
milp_time = toc; 

    %Solve NLP
    if finish~=1    
        if d(max(find(str2num(cur_dec)~=0))) ~= 1 %skip the nlp if the current variable has only one decision
            skip = 0;
            for i = 1:N
               lb_nlp{i} = [lbx{i};lbz_temp{i}];
               ub_nlp{i} = [ubx{i};ubz_temp{i}];
            end
            if strcmp(opts.NLP_solver,'ALADIN')
                [sub_sol, nlp_log] = ALADIN_subAlg( f_fun,g_fun,AB,c,init,lam0,lb_nlp,ub_nlp,Sig,dist_alg_opts, U, L);
            elseif strcmp(opts.NLP_solver,'ADMM')
                [sub_sol, nlp_log] = run_ADMM( f_fun,g_fun,AB,c,init,lam0,lb_nlp,ub_nlp,Sig,dist_alg_opts);
            elseif strcmp(opts.NLP_solver,'ADMM_exp_lasso')
                [sub_sol, nlp_log] = distr_l1_exp(opts.exp_lasso_opts.A0, opts.exp_lasso_opts.b0, lb_nlp, ub_nlp, opts.exp_lasso_opts.mu, N, 1, 1,dist_alg_opts);
            elseif strcmp(opts.subalg,'IPOPT')
                [sub_sol, nlp_log] = solveNLP( f_glob_fun,g_glob_with_EQcons,[vertcat(lbx{:});vertcat(lbz_temp{:})],[vertcat(ubx{:});vertcat(ubz_temp{:})],[vertcat(x_sol{:});vertcat(z_sol{:})] );
            else
                error('Invalid solver selected. Options are: ALADIN, ADMM, ADMM_exp_lasso, IPOPT')
            end
            nlp_time = nlp_log.time;
            if or(strcmp(nlp_log.status,'INFEASIBLE'),strcmp(nlp_log.status,'Infeasible_Problem_Detected'))
               subprob_feas = 0; % the branch will get removed in the fathoming rules
               disp("Infeasible NLP")
            else              
                obj = full(f_glob_fun(sub_sol));
                [x_sol,z_sol,init,lam0,subprob_feas] = organize_DBB_subprob_res(opts.NLP_solver,sub_sol,nlp_log,g_glob,A,B,c,lbx,lbz,ubx,ubz);
            end
        else
            disp("Skipping NLP...")
            skip = 1;
            cur_dec_vec = str2num(cur_dec);
            cur_dec_vec(max(find(cur_dec_vec~=0)))=0;
            prev_dec = num2str(cur_dec_vec);
            tree_ind = find_tree_ind(tree,prev_dec);
            obj = tree_vals(tree_ind(1));
            if and(obj<inf,~isempty(obj))
                subprob_feas = 1;
            else
                subprob_feas = 0;
            end
            node_stack = [node_stack; getChildren(d,cur_dec)];
        end       

    %remove the current node from the candidate list
    node_stack(find_tree_ind(node_stack,cur_dec),:) = [];
    %update node values
    if subprob_feas == 1
       disp(["Node value: ",obj])
       tree = [tree;cur_dec];
       tree_vals = [tree_vals;obj];
    end

tic
    if skip ~= 1
        %check fathoming rules
        if subprob_feas == 0
            disp('Infeasible subproblem.')
        elseif obj>U
            %objective above previous upper bound?
            disp('Previous best upper bound exceeded.')
        elseif sum(abs(mod(vertcat(z_sol{:}),1)))==0
            %integer solution?
            disp('Integer solution detected.') 
            if obj<U
               disp(['New UB: ',num2str(obj)])
               U = obj; 
               XZsol = [vertcat(x_sol{:});vertcat(z_sol{:})];
               node_stack = pruneFromU(U,d,tree,tree_vals,node_stack);
            end
        else
            %add children to node stack
            node_stack = [node_stack; getChildren(d,cur_dec)];
            %update bounds 
            new_U = full(get_U(str2num(cur_dec),x_sol,z_sol,f_fun));
            if new_U<U
               U = new_U;
               disp(['New UB: ',num2str(U)])
               XZsol = [vertcat(x_sol{:});vertcat(z_sol{:})]; 
               node_stack = pruneFromU(U,d,tree,tree_vals,node_stack);
            end
            % try a projection onto the integers
            proj_feas = check_feas_MIP(eps,g_glob,horzcat(A{:}),horzcat(B{:}),c,vertcat(lbx{:}),vertcat(lbz{:}),vertcat(ubx{:}),vertcat(ubz{:}),vertcat(x_sol{:}),round(vertcat(z_sol{:})));
            if proj_feas == 1
                disp("Projection onto integers is feasible.")
                if full(f_glob_fun([vertcat(x_sol{:});round(vertcat(z_sol{:}))]))<U
                   U = full(f_glob_fun([vertcat(x_sol{:});round(vertcat(z_sol{:}))]));
                   disp(['New UB: ',num2str(U)])
                   XZsol = [vertcat(x_sol{:});round(vertcat(z_sol{:}))];  
                   node_stack = pruneFromU(U,d,tree,tree_vals,node_stack);
                end
            end
            L = min(get_Ls(tree,tree_vals,node_stack,d));
            disp(["Current lower bound: ",L])
        end
    end
    end
branch_time = toc;

    log.milp_time = [log.milp_time, milp_time];
    log.nlp_time = [log.nlp_time, nlp_time];
    log.branch_time = [log.branch_time, branch_time];
    log.U = [log.U, U];
    log.L = [log.L, L];
    
    iter = iter + 1; 
    
    %termination check
    [finish,log.status] = terminationCheckBnB(node_stack,U,L,10^-3);
    
    if iter>maxiter
       finish = 1;
       disp('Iteration limit reached.')
       log.status = "TIMEOUT";
    end
    if sum(log.milp_time)+sum(log.nlp_time)+sum(log.branch_time)>timeout
       finish = 1;
       disp('Time limit reached.')
       log.status = "TIMEOUT";
    end
    
end

end
