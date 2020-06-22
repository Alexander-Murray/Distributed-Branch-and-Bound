close all;
clear all;
clc;

import casadi.*

%% Problem parameters
N = 10; % number of partitions
NC = 5; % number of consensus constraints
nx = 10*ones(N,1); % number of real vars in each partition
nz = 1*ones(N,1); % number of discrete vars in each partition
nxz = nx+nz;

%% Problem setup
for i = 1:N
    x = SX.sym(strcat('x',num2str(i),'_'),nx(i),1); % local real vars
    z = SX.sym(strcat('z',num2str(i),'_'),nz(i),1); % local discrete vars
    
    % objective function
    obj{i} = rand(1,nx(i))*x + rand(1,nz(i))*z;
    obj_funs{i} = Function('f',{[x;z]},{obj{i}});
    
    % local constraints
    cons{i} = x.'*x + i*z - 5*rand;
    cons_funs{i} = Function('f',{[x;z]},{cons{i}});
    
    % consensus constraint matrices
    A{i} = rand(NC,nx(i));
    B{i} = zeros(NC,nz(i));
    AB{i} = [A{i}, B{i}];
    
    % initial guess
    x0{i} = zeros(nx(i),1);
    z0{i} = zeros(nz(i),1);
    
    % box constraints
    lbx{i} = -2*ones(nx(i),1);
    ubx{i} = 2*ones(nx(i),1);
    lbz{i} = -2;
    ubz{i} = 2;
    ints{i} = lbz{i}:ubz{i};
    
    % vector denoting which variables are integer valued
    discrete{i}=[zeros(1,nx(i)),ones(1,nz(i))];
end
c = zeros(NC,1); % RHS of consensus constraints
lam0 = zeros(NC,1); % initial guess of Lagrange multiplier of consensus constraints

%% Centralized solution

% Centralized problem formulation
x_glob =  SX.sym('x',[sum(nxz),1]);

obj_cent_temp = 0;
cons_cent_temp= horzcat(AB{:})*x_glob-c;
for i = 1:N
    x_part = x_glob(sum(nx(1:i-1))+1:sum(nx(1:i)));
    z_part = x_glob(sum(nx)+sum(nz(1:i-1))+1:sum(nx)+sum(nz(1:i)));

    obj_cent_temp = obj_cent_temp + obj_funs{i}([x_part;z_part]);
    
    cons_cent_temp = [cons_cent_temp;cons_funs{i}([x_part;z_part])];
end
obj_cent = Function('f',{x_glob},{obj_cent_temp});

cons_cent = Function('g',{x_glob},{cons_cent_temp});

discrete_cent = [zeros(1,sum(nx)),ones(1,sum(nz))];

lbg = [zeros(NC,1);-Inf(length(cons_cent_temp)-NC,1)];
ubg = [zeros(NC,1);zeros(length(cons_cent_temp)-NC,1)];

% solve with bonmin
tic
cent_opts.verbose_init = 0;
cent_opts.bonmin.print_level = 0;
cent_opts.bonmin.bb_log_level = 0;
cent_opts.bonmin.fp_log_level = 0;
cent_opts.bonmin.lp_log_level = 0;
cent_opts.bonmin.milp_log_level = 0;
cent_opts.bonmin.nlp_log_level = 0;
cent_opts.bonmin.oa_log_level = 0;
cent_opts.print_time = 0;
cent_opts.discrete = discrete_cent;

nlp =   struct('x',x_glob,'f',obj_cent(x_glob),'g',cons_cent(x_glob));
cas =   nlpsol('solver','bonmin',nlp,cent_opts);
sol =   cas('lbx', [vertcat(lbx{:});vertcat(lbz{:})],...
            'ubx', [vertcat(ubx{:});vertcat(ubz{:})],...
            'lbg', lbg,...
            'ubg', ubg,...
            'x0', [vertcat(x0{:});vertcat(z0{:})]);
Cent_time = toc;

cas.stats.return_status
Cent_sol = full(sol.x);

%% Branch and Bound solution
params.timeout = 300;
params.eq_relax = 0;
params.cont_relax = 1;
params.nav_strat_U = 'best-first';
params.nav_strat_L = 'best-first';
opts.custom_weights = [];
opts.NLP_solver = 'ipopt';
opts.NLPsolver_opts.print_time = 0;
opts.NLPsolver_opts.print_time = 0;
opts.NLPsolver_opts.ipopt.print_level = 0;
opts.NLPsolver_opts.verbose = 0;
opts.NLPsolver_opts.verbose_init = 0;

SBB_obj{1} = obj_cent;
SBB_con{1} = cons_cent;
SBB_lbx{1} = vertcat(lbx{:});
SBB_ubx{1} = vertcat(ubx{:});
SBB_lbz{1} = vertcat(lbz{:});
SBB_ubz{1} = vertcat(ubz{:});
SBB_A{1} = horzcat(A{:});
SBB_B{1} = horzcat(B{:});
SBB_x0{1} = vertcat(x0{:});
SBB_z0{1} = vertcat(z0{:});
[XZsol,log] = Standard_BnB(SBB_obj,SBB_con,SBB_lbx,SBB_lbz,SBB_ubx,SBB_ubz,SBB_A,SBB_B,c,SBB_x0,SBB_z0,ints,params,opts);

SBB_time = log.root_time + sum(log.milp_time) + sum(log.nlp_time) + sum(log.branch_time);

%% display solution
disp('Bonmin time')
disp(Cent_time)
disp('SBB time')
disp(SBB_time)
disp('SBB obj. / Bonmin obj.')
disp(full(obj_cent(XZsol)/sol.f))
disp('|SBB x - Bonmin x|')
disp(sum(abs(XZsol-full(sol.x))))
