close all;
clear all;
clc;

import casadi.*

rand('seed', 0);
randn('seed', 0);

n = 4; %number of features
m = 10; %number of training examples
N = 10; %nuber of subsystems

if mod(m,N)~=0
   error('m must be divisible by N!') 
end

% algorithm parameters
eps = 10^-1; % Termination threshold
max_iter = 500; % Iteration limit
max_time = 3000; % time limit;

w = sprandn(n, 1, 100/n);       % N(0,1), 10% sparse
v = randn(1);                  % random intercept

X0 = sprandn(m*N, n, 10/n);           % data / observations
btrue = sign(X0*w + v);

% noise is function of problem size use 0.1 for large problem
b0 = sign(X0*w + v + sqrt(0.1)*randn(m*N, 1)); % labels with noise

% packs all observations in to an m*N x n matrix
A0 = spdiags(b0, 0, m*N, m*N) * X0;

x_true = [v; w];

% objective and constraints
ratio = sum(b0 == 1)/(m*N);
mu = 0.1*1/(m*N) * norm((1-ratio)*sum(A0(b0==1,:),1) + ratio*sum(A0(b0==-1,:),1), 'inf');
C = [-b0 -A0]';
x = SX.sym('x',[n+1,1]);
A_cons = eye((n+1)*N) - diag(ones((n+1)*(N-1),1),n+1); A_cons((n+1)*(N-1)+1:(n+1)*N,1:n+1)=-eye(n+1);
obj_cent = 0;
x_glob = SX.sym('x',[(n+1)*N,1]);
for part = 1:N
    c_loc = C(:,1+(part-1)*m:part*m)';
    obj{part} = sum(exp(c_loc*x))+m*mu*norm(x(2:end),1);
    obj_funs{part} = Function('f',{x},{obj{part}});
    lbx{part} = -3*ones(n,1); ubx{part}=3*ones(n,1);
    lbz{part} = -3*ones(1,1); ubz{part}=3*ones(1,1);
    lb{part} = [lbx{part};lbz{part}];
    ub{part} = [ubx{part};ubz{part}];
    loc_cons = [x-ub{part};lb{part}-x];
    cons_funs{part} = Function('g',{x},{loc_cons});
    A_con{part} = A_cons(:,(part-1)*(n+1)+1:part*(n+1)-1); %consensus in real vars
    B_con{part} = A_cons(:,part*(n+1)); %consensus in discrete vars
    AB{part} = [A_con{part},B_con{part}];
    discrete{part} = zeros(1,n+1); discrete{part}(n+1)=1;
    x0{part} = zeros(n,1);
    z0{part} = 0;
    u0{part} = [x0{part};z0{part}];
    
    x_loc = x_glob((part-1)*(n+1)+1:part*(n+1));
    obj_cent = obj_cent + obj_funs{part}(x_loc);
end
b_con = zeros((n+1)*N,1);

lbz_glob = vertcat(lbz{:});
ubz_glob = vertcat(ubz{:});

int_parts = [];
p_count = 1;
v_count = 1;
for i = 1:length(lbz_glob)
   ints{i} = lbz_glob(i):ubz_glob(i); 
   int_parts = [int_parts;[i,p_count,v_count]]; %ints{i} belongs to partition p_count and is the v_count^th integer of that partition
   if i==length(vertcat(lbz{1:p_count}))
      p_count = p_count + 1;
      v_count = 1;
   else
       v_count = v_count + 1;
   end
end

%% Solve problem
obj_fun_cent = Function('f',{x_glob},{obj_cent});
for i = 1:N
   AB{i} = [A_con{i} B_con{i}]; 
end
cons_cent = [horzcat(AB{:})*x_glob-b_con; -horzcat(AB{:})*x_glob+b_con];
cons_fun_cent = Function('g',{x_glob},{cons_cent});
opts.discrete = horzcat(discrete{:});
opts.print_time = 0;% opts.verbose_init = 0;
tic
nlp =   struct('x',x_glob,'f',obj_cent,'g',cons_cent);
cas =   nlpsol('solver','bonmin',nlp,opts);
sol =   cas('lbx', vertcat(lb{:}),...
            'ubx', vertcat(ub{:}),...
            'lbg', -inf(length(cons_cent),1),...
            'ubg', zeros(length(cons_cent),1),...
            'x0', vertcat(u0{:}));
Cent_time = toc;

%% Dist sol
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

opts.exp_lasso_opts.A0 = A0;
opts.exp_lasso_opts.b0 = b0;
opts.exp_lasso_opts.mu = mu;
opts.exp_lasso_opts.QUIET    = 1;
opts.exp_lasso_opts.MAX_ITER = 35;
opts.exp_lasso_opts.ABSTOL   = 1e-4;
opts.exp_lasso_opts.RELTOL   = 1e-2;
opts.NLP_solver = 'ADMM_exp_lasso';
[XZsol2,log2] = Dist_BnB(obj_funs,cons_funs,lbx,lbz,ubx,ubz,A_con,B_con,b_con,x0,z0,ints,int_parts,params,opts);
DBB_time = log2.root_time + sum(log2.milp_time) + sum(log2.nlp_time) + sum(log2.branch_time);

%% display results
disp('Bonmin time   DBB time')
disp([Cent_time  DBB_time])
disp('  Bonmin obj.    DBB obj.')
disp([full(sol.f)  log2.U(end)])
