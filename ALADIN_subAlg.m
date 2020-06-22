% Solves a non-convex optimization problem in consensus form via the ALADIN
function [ xopt, logg , tot_time ] = ALADIN_subAlg( obj_funs,cons_funs,A_mats,b,x0,...
                                    lam0,lbx,ubx,Sig,opts, U, L)
%get solver for QP step
QPsolver = opts.solveQP;
% get dimensions
NsubSys    = length(obj_funs); %number of partitions
Ncons      = size(A_mats{1},1); %number of equality constraints
% build global A matrix
A          = horzcat(A_mats{:});

%% build local subproblems and CasADi functions
import casadi.*
rhoCas      = SX.sym('rho',1,1);
lamCas      = SX.sym('lam',length(lam0),1);

tic
nx = zeros(NsubSys,1);
for i=1:NsubSys  
%     discrete{i} = opts.discrete{i}; %1 denotes integer, 0 denotes real-valued variable
    discrete{i} = zeros(1,length(x0{i}));
    nx(i)    = length(x0{i}); %number of variables in partition i
    xCas    = SX.sym('y',nx(i),1); %local decision variables
    diff_x   = find(ones(1,nx(i))-discrete{i}); %identify real variables for differentiation
    yCas    = SX.sym('x',nx(i),1); %parameterize QP solution
    x_opt{i} = x0{i}; %initial point given to local solvers
    
    % parameter vector of CasADi
    pCas        = [ rhoCas;
                    lamCas;
                    yCas];
                
    % objective function for local NLP's
    obj_fun_Loc_Cas = obj_funs{i}(xCas) + lamCas'*A_mats{i}*xCas ...
                + rhoCas/2*(xCas - yCas)'*Sig{i}*(xCas - yCas);

    
    %Gradient of local objective
    gCas    = gradient(obj_funs{i}(xCas),xCas(diff_x));
    g_fun{i}    = Function(['g' num2str(i)],{xCas},{gCas});
    
    
    % local inequality constraints
    cons_Cas  = cons_funs{i}(xCas);

    % Jacobian of local constraints
    Jac_Cas   = jacobian(cons_Cas,xCas(diff_x));   
    Jac_Fun{i}   = Function(['Jac' num2str(i)],{xCas},{Jac_Cas});
    Nloc_cons{i}=size(Jac_Cas,1); %number of local inequality constraints

    %Hessian  of local objective 
    kappa_Cas = SX.sym('kapp',Nloc_cons{i},1);
    rho_Cas = SX.sym('rho',1,1);
    %The third term of the Hessian is to ensure positive-definiteness.
    %Otherwise, the QP will be non-convex
    Hess_Loc_Cass    = hessian(obj_funs{i}(xCas)+kappa_Cas'*cons_Cas,xCas(diff_x)) + rho_Cas*Sig{i}(diff_x,diff_x); 
    Hess_Loc_Fun{i}    = Function(['H' num2str(i)],{xCas,kappa_Cas,rho_Cas},{Hess_Loc_Cass}); 

    % set up local solvers
    solver_struct = struct('x',xCas,'f',obj_fun_Loc_Cas,'g',cons_Cas,'p',pCas);
    nlp_opts.print_time = 0;
    nlp_opts.ipopt.print_level = 0;
    if sum(discrete{i}) == 0
        nnlp{i} = nlpsol('solver','ipopt',solver_struct,nlp_opts);
    else
        MINLPsolver_opts.discrete = discrete{i};
        minnlp{i} = nlpsol('solver','bonmin',solver_struct,opts.MINLPsolver_opts);
        nnlp{i} = nlpsol('solver','ipopt',solver_struct,nlp_opts); %needed for obtaining langrange multipliers from MIP subprobs
    end
           
end
toc

%% Initialization
i = 1;
yy      = x0;
lam     = lam0;
delx = inf;
rho             = opts.rho0;
mu              = opts.mu0;
tot_time = 0;
stop = 0;

%setup a log for some key vectors
logg         = struct();
logg.X       = vertcat(x0{:});
logg.delY    = [];   
logg.Kappa   = [];
logg.lambda  = [];
logg.obj  = [];
logg.ConViol = [];
logg.status = 'timeout';

%track objective function compared to centralized solution

xx = SX.sym('x', [sum(nx), 1]);
cent_obj = 0;
var_count = 0;
for fun = 1:NsubSys
cent_obj = cent_obj + obj_funs{fun}(xx(var_count+1:var_count+nx(fun)));
var_count = var_count+nx(fun);
end  
cent_obj_fun = Function('f',{xx},{cent_obj});
if ~isempty(opts.cent_sol) 
    figure(1);
    set(gcf,'position',[10 30 400 300])
        hold on;
        grid on;
        plot(opts.maxiter,full(cent_obj_fun(opts.cent_sol)),'-mo','MarkerEdgeColor','k','MarkerFaceColor',[.49 0 .63]);
        title('Objective Function Value')
        xlabel('Iterations')   
    logg.ObjVal  = full(cent_obj_fun(logg.X(:,1)));
end
%track consensus constraint violation
ConViol = sum(abs(A*vertcat(yy{:}) - b)); 
figure(2);
set(gcf,'position',[420 30 400 300])
    hold on;
    grid on;
    plot(1,ConViol,'-mo','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63]);
    title('Equality Constraint Violation')
    xlabel('Iterations')
    axis([1 opts.maxiter 0 1])

%%

while i <= opts.maxiter  && stop == 0
    
    sub_time = zeros(1,NsubSys);

    for j=1:NsubSys
        tic
        % set up parameter vector for local NLP's
        pNum = [ rho;
                 lam;
                 yy{j}]; 
             
        % solve local MINLP's   
        if sum(discrete{j})>0
            sol = minnlp{j}('x0' , x_opt{j},...
                          'p',   pNum,...
                          'lbx', lbx{j},...
                          'ubx', ubx{j},...
                          'lbg', -inf*ones(Nloc_cons{j},1), ...
                          'ubg', zeros(Nloc_cons{j},1));

            %resolve with fixed integer vars to obtain kappas
            relbx{j}=lbx{j};
            reubx{j}=ubx{j};
            relbx{j}(find(discrete{j}))=full(sol.x(find(discrete{j})));
            reubx{j}(find(discrete{j}))=full(sol.x(find(discrete{j})));

            resol= nnlp{j}('x0' , full(sol.x).',...
                          'p',   pNum,...
                          'lbx', relbx{j},...
                          'ubx', reubx{j},...
                          'lbg', -inf*ones(Nloc_cons{j},1), ...
                          'ubg', zeros(Nloc_cons{j},1));

            if ~or(strcmp(minnlp{j}.stats.return_status,'SUCCESS'),strcmp(minnlp{j}.stats.return_status,'OPTIMAL'))
                %this may make the solution jump around
                sol.x = resol.x;
    %             keyboard
               warning(  '%s, Iteration: %i, Partition: %i',minnlp{j}.stats.return_status,i,j)
            end

            %Store local solutions
            x_opt{j}           = full(sol.x);
            kappa_opt{j}       = full(resol.lam_g);

       else
                sol = nnlp{j}('x0' , x_opt{j},...
                          'p',   pNum,...
                          'lbx', lbx{j},...
                          'ubx', ubx{j},...
                          'lbg', -inf*ones(Nloc_cons{j},1), ...
                          'ubg', zeros(Nloc_cons{j},1));

            if ~or(strcmp(nnlp{j}.stats.return_status,'Solve_Succeeded'),...
                   strcmp(nnlp{j}.stats.return_status,'Solved_To_Acceptable_Level'))
               warning( nnlp{j}.stats.return_status )
            end

            %Store local solutions
            x_opt{j}           = full(sol.x);
            kappa_opt{j}       = full(sol.lam_g);
         
        end
      
        % active set detection
        subprob_cons = full(cons_funs{j}(x_opt{j}));
        inact = subprob_cons<-1e-5;
        kappa_opt{j}(inact)= 0;

        % Jacobians of inequality constraints 
        %(lower left component of AQP)
        Ci         = full(Jac_Fun{j}(x_opt{j}));
        for d = 1:length(discrete{j})
            if discrete{j}(d)==1
                Ci(:,d) = 0;
            end
        end
        C{j}       = Ci(~inact,:); % eliminate inactive entries

        % evaluate gradients and hessians for the QP
        Hess_Loc_FunEval{j} = Hess_Loc_Fun{j}(x_opt{j},kappa_opt{j},rho);
        g_fun_eval{j}       = g_fun{j}(x_opt{j});
        
        % eigenvalue decomposition of the hessian
        [V,D]       = eig(full(Hess_Loc_FunEval{j}));
        e           = diag(D);

%         % flip the eigenvalues 
%         e           = abs(e);                       

        % modify zero and negative eigenvalues (regularization)
        reg         = 1e-4;
        e(e<reg)    = reg; % 1e-4

        Hess_Loc_FunEval{j}  = V*diag(e)*transpose(V);

         sub_time(i) = toc;
    end       
    tot_time = tot_time + max(sub_time);
    
    %Termination Check
    %%%% Add check on objective value %%%%
    if max(abs(horzcat(A_mats{:})*vertcat(x_opt{:})-b))<opts.eps
       xopt = vertcat(x_opt{:});
       disp("Terminating due to EQ constraint satisfaction...")
       logg.status = 'SUCCESS';
       break 
    end
    
    % build the QP
    rhsQP  = b-A*vertcat(x_opt{:});    
    
    C_QP   = rref(blkdiag(C{:})); 

    %remove zero rows from C_QP
    CQP_rows=size(C_QP,1);
    r=1;
    while r <= CQP_rows
        if sum(abs(C_QP(r,:)))==0
            C_QP(r,:)=[];
            CQP_rows = CQP_rows-1;
        else
            r=r+1;
        end
    end
    
    %linsolve and pinv don't seem to allow for sparse matrices
    if or(strcmp(QPsolver,'linsolve'),strcmp(QPsolver,'pinv'))
        HQP     = full(blkdiag(Hess_Loc_FunEval{:},mu*eye(Ncons)));     
    else
        HQP     = sparse(blkdiag(Hess_Loc_FunEval{:},mu*eye(Ncons)));    
    end
        gQP     = full(vertcat(g_fun_eval{:},lam));
    
        
    AQP     = [A  -eye(Ncons); C_QP  zeros(CQP_rows,Ncons)];
    bQP     = [rhsQP;zeros(CQP_rows,1)];
    
    %check positive-definiteness of HQP. TURN OFF FOR SPEED
%     [~,p] = chol(HQP);
%     if and(p>0,or(strcmp(QPsolver,'linsolve'),strcmp(QPsolver,'backslash')))
%        warning('linsolve may not return optimal value for non-positive definite hessians')
%        keyboard
%     end
    
    %remove ints from AQP
    ints = find(horzcat(discrete{:}));
    AQP(:,ints)=[];
    
    % solve global QP
    tic
    [QP_sol, lamges] = solveQP(HQP,gQP,AQP,bQP,QPsolver);
    tot_time = tot_time + toc;
    delx = QP_sol(1:(end-Ncons)); 
    
 
alpha_1 = 1; alpha_2 = 1; alpha_3 = 1;

%%
    % lambda only for cons. constr.
    lam = lam + alpha_3*(lamges(1:Ncons) - lam);

    % primal variable update
     ctr = 1;
    for j=1:NsubSys
        delta_x{j} = zeros(nx(j),1);
        for var = 1:nx(j)
                delta_x{j}(var) = delx(ctr);
                ctr = ctr + 1;
        end
        
        yy{j} = yy{j} + alpha_1*(x_opt{j}-yy{j}) + alpha_2*delta_x{j}; 
        
    end 
  
    
    % logging
    logg.X          = [logg.X vertcat(x_opt{:})];
    logg.delY       = [logg.delY vertcat(delta_x{:})];
    logg.Kappa      = [logg.Kappa vertcat(kappa_opt{:})];
    logg.lambda     = [logg.lambda lam];
    logg.obj     = [logg.obj full(cent_obj_fun(vertcat(x_opt{:})))];
    logg.ConViol    = [logg.ConViol A*vertcat(x_opt{:}) - b];
    
    % early termionation check
    if i>opts.min_iter
        if max(abs(horzcat(A_mats{:})*vertcat(x_opt{:})-b))<logg.ConViol(end) %convergence has started
        if full(logg.ObjVal(i-1)) > full(obj) %converge from above
            if full(obj) < L
                logg.status = "Terminated by L"; 
                stop = 1;
            end
        else %converge from below
            if full(obj) > U
                logg.status = "Terminated by U"; 
                stop = 1;
            end
        end
        end
    end
    
    % plotting
    InOpt = cent_obj_fun(logg.X(:,i+1)+logg.delY(:,i));
    ConViol = max(abs(A*vertcat(x_opt{:}) - b));

    x = i;

    figure(1);
    y = full(InOpt);
    hold on;
    grid on;
    plot(x,y,'-mo','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63]);

    figure(2);
    z = full(ConViol);
    hold on;
    grid on;
    plot(x+1,z,'-mo','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63]);

    
    %next iteration
    i = i+1;
    
    %updating rho and mu slightly at each iteration can help improve
    %convergence, although this is a heuristic tool
    if rho < opts.rhoMax
       rho = rho*opts.rhoUpdate;
       sprintf('rho: %d',rho)
    end
    if mu < opts.muMax
       mu = mu*opts.muUpdate;
       sprintf('mu: %d',mu)
    end
    
     if i > opts.maxiter
        disp("Terminating due to iteration limit...") 
     end
     if stop==1
        disp("Terminating due to B&B bounds...") 
     end
    

end

xopt = vertcat(x_opt{:});
logg.time = tot_time;

end
