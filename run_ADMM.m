function [ xopt, logg ] = run_ADMM( ffi,hhi,AA,b,yy0,llam0,llbx,uubx,Sig,dist_alg_opts )
tic
good_statuses = {'Search_Direction_Becomes_Too_Small','Solved_To_Acceptable_Level','Solve_Succeeded'};
   
rho=dist_alg_opts.rho0;
rho_update=dist_alg_opts.rhoUpdate;
eps=dist_alg_opts.maxiter;
maxiter=dist_alg_opts.maxiter;

NsubSys = length(ffi);
Ncons   = size(AA{1},1);

%% built local subproblems and CasADi functions
import casadi.*
rhoCas      = SX.sym('rho',1,1);
kapp{NsubSys,1}=[];
nx = zeros(NsubSys,1);
nnlp{NsubSys,1}=[];
lamCas   = SX.sym('lam',Ncons,1);
for i=1:NsubSys
    nx(i)    = length(yy0{i});
    yyCas    = SX.sym('z',nx(i),1);
    xxCas    = SX.sym('x',nx(i),1);
    nnhi{i} = length(hhi{i}(yyCas));
    
    % parameter vector of CasADi
    pCas        = [ rhoCas;
                    lamCas;
                    xxCas];
                
    % objective function for local NLP's
    ffiLocCas = ffi{i}(yyCas) + lamCas'*AA{i}*yyCas ...
                + rhoCas/2*(AA{i}*(yyCas - xxCas))'*Sig{i}*(AA{i}*(yyCas - xxCas));
  
    
    % local inequality constraints
    hhiCas  = hhi{i}(yyCas);
                
    % set up local solvers
    nlp_opts.print_time = 0;
    nlp_opts.ipopt.print_level = 0;
    nlp     = struct('x',yyCas,'f',ffiLocCas,'g',hhiCas,'p',pCas);
    nnlp{i} = nlpsol('solver','ipopt',nlp,nlp_opts);
end
A   = horzcat(AA{:});

xCas = SX.sym('x', [sum(nx), 1]);
cent_obj = 0;
var_count = 0;
for fun = 1:NsubSys
cent_obj = cent_obj + ffi{fun}(xCas(var_count+1:var_count+nx(fun)));
var_count = var_count+nx(fun);
end 
cent_obj_fun = Function('f',{xCas},{cent_obj});

logg        = struct();
logg.Y      = vertcat(yy0{:});
logg.lambda = vertcat(llam0{:});
logg.Kappa  = [];
logg.obj  = full(cent_obj_fun(vertcat(yy0{:})));
logg.status = 'timeout';
logg.EqViol = sum(abs(A*vertcat(yy0{:}) - b)); 
logg.time = 0;
logg.par_time = [];
logg.seq_time = [];

figure(1);
set(gcf,'position',[10 30 400 300])
    hold on;
    grid on;
    plot(1,logg.obj(1),'-mo','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63]);
    title('Objective Function Value')
    xlabel('Iterations')   

%track consensus constraint violation
figure(2);
set(gcf,'position',[420 30 400 300])
    hold on;
    grid on;
    plot(1,logg.EqViol(1),'-mo','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63]);
    title('Equality Constraint Violation')
    xlabel('Iterations')
%     axis([1 maxiter 0 1])

figure(3)
ylim([0,1])

%% ADMM iterations
% initialization
i       = 1;
yy      = yy0;
xx      = yy0;
llam    = llam0;
infeas = 0;

while i<=maxiter% && norm(delx,inf)>eps   
    subtime = zeros(1,NsubSys);
    for j=1:NsubSys
        % set up parameter vector for local NLP's
        pNum = [ rho;
                 llam{j};
                 xx{j}];
        tic                           
        % solve local NLP's
        sol = nnlp{j}('x0' , yy{j},...
                      'p',   pNum,...
                      'lbx', llbx{j},...
                      'ubx', uubx{j},...
                      'lbg', -inf*ones(nnhi{j},1), ...
                      'ubg', zeros(nnhi{j},1));     
        if ~any(strcmp(good_statuses,nnlp{j}.stats.return_status))          
            warning( nnlp{j}.stats.return_status )
            if i > 1
                infeas = 1;
                break  
            end
        end
        subtime(i) = toc;
        
        yy{j}           = full(sol.x);
        kapp{j}         = full(sol.lam_g);

        % multiplier update
        llam{j} = llam{j} + rho*AA{j}*(yy{j}-xx{j});
    end
    if infeas == 1
       logg.status = 'INFEASIBLE';
       y = vertcat(yy0{:});
       break
    end
    logg.par_time = [logg.par_time, max(subtime)];
    % global x vector
    y = vertcat(yy{:});
    
    % Termination check
    InOpt = full(cent_obj_fun(y));
    ConViol = sum(abs(A*y-b).^2);%full(max(abs(A*y - b)));
    if ConViol<eps && abs(InOpt-logg.obj(i))<eps
       logg.status = 'SUCCESS';
       break 
    end

    if and(max(abs(y-logg.Y(:,i)))<10^-5,abs(ConViol-logg.EqViol(i))<10^-5)
        disp('ADMM is stuck. Terminating...')
        break
    end
  
    % Solve ctr. QP
    hQP=[];
    for j=1:NsubSys
       hQP  = sparse([hQP, -rho*yy{j}'*AA{j}'*AA{j}-llam{j}'*AA{j}]);
    end   
    
    %% build Ha nd A for ctr. QP
    HQP = [];
    for j=1:NsubSys
        HQP = sparse(blkdiag(HQP, rho*AA{j}'*AA{j} + eye(nx(j),nx(j))*1e-9)) ; % regularization, because some states aren't involved in the solution
    end   
    
    % solve QP
    tic
    [x, ~] = solveQP(HQP,hQP',A,b,'backslash');
    QP_time = toc;
    logg.seq_time = [logg.seq_time,QP_time];

    % divide into subvectors
    ctr = 1;
    for j=1:NsubSys
        ni          = length(xx{j});
        xx{j}       = x(ctr:(ctr+ni-1)); 
        ctr = ctr + ni;
    end
    
    i = i+1;
    
    % plotting
    dual_res = sum(abs(rho*(vertcat(llam{:})-logg.lambda(:,i-1))).^2);
    primal_res = ConViol;

    figure(1);
    hold on;
    grid on;
    plot(i+1,InOpt,'-mo','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63]);

    figure(2);
    hold on;
    grid on;
    plot(i+1,ConViol,'-mo','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63]);
    
    % logging
    logg.Y          = [logg.Y y];
    logg.lambda     = [logg.lambda vertcat(llam{:})];
    logg.Kappa      = [logg.Kappa vertcat(kapp{:})];
    logg.obj        = [logg.obj InOpt];
    logg.EqViol     = [logg.EqViol ConViol];
    
    figure(3)
    hold on
    grid on;
    plot(i,dual_res,'-mo','MarkerEdgeColor','k','MarkerFaceColor',[.49 0 .63])
    
    if ConViol-logg.EqViol(i-1)>-10^-3 % increase rho if decrease is not fast enough
       rho = rho*rho_update; 
    end
    if i>2 && logg.obj(i-1)-logg.obj(i-2)>0
    if InOpt-logg.obj(i-1)>4*(logg.obj(i-1)-logg.obj(i-2))
       rho = rho/rho_update; 
    end
    end
    rho = rho*1.02;

    disp(["rho = ",rho])

end
close(1); close(2); close(3);
hold on
figure(4)
plot(1:i,logg.obj);
xlim([0 maxiter])
hold off
hold on
figure(5)
plot(1:i,logg.EqViol);
xlim([0 maxiter])
hold off

% Comment: The solution y of the local NLP's is used, because not all
% decision variables are considered in the global QP. Because of the
% regualrization step, they are set to zero what makes the solution
% unusable.
xopt = y;
logg.time = sum(logg.par_time)+sum(logg.seq_time);
logg.total_time = toc; %includes time needed to plot figures, etc.
end

