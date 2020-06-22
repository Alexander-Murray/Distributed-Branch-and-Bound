% Based on the ADMM algorithm from 
% "Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers"
% original MATLAB code available at:
% https://web.stanford.edu/~boyd/papers/admm/logreg-l1/logreg.html
function [z, history] = distr_l1_exp(A, b, lb, ub, mu, N, rho, alpha,opts)

t_start = tic;

% Preprocessing
[m, n] = size(A);
m = m / N;  % should be divisible
% ADMM solver
x = zeros(n+1,N);
z = zeros(n+1,N);
u = zeros(n+1,N);


if ~opts.QUIET
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', '# bfgs', ...
      'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
end
p = size(z,1);
C = [-b -A]';

history.status = 'in progress';
history.par_time=[];
history.seq_time=[];
for k = 1:opts.MAX_ITER
    par_time = zeros(N,1);
    % serial x-update
    for i = 1:N
        if ~opts.QUIET
            disp(['iter: ', 'subprob: ']); disp([k,i]);
        end
        
        K = C(:,1+(i-1)*m:i*m)';
        [x(:,i), par_time(i)] = x_update(lb{i},ub{i},K, u(:,i), z(:,i), rho, x(:,i));
    end
    history.par_time = [history.par_time, par_time];
    
    tic
    % z-update with relaxation
    zold = z;
    x_hat = alpha*x + (1-alpha)*zold;
    ztilde = mean(x_hat + u,2);
    ztilde(2:end) = shrinkage( ztilde(2:end), (m*N)*mu/(rho*N) );

    z = ztilde*ones(1,N);

    % u-update
    u = u + (x_hat - z);
    
    seq_time = toc;
    history.seq_time = [history.seq_time, seq_time];

    % diagnostics, reporting, termination checks
    history.objval(k)  = objective(A, b, mu, x, z(:,1));

    history.r_norm(k)  = norm(x - z, 'fro');
    history.s_norm(k)  = norm(rho*(z - zold),'fro');
    history.eps_pri(k) = sqrt(p*N)*opts.ABSTOL + opts.RELTOL*max(norm(x,'fro'), norm(z,'fro'));
    history.eps_dual(k)= sqrt(p*N)*opts.ABSTOL + opts.RELTOL*norm(rho*u,'fro');

    if or(history.r_norm(k) < history.eps_pri(k) && history.s_norm(k) < history.eps_dual(k), k==opts.MAX_ITER)
        z = reshape(z,(n+1)*N,1);
        history.time = sum(history.seq_time) + sum(max(history.par_time));
        history.status = 'SUCCESS';
        history.obj = history.objval(end);
        if ~opts.QUIET
            K = length(history.objval);

            h = figure;
            plot(1:K, history.objval, 'k', 'MarkerSize', 10, 'LineWidth', 2);
            ylabel('f(x^k) + g(z^k)'); xlabel('iter (k)');

            g = figure;
            subplot(2,1,1);
            semilogy(1:K, max(1e-8, history.r_norm), 'k', ...
                1:K, history.eps_pri, 'k--',  'LineWidth', 2);
            ylabel('||r||_2');

            subplot(2,1,2);
            semilogy(1:K, max(1e-8, history.s_norm), 'k', ...
                1:K, history.eps_dual, 'k--', 'LineWidth', 2);
            ylabel('||s||_2'); xlabel('iter (k)');

            disp(['par_time: ','seq_time: '])
            disp([sum(max(history.par_time)), sum(history.seq_time)])
        end
        
        break;
    end
end

if ~opts.QUIET
    toc(t_start);
end
z = z(:,1);
end

function obj = objective(A, b, mu, x, z)
    m = size(A,1);
    obj = sum(exp(-A*z(2:end) -b*z(1))) + m*mu*norm(z(2:end),1);
end

function [x t] = x_update(lb,ub,C, u, z, rho, x0)
    % solve the x update
    auxdata{1} = C;
    auxdata{2} = z;
    auxdata{3} = u;
    auxdata{4} = rho;   

    loc_obj = @(x)sum(exp(C*x)) + (1/2)*(x - z + u)'*rho*(x - z + u);

tic
opts = optimset('Display','off','Algorithm','interior-point', 'MaxIter', 10000, 'MaxFunEvals', 10000);
x = fmincon(loc_obj,x0,[],[],[],[],lb,ub,[],opts);
t=toc;
end


function z = shrinkage(a, kappa)
    z = max(0, a-kappa) - max(0, -a-kappa);
end