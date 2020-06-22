
function [ xopt, lam ] = solveQP( H, g, A, b, solver )
if strcmp(solver, 'ipopt')
    import casadi.*
    neq = size(A,1);
    nx   =  size(H,1);
    x   =   SX.sym('x',nx,1);
    f   =   0.5*x.'*H*x+g'*x;
    gc   =   A*x-b;
    lb = -10*ones(nx,1);
    ub = 10*ones(nx,1);
    opts.ipopt =struct('max_iter', 250);
    nlp =   struct('x',x,'f',f,'g',gc);
    cas =   nlpsol('solver','ipopt',nlp,opts);
    sol =   cas('lbx', lb,...
                'ubx', ub,...
                'lbg', zeros(neq,1), ...
                'ubg', zeros(neq,1));
    xopt    =   full(sol.x);
    
    lam     =   full(sol.lam_g);  
    
end

%%
if strcmp(solver, 'linsolve')
     
    neq = size(A,1);
    nx  = size(H,1);
    
    LEQS_A  =   [H A';
                 A zeros(neq)];
    LEQS_B  =   [-g; b];
    
    if abs(det(LEQS_A))<10^-100
        Ab = unique([A,b],'rows');
        A=Ab(:,1:size(Ab,2)-1);
        b=Ab(:,size(Ab,2));
        neq = size(A,1);
        LEQS_A  =   [H A';
                    A zeros(neq)];
        LEQS_B  =   [-g; b];
    end

    if abs(det(LEQS_A))<10^-100
    %    keyboard;
     LEQS_x  = pinv(LEQS_A)*LEQS_B; 
    else
     LEQS_x  = linsolve(LEQS_A,LEQS_B);

        if sum(isnan(LEQS_x)) > 0 % regularization if no solution
           LEQS_A    = LEQS_A + abs(min(eig(LEQS_A)))*eye(size(LEQS_A));

           LEQS_x  = linsolve(LEQS_A,LEQS_B);
        end
    end

    xopt    = LEQS_x(1:nx);
    lam     = LEQS_x((nx+1):end); 

end

%%
if strcmp(solver, 'backslash')

    if cond(A) > 10^4
        Ab = unique([A,b],'rows');
        A=Ab(:,1:size(Ab,2)-1);
        b=Ab(:,size(Ab,2));
    end
     
    neq = size(A,1);
    nx  = size(H,1);
    
    LEQS_A  =   [H A';
                 A zeros(neq)];
    LEQS_B  =   [-g; b];
 
    
LEQS_x = LEQS_A\LEQS_B;

    xopt    = LEQS_x(1:nx);
    lam     = LEQS_x((nx+1):end); 

end


%%
if strcmp(solver, 'quadprog')
    % introduced lb,ub, to deal with indefiniteness of H 
    nx   =  size(H,1);
    opt = optimoptions('quadprog','Algorithm','interior-point-convex');
    [xopt,~,~,~,lam_str]=...
        quadprog(H,g,[],[],A,b,-10*ones(nx,1),10*ones(nx,1),[],opt);
    lam  = lam_str.eqlin;
end

%%
if strcmp(solver, 'pinv')
     
    neq = size(A,1);
    nx  = size(H,1);
    
    LEQS_A  =   [H A';
                 A zeros(neq)];
    LEQS_B  =   [-g; b];
 
    %This check is often unnecessary and wastes time if the det is
    %difficult to compute
    if abs(det(LEQS_A))<10^-100

        Ab = unique([A,b],'rows');
        A=Ab(:,1:size(Ab,2)-1);
        b=Ab(:,size(Ab,2));
        neq = size(A,1);
        LEQS_A  =   [H A';
                    A zeros(neq)];
        LEQS_B  =   [-g; b];
    end

    LEQS_x = pinv(LEQS_A)*LEQS_B;

    xopt    = LEQS_x(1:nx);
    lam     = LEQS_x((nx+1):end); 
end
end