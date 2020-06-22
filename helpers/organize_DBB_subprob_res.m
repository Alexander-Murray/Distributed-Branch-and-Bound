function [x_sol,z_sol,init,lam0,subprob_feas] = organize_DBB_subprob_res(subalg,sub_sol,nlp_log,g_glob,A,B,c,lbx,lbz,ubx,ubz)
    N = length(lbx);
    for i = 1:N
       rvars(i) = length(lbx{i});
       ivars(i) = length(lbz{i});
    end
    
    if strcmp(subalg,'ALADIN')
        for n = 1:N
            loc_x_cur_decex = 1+sum(rvars(1:n-1))+sum(ivars(1:n-1)):sum(rvars(1:n-1))+sum(ivars(1:n-1))+rvars(n);
            loc_z_cur_decex = 1+sum(rvars(1:n))+sum(ivars(1:n-1)):sum(rvars(1:n))+sum(ivars(1:n-1))+ivars(n);
            x_sol{n} = full(sub_sol(loc_x_cur_decex));
            z_sol{n} = full(sub_sol(loc_z_cur_decex));
        end
        lam0 = nlp_log.lambda(:,end); % warmstart lambda using prev sol
        if strcmp(nlp_log.status,'timeout')
            proj_feas = check_feas_MIP(g_glob,horzcat(A{:}),horzcat(B{:}),c,vertcat(lbx{:}),vertcat(lbz{:}),vertcat(ubx{:}),vertcat(ubz{:}),vertcat(x_sol{:}),round(vertcat(z_sol{:})));
            if proj_feas == 1
                subprob_feas = 1;
                disp("Algorithm timeout. Projection onto integers is feasible.")
            else
                subprob_feas = 0;
                disp("NLP failed to converge!")
            end
        else
           subprob_feas = 1; 
        end
    elseif or(strcmp(subalg,'ADMM'),strcmp(subalg,'ADMM2'))
        for n = 1:N
            loc_x_cur_decex = 1+sum(rvars(1:n-1))+sum(ivars(1:n-1)):sum(rvars(1:n-1))+sum(ivars(1:n-1))+rvars(n);
            loc_z_cur_decex = 1+sum(rvars(1:n))+sum(ivars(1:n-1)):sum(rvars(1:n))+sum(ivars(1:n-1))+ivars(n);
            x_sol{n} = full(sub_sol(loc_x_cur_decex)); % warmstart lambda using prev sol
            z_sol{n} = full(sub_sol(loc_z_cur_decex)); % warmstart lambda using prev sol
            if strcmp(subalg,'ADMM')
                lam0{n} = zeros(length(c),1); 
            elseif strcmp(subalg,'ADMM2')
                lam0{n} = nlp_log.lambda(sum(rvars(1:n-1)+ivars(1:n-1))+1:sum(rvars(1:n)+ivars(1:n)),end); % warmstart lambda using prev sol
            end
        end
        if strcmp(nlp_log.status,'timeout')
            proj_feas = check_feas_MIPv2(10^-3,g_glob,horzcat(A{:}),horzcat(B{:}),c,vertcat(lbx{:}),vertcat(lbz{:}),vertcat(ubx{:}),vertcat(ubz{:}),vertcat(x_sol{:}),vertcat(z_sol{:}));
            if proj_feas == 1
                subprob_feas = 1;
                disp("Algorithm timeout. NLP solution is feasible.")
            else
                subprob_feas = 0;
                disp("NLP failed to converge!")
            end
        else
           subprob_feas = 1; 
        end
   elseif strcmp(subalg,'IPOPT')
        lam0 = 0;
        subprob_feas = 1;
        for n = 1:N
            loc_x_cur_decex = 1+sum(rvars(1:n-1)):sum(rvars(1:n-1))+rvars(n);
            loc_z_cur_decex = 1+sum(ivars(1:n-1)):sum(ivars(1:n-1))+ivars(n);
            x_sol{n} = full(sub_sol(loc_x_cur_decex));
            z_sol{n} = full(sub_sol(sum(rvars)+loc_z_cur_decex));
        end
    elseif strcmp(subalg,'ADMM_QP')
        lam0 = 0;
        subprob_feas = 1;
        for n = 1:N
            loc_x_cur_decex = 1+sum(rvars(1:n-1))+sum(ivars(1:n-1)):sum(rvars(1:n-1))+sum(ivars(1:n-1))+rvars(n);
            loc_z_cur_decex = 1+sum(rvars(1:n))+sum(ivars(1:n-1)):sum(rvars(1:n))+sum(ivars(1:n-1))+ivars(n);
            x_sol{n} = full(sub_sol(loc_x_cur_decex)); % warmstart lambda using prev sol
            z_sol{n} = full(sub_sol(loc_z_cur_decex)); % warmstart lambda using prev sol
        end
    elseif strcmp(subalg,'ADMM_exp_lasso')
        lam0 = 0;
        subprob_feas = 1;
        for n = 1:N
            loc_x_cur_decex = 1+sum(rvars(1:n-1))+sum(ivars(1:n-1)):sum(rvars(1:n-1))+sum(ivars(1:n-1))+rvars(n);
            loc_z_cur_decex = 1+sum(rvars(1:n))+sum(ivars(1:n-1)):sum(rvars(1:n))+sum(ivars(1:n-1))+ivars(n);
            x_sol{n} = full(sub_sol(loc_x_cur_decex)); % warmstart lambda using prev sol
            z_sol{n} = full(sub_sol(loc_z_cur_decex)); % warmstart lambda using prev sol
        end
    end
    for n = 1:N
       init{n} = [x_sol{n};z_sol{n}]; % warmstart using prev sol (improve this by checking node value, relaxed ints, etc.) 
    end
    
end