function [finish,status] = terminationCheckBnB(node_stack,U,L,eps)

    if isempty(node_stack)
        finish = 1;
        if U == inf
            disp('Integer decision tree fully searched. No feasible solution detected.')
            status = "INFEASIBLE";
        elseif U-L < eps
            disp("Integer decision tree fully searched.")
           status = "OPTIMAL";
        else
            disp(["Integer decision tree fully searched. Optimality gap: ",U-L])
            status = "SUBOPTIMAL";
        end
    elseif U-L < eps
           finish = 1;
           disp(["Terminated with U-L < ",eps])
           status = "OPTIMAL";
    else
        finish = 0;
        status = "IN PROGRESS. IF YOU ARE SEEING THIS THEN SOMETHING WENT WRONG!";
    end
end