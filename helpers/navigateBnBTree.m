function [finish,dec_ind,lbz,ubz,status] = navigateBnBTree(node_info,U,L,d,lbz,ubz,ints,int_parts,eps)
    node_stack = node_info(:,1);
    nodal_weights = node_info(:,2);
    [finish,status] = terminationCheckBnB(node_stack,U,L,eps);
    ind = node_stack(find(str2double(nodal_weights)==max(str2double(nodal_weights))));
    dec_ind = ind(1); % in case of ties
    disp(['node: ',dec_ind]);
    
    dec_num = str2num(char(dec_ind));
    
    n = length(d);
    for i = 1:n %loop over all integer decisions
        if dec_num(i)~=0 % check if integer var must be fixed
           lbz{int_parts(i,2)}(int_parts(i,3)) = ints{i}(dec_num(i));
           ubz{int_parts(i,2)}(int_parts(i,3)) = ints{i}(dec_num(i));
        end 
    end
end