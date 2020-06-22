function node_stack = pruneFromU(U,d,tree,tree_vals,node_stack)
    n = length(node_stack);
    prune = zeros(n,1);
    for i = 1:n
       parent = getParent(d,str2num(node_stack(i))); 
       tree_ind = find_tree_ind(tree,parent);
       if tree_vals(tree_ind)>U
           prune(i) = 1;
       end
    end
    prune = logical(prune);
    node_stack(prune,:) = [];
end