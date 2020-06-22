function tree_ind = find_tree_ind(tree,cur_dec)
    tree_size = length(tree);
    tree_ind = [];
    
    for i = 1:tree_size
       if strcmp(tree(i),cur_dec)
           tree_ind = i;
       end
    end
end