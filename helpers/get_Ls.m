% gets lower bounds from nodes with unexplored children
function [Ls,indices] = get_Ls(tree,tree_vals,node_stack,d)
   nExplored_nodes = length(tree);
   Ls = [];
   indices = [];
   for i = 1:nExplored_nodes
       children = getChildren(d,tree(i));
       for j = 1:length(children)
           if ~isempty(find_tree_ind(node_stack,children(j))) % check if node has unexplored children
                Ls = [Ls;tree_vals(i)];
                indices = [indices; i];
                break
           end
       end
   end
end