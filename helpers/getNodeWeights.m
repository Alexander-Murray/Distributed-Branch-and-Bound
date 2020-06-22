function weights = getNodeWeights(d,tree,tree_vals,node_stack,strategy,custom_weights)
    nNodes = length(node_stack); 
    weights = zeros(nNodes,1);
    NS_dub = str2num(char(node_stack));
    nd = length(NS_dub(1,:));
    
    if isstring(custom_weights) %assumes as input custom_weights as n x 2 string array with dec_vec in column 1 and weights in column 2
        cw_numeric = zeros(size(custom_weights,1),nd+1);
        for i = 1:size(custom_weights,1)
            cw_numeric(i,:) = [str2num(custom_weights(i,1)), str2num(custom_weights(i,2))];
        end
        custom_weights = cw_numeric;
    end
    
    % pre-defined bonuses
    if ~isempty(custom_weights)
        for i = 1:nNodes
            di = 1;
            while NS_dub(i,di) ~= 0
                weights(i) = weights(i) + sum(custom_weights(NS_dub(i,di)==custom_weights(:,di),nd+1));
                di = di+1;
                if di > nd
                   break 
                end
            end
        end
    end
    
    if strcmp(strategy,'depth-first')
        for i = 1:nNodes
            weights(i) = weights(i)+sum(NS_dub(i,:)~=0);
        end
    elseif strcmp(strategy,'breadth-first')
        for i = 1:nNodes
            weights(i) = weights(i)+sum(NS_dub(i,:)==0);
        end
    elseif (strcmp(strategy,'best-first'))
        % augment weighting based on best lower bounds
        [Ls,indices] = get_Ls(tree,tree_vals,node_stack,d);
        nInd = length(indices);
        % give boosts according to obj. val.
        for ind = 1:nInd
            children = getChildren(d,tree(indices(ind)));
            for i = 1:length(children)
                weight_ind = find_tree_ind(node_stack,children(i)); % small bonuses given based on lower bounds
                if ~isempty(weight_ind)    
                    if max(Ls)-min(Ls) == 0
                        weights(weight_ind) = weights(weight_ind) + (Ls(ind)-min(Ls))/(max(Ls)-min(Ls)+1);
                    else
                        weights(weight_ind) = weights(weight_ind) + (Ls(ind)-min(Ls))/(max(Ls)-min(Ls));
                    end
                end
            end
        end
    elseif strcmp(strategy,'improve L')
        [Ls,indices] = get_Ls(tree,tree_vals,node_stack,d);
        sorted_Ls = sortrows([Ls,indices]);
        nInd = length(indices);
        minL = min(Ls);
        % give boosts according to obj. val.
        for ind = 1:nInd
            children = getChildren(d,tree(indices(ind)));
            for i = 1:length(children)
                weight_ind = find_tree_ind(node_stack,children(i)); 
                if ~isempty(weight_ind)    
                    % breadth-first on nodes that can improve lower bound
                      weights(weight_ind) = nNodes*sorted_Ls(ind,1)/minL;
                end
            end
        end
    elseif strcmp(strategy,'best bounds')
        % augment weighting based on best lower bounds
        [Ls,indices] = get_Ls(tree,tree_vals,node_stack,d);
        nInd = length(indices);
        best_ind = indices(find(Ls == min(Ls))); % best lower bound gets a bonus
        best_ind = min(best_ind); %in case there are multiple minima
        children = getChildren(d,tree(best_ind));
        % give children of best LB a boost
        for i = 1:length(children)
            weight_ind = find_tree_ind(node_stack,children(i));
            weights(weight_ind) = weights(weight_ind) + length(d) + 1;
        end
        % give boosts according to obj. val.
        for ind = 1:nInd
            children = getChildren(d,tree(indices(ind)));
            for i = 1:length(children)
                weight_ind = find_tree_ind(node_stack,children(i)); % small bonuses given based on lower bounds
                if ~isempty(weight_ind)    
                    if max(Ls)-min(Ls) == 0
                        weights(weight_ind) = weights(weight_ind) + (Ls(ind)-min(Ls))/(max(Ls)-min(Ls)+1);
                    else
                        weights(weight_ind) = weights(weight_ind) + (Ls(ind)-min(Ls))/(max(Ls)-min(Ls));
                    end
                end
            end
        end
    end
 
end