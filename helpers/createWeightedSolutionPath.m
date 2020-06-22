function weights = createWeightedSolutionPath(sol,ints,s,e)
    if e<s
       error('s must be less than or equal to e') 
    end

    weights = [];
    nInts = length(sol);
    dec = zeros(1,nInts);
    for i = 1:nInts
       % number of decisions for variable i
       d(i) = length(ints{i});
    end
    
    % create weighting for given sol
    for i = s:e   
       dec(i) = find(ints{i}==sol(i));
       weights = [weights;[string(num2str(dec)),4*length(d)-i]];
    end 
    
    % create weighting for children of sol
    nodes = getChildrenv2(d,string(num2str(dec)));
    while ~isempty(nodes)
        current_node = nodes(1);
        weights = [weights;[current_node,3*length(d)]];
        nodes = [nodes; getChildrenv2(d,current_node)];
        nodes(1) = [];
    end
end