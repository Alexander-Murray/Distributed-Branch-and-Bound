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
end