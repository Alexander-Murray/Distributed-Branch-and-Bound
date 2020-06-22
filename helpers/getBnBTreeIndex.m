function ind = getBnBTreeIndex(d,x)
    % d is the vector whose entries are the number of decisions
    % corresponding to each integer variable.
    % x is the index of the chosen integer decision. A value of 0 means
    % that no value is chosen.
    % ind is the index of the array which stores the nodes of the Branch
    % and Bound decision tree
    
    if length(d)~=length(x)
       error("d and x must be the same length!") 
    end
    
    ind = uint64(x(1)+1);
    p = uint64(1);
    i = 2;
    x_i = x(i);
    while and(x_i~=0, i <=length(d))
        x_i = x(i);

        p = p*(d(i-1)+1);

        ind = ind + x_i*p;

        i = i+1;
    end
end