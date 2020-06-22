function [ints,int_parts] = get_int_parts(lbz,ubz)
    lbz_glob = vertcat(lbz{:});
    ubz_glob = vertcat(ubz{:});
    int_parts = [];
    p_count = 1;
    v_count = 1;
    for i = 1:length(lbz_glob)
       ints{i} = lbz_glob(i):ubz_glob(i); 
       int_parts = [int_parts;[i,p_count,v_count]];
       if i==length(vertcat(lbz{1:p_count}))
          p_count = p_count + 1;
          v_count = 1;
       else
           v_count = v_count + 1;
       end
    end
end