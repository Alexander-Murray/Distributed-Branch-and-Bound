function children = getChildren(d,cur_dec_string)
    children = [];
    cur_dec = str2num(cur_dec_string);
    nextUnfixedVar = min(find(cur_dec==0));
    nChildren = d(nextUnfixedVar);
    for i = 1:nChildren
       dec_temp =  cur_dec;
       dec_temp(nextUnfixedVar)=i;
       children = [children; string(num2str(dec_temp))];
    end
end