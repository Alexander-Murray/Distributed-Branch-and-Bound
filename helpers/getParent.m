% no node indices
function parent = getParent(d,cur_dec)
    lastFixedVar = max(find(cur_dec~=0));
    parent_dec = cur_dec;
    parent_dec(lastFixedVar) = 0;
    parent = getBnBTreeIndex(d,parent_dec);
end