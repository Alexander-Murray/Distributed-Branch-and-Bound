function mat = remove_rows_of_zeros(mat)
    rows = size(mat,1);
    rsum = zeros(rows,1);
    for i = 1:rows
       rsum(i) = sum(abs(mat(i,:))); 
    end
    mat = mat(rsum>0,:);
end