function [A,new_indices] = rrefFixer(B)
    [rows,cols] = size(B);
    index = 1;
    ident_cols = zeros(1,cols);
    
    for c = 1:cols
       if and(sum(B(:,c)==1)==1,sum(B(:,c))==1) %exactly one one in the column
           if B(index,c) == 1 %skip repeated columns
              ident_cols(c) = 1; 
              index = index +1;
           end
       end
       if index > rows
          break 
       end
    end
    A = [B(:,logical(ident_cols)), B(:,~ident_cols)];
    new_indices = [find(ident_cols),find(~ident_cols)];
end