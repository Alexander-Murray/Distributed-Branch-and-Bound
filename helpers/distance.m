function d = distance(z,lbz,ubz)
    ind1 = z<lbz;
    ind2 = z>ubz;
    d = sum(abs(z(ind1)-lbz(ind1))) + sum(abs(z(ind2)-ubz(ind2)));
end