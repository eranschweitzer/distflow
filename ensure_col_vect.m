function x = ensure_col_vect(x)
%%% make sure x is a column vector

m = size(x,1); n = size(x,2);
if (m ~= 1) && (n ~= 1)
    error('input x must be a vector! size = (%d,%d)\n', m, n)
elseif n ~= 1
    x = x.';    
end

end