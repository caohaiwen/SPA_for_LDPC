function [A] = H_sf(A)
%Transform H into the (approximate) standard form (I|P);

[m,n] = size(A);

% Loop over the entire matrix.
i = 1;
j = 1;
jb = [];
zero_rows = [];
num_zeros = 0;
while (i <= m) && (j <= n)
   % Find value and index of largest element in the remainder of column j.
   [p,k] = max(A(i:m,j));   
   if (p == 0) 
       num_zeros = num_zeros + 1;
       zero_rows(num_zeros) = i;
   end
   k = k+i-1;
   % Remember column index
   jb = [jb j];
   % Swap i-th and k-th rows.
   A([i k],j:n) = A([k i],j:n);
   % Divide the pivot row by the pivot element.
   A(i,j:n) = A(i,j:n);
   % Subtract multiples of the pivot row from all the other rows.
   for k = [1:i-1 i+1:m]
        A(k,j:n) = mod(A(k,j:n) + A(k,j)*A(i,j:n), 2);
   end
   i = i + 1;
   j = j + 1;
end

for num = 1 : num_zeros
    current_row = zero_rows(num);
    [pp, kk] = max(A(current_row, :));
    if (pp == 1)
        for k = [1:current_row-1 current_row+1:m]
            A(k,kk:n) = mod(A(k,kk:n) + A(k,kk)*A(current_row,kk:n), 2);
        end
    end
end

    


