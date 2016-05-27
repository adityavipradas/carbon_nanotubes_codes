function istrue = any_equal(A,B)
%Outputs Boolean 1 ("True") if any element in 1D Array A is equal in value to any element in 1D Array B (regardless of position). 
%(only to be used for 1D arrays
istrue = 0;

for i = 1:length(A)
    for j = 1:length(B)
        if A(i) == B(j)
            istrue = 1;
        else
        end
    end
end

        