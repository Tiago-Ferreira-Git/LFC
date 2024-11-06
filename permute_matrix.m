function P = permute_matrix(A,rows_to_switch)
%permute_matrix   Given a state-space representation, this function obtains the transformation matrix
%   for informed one-step.
%
%   Inputs:                          
%       
%       A - Global transition matrix. 
%       rows_to_switch - RESs index position.
%
%
%   Outputs:
%
%       P - Tranformation matrix.

    if size(rows_to_switch,1) > size(A,1)
        error 'Number of non-controllable states greater than matrix dimension'
    end
    n_NC = size(rows_to_switch,2);
    n_C = size(A,1) - n_NC;
    

    
    
    rows_NC = n_C+1:size(A,1);
    mask = ismember(rows_to_switch,rows_NC);
    
    %there is no need to switch non controlable modes if they are in the
    %last rows (remove them from the "to switch rows"))
    if any(mask)
        mask_rows = ismember(rows_NC,rows_to_switch);

        rows_to_switch(mask) = [];
        rows_NC(mask_rows) = [];

    end

    P = eye(size(A));
    for i = 1:size(rows_to_switch,2)
        P(rows_to_switch(i),rows_to_switch(i)) = 0; 
        P(rows_to_switch(i),rows_NC(i)) = 1;

        P(rows_NC(i),rows_NC(i)) = 0;
        P(rows_NC(i),rows_to_switch(i)) = 1;

    end

    
    
end

