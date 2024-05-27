function [W,ss] = permute_matrix(A,rows_to_switch)
    

    
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
    ss = 1:(n_C+n_NC);
    ss(:,rows_to_switch') = rows_NC;
    ss(:,rows_NC) = rows_to_switch;
    % 
    % rows_to_switch = sort(rows_to_switch);
    % rows_NC = sort(rows_NC);

    W = eye(size(A));
    for i = 1:size(rows_to_switch,2)
        %W(rows_to_switch(i),:) = zeros(1,n_C+n_NC);
        W(rows_to_switch(i),rows_to_switch(i)) = 0; 
        W(rows_to_switch(i),rows_NC(i)) = 1;

        W(rows_NC(i),rows_NC(i)) = 0;
        W(rows_NC(i),rows_to_switch(i)) = 1;

    end

    W = W;

    % [ss_sorted,I] = sort(ss(1:n_C));
    % if ~all(all(ss_sorted == ss(1:n_C)))
    %     W_sort = eye(size(A));
    %     W_sort(1:n_C,1:n_C) = zeros(n_C);
    % 
    % 
    %     for i = 1:n_C
    %         W_sort(I(i),i) = 1;
    %         %W_sort(I(i,2),I(i,1)) = 1;
    %     end
    %     W = W_sort\W; 
    % end
    % ss(1:n_C) = ss_sorted;
    % 
    % [ss_sorted,I] = sort(ss(n_C+1:end));
    % if ~all(all(ss_sorted == ss(n_C+1:end)))
    %     I = I + n_C;
    %     W_sort = eye(size(A));
    % 
    %     W_sort(n_C+1:end,n_C+1:end) = zeros(n_NC);
    %     for i = 1:size(I,1)
    %         W_sort(I(i),i) = 1;
    %     end
    %     W = W_sort\W;  
    % end
    % ss(n_C+1:end) = ss_sorted;
    
    
    
end

