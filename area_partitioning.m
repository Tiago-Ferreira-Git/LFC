function areas = area_partitioning(lines,k,machine_buses)
    %rng(1);
    if size(lines,1) ~= size(unique(lines(:,1:2),'rows'),1)
        %error 'Multiple Defined lines'
        lines = unique(lines,'rows');
    end
    n_lines = size(lines,1);
    
    W = zeros(max([lines(:,1);lines(:,2)]));
    D = zeros(max([lines(:,1);lines(:,2)]));
    for i = 1:n_lines
        W(lines(i,1) , lines(i,2)) = 1;
        W(lines(i,2) , lines(i,1)) = 1;


        D(lines(i,1),lines(i,1)) = D(lines(i,1),lines(i,1)) + 1;
        D(lines(i,2),lines(i,2)) = D(lines(i,2),lines(i,2)) + 1;
    end

    L = D-W;
    
    [V,Eig] = eig(L,"vector");
    %d = eigs(L,3,'lr');
    [Eig,I] = sort(Eig);
    V = V(:,I);

    j = 0;
    flag = false;
    while ~flag
        [areas,C] = kmeans(V(:,1:k),k,'Distance','cosine');
        for i=1:k
            area_bus = find(areas==i);
            if ~any(ismember(area_bus,machine_buses))
                flag = true;
                %warning 'At least one area without machines'
                break
                %error 'At least one area without machines'
            end
        end

        if flag == false
            break
        end
        flag = false;
        j = j +1;
        if rem(j,1000) == 0
            j
        end
    end


    
end

