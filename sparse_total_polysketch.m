function [S] = sparse_total_polysketch(P,p_biindex,dinindex,doutindex)
% the code is to compute sketch1*P*sketch2 in an efficient way

num_effsample = size(P,1);
Pcolumn = p_biindex(:,doutindex) - 1; 

Pcolumn_unique = unique(Pcolumn,'rows');

indexlen_column = size(Pcolumn_unique,1);


Prow = p_biindex(:,dinindex) - 1; %compute corresponding index
Prow_unique = unique(Prow,'rows');


if size(Prow_unique,1) == 1
    Prow_unique = [Prow_unique;1-Prow_unique];
end


indexlen_row = size(Prow_unique,1);


%% Kleft (sketch matrix)
din = length(dinindex);
left_dim = din;
poly_order = 4;
if left_dim == 1
    % epss = 1e-6;
    % sketchl = [1,1-epss;1,1-2*epss];
    sketchl = [1,0;0,1];
else
    cnt = 0;
    for i = 0:poly_order
        if left_dim >= i 
            cnt = cnt + nchoosek(left_dim, i);
        end
    end
    % Construct sketch
    sketchl = zeros(indexlen_row, cnt);
    count = 1;
    for i = 0:poly_order
        if i == 0
            sketchl(:, 1) = 1;
            count = count + 1;
            continue
        end
        tmp = nchoosek((1:left_dim), i);
        for j = 1:size(tmp, 1)
            tmp_sketch = ones(indexlen_row, 1);
            for k = 1:i
                tmp_sketch = tmp_sketch .* (2*Prow_unique(:, tmp(j, k))-1);
            end
            sketchl(:, count) = tmp_sketch;
            count = count + 1;
        end
    end

end
Kleft = sketchl;










%% Kright
dout = length(doutindex);
right_dim = dout;
poly_order = 4;
if right_dim == 1
    % epss = 1e-6;
    % sketchr = [1,1-epss;1,1-2*epss];
    sketchr = [1,0;0,1];
else
    cnt = 0;
    for i = 0:poly_order
        if right_dim >= i 
            cnt = cnt + nchoosek(right_dim, i);
        end
    end
    % Construct sketch
    sketchr = zeros(indexlen_column, cnt);
    count = 1;
    for i = 0:poly_order
        if i == 0
            sketchr(:, 1) = 1;
            count = count + 1;
            continue
        end
        tmp = nchoosek((1:right_dim), i);
        for j = 1:size(tmp, 1)
            tmp_sketch = ones(indexlen_column, 1);
            for k = 1:i
                tmp_sketch = tmp_sketch .* (2*Pcolumn_unique(:, tmp(j, k))-1);
            end
            sketchr(:, count) = tmp_sketch;
            count = count + 1;
        end
    end

end
Kright = sketchr;

%% Pcen
Pcen = zeros(indexlen_row,indexlen_column);

for pind = 1 : num_effsample     
    pind_row = p_biindex(pind,dinindex) - 1; 
    pind_column = p_biindex(pind,doutindex) -1;

    index1 = find(sum(abs(Prow_unique-pind_row),2)==0);
    index2 = find(sum(abs(Pcolumn_unique-pind_column),2)==0);
    Pcen(index1,index2) = P(pind);

end


S = (Kleft'*Pcen)*Kright;


end