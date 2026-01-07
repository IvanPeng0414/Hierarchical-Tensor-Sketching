function [samples,normal_const] = tt_sample(TT, num_sample)

% Conditional sampling for a density TT, TT is compatible with TTtoolbox
% Since the rank cannot be 1, we add a tt_zeros if needed
if min(TT.r) == 1
    TT = TT + tt_zeros(TT.n,TT.d);
end

% Precomputation
samples = zeros(num_sample, TT.d);

integrated = cell(1,TT.d);
for i=1:TT.d
    integrated{i} = squeeze(sum(TT{i},2));

    if i==1
        integrated{i} = transpose(squeeze(sum(TT{i},2)));
    end
end
prod_integrated_from_right = cell(1,TT.d);
res = 1;
for i=1:TT.d
    res = integrated{numel(integrated)-i+1} * res;
    prod_integrated_from_right{i} = res;
end
prod_integrated_from_right = fliplr(prod_integrated_from_right);
normal_const = prod_integrated_from_right{1};

for sample_id = 1:num_sample
    tmp_sample = zeros(1, TT.d);
    left       = 1;

    for i=1:TT.d
        if i==1
            density = squeeze(TT{i}) * prod_integrated_from_right{i+1};
            sp      = randsample(numel(density),1,true,density);
            tmp_sample(i) = sp;

            left    = left * squeeze(TT{i}(:,sp,:));
            left    = left.';
        elseif i==TT.d
            density = squeeze(tns_mult(left,2,TT{i},1));
            sp      = randsample(numel(density),1,true,density);
            tmp_sample(i) = sp;
        else
            density = squeeze(tns_mult(left,2,TT{i},1)) * prod_integrated_from_right{i+1};
            sp      = randsample(numel(density),1,true,density);
            tmp_sample(i) = sp;

            left    = left * squeeze(TT{i}(:,sp,:));
        end
    end

    samples(sample_id,:) = tmp_sample;
end

end