function [samples,sample_tmp,normal_const] = dmrg_sample(d,n,den_fxn,num_sample,eps)

    [ytmp] = dmrg_cross(d,n,den_fxn,eps);
    y = ytmp.*ytmp;
    %[y]= round(y,1e-6);
    [sample_tmp,normal_const] = tt_sample(y, num_sample);
    %update the same one and update the weights.
    outputs = containers.Map('KeyType','char','ValueType','double');
    for index = 1 : num_sample
        sample_code = mat2str(sample_tmp(index,:));
        if isKey(outputs, sample_code)
            outputs(sample_code) = outputs(sample_code) + 1/num_sample;
        else
            outputs(sample_code) = 1/num_sample;
        end    
    end
    mykeys = keys(outputs);
    samples = zeros(numel(mykeys),d+1);
    for i=1:numel(mykeys)
        out = str2num(mykeys{i});
        samples(i,:) = [out outputs(mykeys{i})];
    end
    
end