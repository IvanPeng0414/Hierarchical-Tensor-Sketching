function [y] = apply_tree_one_sample_binary2d_final(C0,C,sample,H,L)

n = 2;
d = n^L; 
U = cell(L,d);
for l=fliplr(1:L)
    
    if l==L       
        for k=1:d            
            U{l,k} = C{l,k}(sample(:,H(L).cluster_idx{k}),:);            
        end        
    else        
        for k = 1:numel(H(l).child)
            
            if ~isempty(H(l).child{k})
                i1 = H(l).child{k}(1);
                i2 = H(l).child{k}(2);
                
                
                A1   = U{l+1,i1};
                A2   = U{l+1,i2};
              
                tmp = tns_mult(A1,2,C{l,k},1);
                tmp = permute(tns_mult(A2,2,tmp,2),[2 1 3]);              
                %U{l,k} = reshape(tmp,[],n^(d/2^(l)));
                U{l,k} = reshape(tmp,[],numel(tmp) );
            end
            
         end
        
    end
        
       
end
y = U{1,1}*C0*U{1,2}';

end