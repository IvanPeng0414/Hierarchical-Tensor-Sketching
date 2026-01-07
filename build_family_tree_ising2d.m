% define relationship tree

[tmp1,tmp2] = meshgrid(1:d1,1:d1);
grid_pos_old    = [tmp1(:) tmp2(:)];
grid_pos = grid_pos_old;

for i = 1 : floor(d1/2)
    grid_pos((2*i-1)*d1+1:2*i*d1,:) = flip(grid_pos_old((2*i-1)*d1+1:2*i*d1,:),1);
end
    
H(1).cluster_pos{1} = grid_pos(grid_pos(:,1)<d1/2+0.01,:);
H(1).cluster_pos{2} = grid_pos(grid_pos(:,1)>d1/2+0.01,:);

H(1).cluster_idx{1} = find_subset(H(1).cluster_pos{1}, grid_pos);
H(1).cluster_idx{2} = find_subset(H(1).cluster_pos{2}, grid_pos);


for l=1:L-1
    
    if mod(l,2)==1
        
        nc = numel(H(l).cluster_pos);
        
        ct = 1;
        for k=1:nc  %follow each direction separately
            idx_curr = H(l).cluster_pos{k};
            i1 = max(idx_curr(:,2));
            i2 = min(idx_curr(:,2));
            
            mid    = i2+(i1-i2)/2;
            
            H(l+1).cluster_pos{ct} = idx_curr(idx_curr(:,2)<mid,:);
            H(l+1).cluster_idx{ct} = find_subset(H(l+1).cluster_pos{ct}, grid_pos);
            
            ct = ct+1;
            H(l+1).cluster_pos{ct} = idx_curr(idx_curr(:,2)>mid,:);
            H(l+1).cluster_idx{ct} = find_subset(H(l+1).cluster_pos{ct}, grid_pos);
            
            ct = ct+1;
            
        end
        
    else
        
        nc = numel(H(l).cluster_pos);
        
        ct = 1;
        for k=1:nc
            idx_curr = H(l).cluster_pos{k};
            i1 = max(idx_curr(:,1));
            i2 = min(idx_curr(:,1));
            mid    = i2+(i1-i2)/2;
            
            H(l+1).cluster_pos{ct} = idx_curr(idx_curr(:,1)<mid,:);
            H(l+1).cluster_idx{ct} = find_subset(H(l+1).cluster_pos{ct}, grid_pos);
            
            ct = ct+1;
            H(l+1).cluster_pos{ct} = idx_curr(idx_curr(:,1)>mid,:);
            H(l+1).cluster_idx{ct} = find_subset(H(l+1).cluster_pos{ct}, grid_pos);
            
            ct = ct+1;
            
        end
        
        
    end
    
end




%parent-child relation

for l=1:L-1
   
   p1 = H(l).cluster_idx;
   p2 = H(l+1).cluster_idx;
   
   H(l).child = {};
   H(l+1).parent = {};
   for i=1:numel(p1)
       ch_tmp = [];
       for j=1:numel(p2)
           
           tmp = find_subset(p2{j},p1{i});
           if ~isempty(tmp)
               
              ch_tmp = [ch_tmp j];
              H(l+1).parent{j} = i;
              
           end
           
       end
       
        H(l).child =  [H(l).child ch_tmp];
   end
   
    
end


%redefine tree
for l=fliplr(1:L-1)  
    for k=1:2^l
        tmp = [];
        for kk=1:numel(H(l).child{k})
            tmp = [tmp; H(l+1).cluster_idx{H(l).child{k}(kk)}(:)];
        end
        %tmp'
        H(l).cluster_idx{k} = tmp;
    end
end



function [idx] = find_subset(A,B)
idx = [];
for i=1:size(A,1)
   idx = [idx; find(sum(abs(A(i,:)-B),2)==0)];
end
end