%% add TT-Toolbox to the path
addpath('~/TT-Toolbox-master/TT-Toolbox-master');

%% parameters for Ising model
L    = 4;
n    = 2;
d    = 2^L;
d1   = sqrt(d);
beta = 0.6; eps = 3e-2;
tol  = 1*10^-6;

%% binary tree partition
build_family_tree_ising2d;

%% ground truth distribution

J1 = 1.0; J2 = 0.0;
den_fxn = @(ind) exp(-beta*ising2dden(J1,J2,d1,n,ind));

%% gain samples
num_sample = 64000; 
den_fxn2 = @(ind) exp(-beta/2*ising2dden(J1,J2,d1,n,ind));
[samples_aug] = dmrg_sample(d1^2,n,den_fxn2,num_sample,eps);
Plist = samples_aug(:,end);

p_biindex = samples_aug(:,1:d1^2);

%% set rank
target_rank = cell(L,d);
for l = 1 : L 
    for k=1:numel(H(l).cluster_idx)   
        target_rank{l,k} = 3; 
    end
end


exact_rank = cell(L,d);

PF   = cell(L,d); PFindex = cell(L,d);
A    = cell(L,d);
B    = cell(L,d);
dtotalcell = cell(L,d);

%% compute core function
for l=1:L
    l
    if l==1
        % set top level first, the core is a 2-d matrix.
        nc = numel(H(l).cluster_idx);        
        
        Utotal = cell(nc,1); indtotal = cell(nc,1);
        for k=1:nc
           
            din = numel(H(l).cluster_idx{k});
            dout = 2^(L-l+1)-numel(H(l).cluster_idx{k});  
            dinindex = sort(reshape(H(l).cluster_idx{k},[din,1])); indtotal{k} = dinindex;
            dtotalcell{l,k} = dinindex;
            doutindex = sort(reshape(setdiff(1:d, H(l).cluster_idx{k}),[dout,1]));
            doutindexcell = cell(1,1); doutindexcell{1} = doutindex;
            % S is sketch1*P*sketch2
            [S] = sparse_total_polysketch(Plist,p_biindex,dinindex,doutindex);
            % do randomized SVD
            [U,s,V] = rsvd(S,2*target_rank{l,k},1,target_rank{l,k});               
            idx      = find(diag(s)/s(1,1)>tol);
            if numel(idx)>target_rank{l,k}
                idx = idx(1:target_rank{l,k});
            end
            
            V = V(:,idx); 
            U = U(:,idx); 
            Utotal{k} = U;
            r = numel(idx);  exact_rank{l,k} = r;                     
            % compute P_k^l, for next level
            [PF{l,k},PFindex{l,k}] = PF_sparse_right_polysketch(Plist,p_biindex,V,dinindex,doutindex);
            % compute A
            A{l,k}    = U'*S*V;  
                             
        end
        % compute B
        dtotalindex = [1:1:d];
        [B0] = B0_sparse_polysketch(Plist,p_biindex,Utotal{1},Utotal{2},indtotal{1},indtotal{2});

    else
        % compute A_k^l and B_k^l in the following levels
        ct = 1;
        for k=1:numel(H(l-1).cluster_idx)
            
            if ~isempty(H(l-1).cluster_idx{k})
                nc = numel(H(l-1).child{k});
                
                Utotal = cell(nc,1); indtotal = cell(nc,1);
                
                for kk = 1 : nc

               
                    din = numel(H(l).cluster_idx{H(l-1).child{k}(kk)});
                    dout = 2^L - din;
                    
                    dinindex = sort(reshape(H(l).cluster_idx{H(l-1).child{k}(kk)},[din,1])); indtotal{kk} = dinindex;
                    
                    doutindex = sort(reshape(setdiff(1:d,H(l).cluster_idx{H(l-1).child{k}(kk)}),[dout,1]));

                    doutindexcell = cell(l,1);
                    doutindexcell{l} = setdiff(H(l-1).cluster_idx{k}, H(l).cluster_idx{H(l-1).child{k}(kk)});
                    douttmp = [doutindexcell{l}(:)'];
                    if l >= 3
                        for lindex = l-1 : -1 : 2
                            doutindexcell{lindex} = setdiff(setdiff(H(lindex-1).cluster_idx{ceil(k/2^(l-lindex))},H(l).cluster_idx{H(l-1).child{k}(kk)}), douttmp);
                            douttmp = [douttmp,doutindexcell{lindex}(:)'];
                        end
                    end
                    doutindexcell{1} = setdiff(setdiff([1:1:d]',H(l).cluster_idx{H(l-1).child{k}(kk)}),douttmp);
                    % compute sketch1*P*sketch2                   
                    [S] = sparse_total_polysketch(Plist,p_biindex,dinindex,doutindex);
                    [U,s,V]  = rsvd(S,2*target_rank{l,k},1,target_rank{l,k}); 

                    idx      = find(diag(s)/s(1,1)>tol);
                    if numel(idx)>target_rank{l,k}
                        idx = idx(1:target_rank{l,k});
                    end
                    % to avoid the 1dim to make the singularity
                    if numel(idx) < 2
                        idx = [1;2];
                    end                    
                    V = V(:,idx);
                    U = U(:,idx);  Utotal{kk} = U;
                    r = numel(idx); 
                    % compute P_k^l for the next level
                    [PF{l,ct},PFindex{l,ct}] = PF_sparse_right_polysketch(Plist,p_biindex,V,dinindex,doutindex);

                    exact_rank{l,ct} = r;  dtotalcell{l,ct} = dinindex; 
                    % compute A_k^l
                    A{l,H(l-1).child{k}(kk)}   = U'*S*V; 

                    ct = ct+1;                    
                end
                    
                % compute B_k^l
                ptmp = PF{l-1,k}; r = size(ptmp,2);
                dtotalindex = dtotalcell{l-1,k}; 
                p_biindex_new = PFindex{l-1,k};

                B0tmp = B0_sparse_polysketch(ptmp(:,1),p_biindex_new,Utotal{1},Utotal{2},indtotal{1},indtotal{2});
                Btotal = zeros(size(B0tmp,1),size(B0tmp,2),exact_rank{l-1,k});
                Btotal(:,:,1) = B0tmp;
                for rind = 2 : exact_rank{l-1,k}
                    B0tmp = B0_sparse_polysketch(ptmp(:,rind),p_biindex_new,Utotal{1},Utotal{2},indtotal{1},indtotal{2});
                    Btotal(:,:,rind) = B0tmp;
                end
                B{l-1,k} = Btotal;
                                
            end
        end
        
    end
        
   
end



%% compute cores C=A_1^+ B A_2^+

C = cell(L,d);
for l=1:L
    
    
    if l==1
        % compute the top level
        C0 = pinv(A{1,1})*B0*pinv(A{1,2}');
        
    else
        % compute the following levels     
        for k = 1:numel(H(l-1).child)
            
            if ~isempty(H(l-1).child{k})
                
                A1   = A{l,H(l-1).child{k}(1)};
                
                A2   = A{l,H(l-1).child{k}(2)};
                
                b    = B{l-1,k};
                
                tmp1 = pinv(A1);
                tmp2 = pinv(A2);
                
                C{l-1,k}   = tns_mult(tmp1,2,b,1);
                C{l-1,k}   = permute(tns_mult(tmp2,2,C{l-1,k},2),[2 1 3]);
                
 
            end
        end
    
    end
end
for k = 1 : d
    C{L,k} = PF{L,k}; 
end
    
%% compute error

tmp      = cell(1,d);
[tmp{:}] = ind2sub(ones(1,d)*2,1:2^d);
all_samples = zeros(2^d,d+1);
for i=1:numel(tmp)
    all_samples(:,i) = tmp{i};
end
for i=1:2^d
    all_samples(i,end) = den_fxn(all_samples(i,1:d));
end
all_samples(:,end) = all_samples(:,end)/sum(all_samples(:,end));

Ptrue = all_samples(:,end);

%% represent the density via the cores
Precon_full = zeros(2^d,1);
for i = 1 : 2^d
    sample = all_samples(i,1:d); 
    [y] = apply_tree_one_sample_binary2d_final(C0,C,sample,H,L);
    Precon_full(i) = y;
end
fprintf("error with full true distribution %.8f.\n",norm(Precon_full-Ptrue)/norm(Ptrue))


