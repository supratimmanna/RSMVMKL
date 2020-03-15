warning off

% load YALE_165n_1024d_15c_zscore_uni1.mat
% load YALE_165n_1024d_15c_zscore_uni_allkernel1.mat

alpha=1e-5;
beta=10;
gamma= 5;
res_R=[];


result=selfweightmkl_Robust(K,y,alpha,beta,gamma)
% for i= 1:10
%     result=selfweightmkl_Robust(K,y,alpha,beta,gamma)
%     res_R=[res_R;result];
% end
