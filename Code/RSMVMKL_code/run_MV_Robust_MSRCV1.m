warning off

alpha=1e-5;
beta=1000;
gamma=5; %2.5;
res_R=[];


result=selfweightmkl_MV_Robust(Kernel,y,alpha,beta,gamma)
% res_R=[res_R;result];
% for i= 1:10
%     result=selfweightmkl_MV_Robust(Kernel,y,alpha,beta,gamma)
%     res_R=[res_R;result];
% end
