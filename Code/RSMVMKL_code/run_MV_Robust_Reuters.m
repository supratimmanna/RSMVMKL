warning off

alpha=1e-4;
beta= 100; %100;%1000;
gamma=5; %5;%.92;
res_R=[];


result=selfweightmkl_MV_Robust(Kernel,y,alpha,beta,gamma)
res_R=[res_R;result];
% for i= 1:10
%     result=selfweightmkl_MV_Robust(Kernel,y,alpha,beta,gamma)
%     res_R=[res_R;result];
% end
