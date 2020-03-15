warning off

alpha=1e-3;
beta= 1000; 
gamma=15; 
res_R=[];

%%% for best result make L=D-S%%%%%%%%%
result=selfweightmkl_MV_Robust(Kernel,y,alpha,beta,gamma)
res_R=[res_R;result];
for i= 1:10
    result=selfweightmkl_MV_Robust(Kernel,y,alpha,beta,gamma)
    res_R=[res_R;result];
end
