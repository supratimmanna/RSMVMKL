alpha_all=[1e-6 1e-5 1e-4 1e-3];
gam_all = [1 3 5 8 10 12 15 18 20 25];
% r_all = [.1 .3 .5];
beta_all=[1000 2000];
% r=0.1;%rate of labeled data
res=[];
plot=[];
plot20=[];


for data=1:3
    Kernel = KK{1,data+1};
    y=yy{1,data+1};
    for b=1:2
        beta=beta_all(b);
        
        
        for ss=1:4
            alpha = alpha_all(ss);
            for tt=1:10
                gam=gam_all(tt);
                res=[];
                for i=1:2
                    
                    result=selfweightmkl_MV_Robust(Kernel,y,alpha,beta,gam)
                    res=[res;result];
                    
                end
                mn2=mean(res);
%                 res=sort(res,'descend');
%                 mn10=mean(res(1:10,:));
                maxacc=max(res(:,1));
                plot=[plot;alpha gam maxacc mn2 beta data+1];
%                 plot20=[plot20;alpha gam mn2 beta];
            end
        end
    end
end

save('ACC_NMI_PUR_Iter_AllData_Clustering_Robust.mat','plot');
