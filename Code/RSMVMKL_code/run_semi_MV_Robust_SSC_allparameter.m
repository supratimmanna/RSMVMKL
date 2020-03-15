alpha_all=[1e-6 1e-5 1e-4 1e-3];
gam_all = [20 25];
% r_all = [.1 .3 .5];
beta_all=[1000 2000];
% r=0.1;%rate of labeled data
res=[];
plot=[];
plot20=[];


for b=1:2
    beta=beta_all(b);
    
    
    for ss=1:4
        alpha = alpha_all(ss);
        for tt=1:2
            gam=gam_all(tt);
            res=[];
            for ii=1:30
                r=0.3;
                [m,n,rr,tt]=size(Kernel);
                c=length(unique(y)); % number of class
                numperc=floor(n/c); % number of data per class
                labelperc=floor(r*numperc); % number of labeled data per class
                labelindperc=sort(randperm(numperc,labelperc)); % index of labeled data selected
                labelind=[]; % labelind: index of known label
                for i=1:c
                    labelind=[labelind labelindperc+(i-1)*numperc];
                end
                
                [result]=selfweightmklsemi_MV_Robust(Kernel,y,labelind,alpha,beta,gam)
                res=[res;result];
            end
            mn20=mean(res);
            res=sort(res,'descend');
            mn10=mean(res(1:10,:));
            plot=[plot;alpha gam mn20 mn10 beta];
        end
    end
end


save('ACC_SSC_Robust_Reuters.mat','plot');
