function [result]=selfweightmkl_MV_Robust(Kernel,s,alpha,beta,gamma)
% s is the true class label.
warning off
[m,n,nn,mm]=size(Kernel);
S_Error=[];
iteration_No = [];
% norm_S=[];
Acc_MV=[];

%%%%% Initialize the Kernel matrix(Kv) for each view
x=zeros(n);
for i=1:mm
    for j=1:nn
        x=x+Kernel(:,:,j,i);
    end
    sumKv(:,:,i)=x;
    x=zeros(n);
end
Kv=sumKv/nn;

%%%%% Initialize the optimal Kernel matrix

sumK=zeros(n);
for p=1:mm  %%%%%%%%%%%%%%%%%%%% changed to nn instead of mm
    sumK=sumK+Kv(:,:,p)
end
K=sumK/mm;

%%%% Initialization of adjacency or similarity matrix
S=eye(n);

% for tt=1:n
%   b(tt)= (K(tt,tt)-2*K(tt,:)*S(:,tt)+S(:,tt)'*K*S(:,tt))^0.5;
%   B(tt) = 1/b(tt);
% end
for tt=1:n
    B(tt)=1/tt;
end
% B=ones(1,n);
for ia=1:n
    E_all = zeros(n);
    E_all(ia,ia)= 1;
    EE(:,:,ia)= E_all;
end
%%%% Calculate the total number of classes
c=length(unique(s));

itt = 0;
%options = optimset( 'Algorithm','interior-point-convex','Display','off');
for i=1:100
    itt = i;
    S_old=S;
%     S=abs(S);
    iteration=i
    
 %%%% Creation of degree matrix D and laplacian matrix L   
    D = diag(sum(S));
    L=eye(n)-sqrt(inv(D))*S*sqrt(inv(D));
    L=(L+L')/2;
%     L = D-S;
 %%%% Eigen value calculation of the Laplacian matrix   
    [F, temp, ev]=eig1(L, c, 0);
%     F=abs(F);
%     opts.disp = 0;
%     [V E] = eigs(L,c,'LA',opts);  
%     F = V(:,1:ceil(c*1));
%     norm_mat = repmat(sqrt(sum(F.*F,2)),1,size(F,2));
%     %%avoid divide by zero
%     for i=1:size(norm_mat,1)
%         if (norm_mat(i,1)==0)
%             norm_mat(i,:) = 1;
%         end
%     end
%     F = F./norm_mat;
    
 %%%% Weight calculation for the optinal kernel of each view    
%     sumHv=zeros(n);
%     sumQ=0;
%     for j=1:mm
%         temp(j)=norm(Kv(:,:,j)-K,'fro');
%         Q(j)=1/(2*temp(j));
%         sumHv=sumHv+Q(j)*Kv(:,:,j);
%         sumQ=sumQ+Q(j);
%     end
    
 %%%% Weight calculation for each kernel of the multiple kernels
%  sumH=zeros(n);
%  sumZ=0;
 for p = 1:mm
     for q = 1:nn
         inst(p,q)= norm((Kernel(:,:,q,p)-K),'fro');
         s1(p,q)= 1/(2*inst(p,q));
     end
 end
 
%  for v=1:mm
%     sumZ(v)=sum(s1(v,:)) %%%% summation of Z(v,i)
%  end 
sumZ=sum(s1(:));

xx=zeros(n);
for tt=1:mm
    for ss=1:nn
        xx=xx+s1(tt,ss)*Kernel(:,:,ss,tt);
    end
    sumZH_v(:,:,tt)=xx; %%%% summation of Z(v,i)*Kernel(:,:,v,i)
    xx=zeros(n);
end

sumZH=zeros(n);
for st=1:mm
    sumZH=sumZH+sumZH_v(:,:,st);
end

%%%%% Calculation of Tr(F'LF) and Similarity matrix S

    for ij=1:n
        for ji=1:n
            all(ji)=norm(F(ij,:)-F(ji,:));
        end
        ab=B(ij);

        S(:,ij)=(ab*K+gamma*eye(n))\(ab*K(ij,:)'-alpha/4*all');
%         S(:,ij)=(K+gamma*eye(n))\(K(ij,:)'-alpha/4*all');
    end
    
    S(find(S<0))=0;
    S= (S+S')/2;
    S=abs(S);
  
 %%%%% Calculation of optimal kernel matrix Kv for vth view
%  for vv=1:mm
%      Kv(:,:,vv)=(beta*sumZH(:,:,vv)+zeta*Q(vv)*K)/(beta*sumZ(vv)+zeta*Q(vv))
%      Kv(find(Kv<0))=0;
%      Kv_transpose=permute(Kv,[2 1 3]);
%      Kv=(Kv+Kv_transpose)/2;
%  end
     
 sumB=zeros(n);
    for st=1:n
        S_hat_all = zeros(n);
        S_hat_all(:,st) = S(:,st);
        S_hat(:,:,st) = S_hat_all';
        sumB = sumB + B(st)*(EE(:,:,st)-2* S_hat(:,:,st)+S(:,st)'*S(:,st)*eye(n));
    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%% Calculation of optimal Kernel Matrix
    %      K=(Z'+beta*sumH)/(Z+beta*sum1*eye(n));
    K=(2*beta*sumZH-sumB)/(2*beta*sumZ);
    K(find(K<0))=0;
    K=(K'+K)/2;
    K=abs(K);
    
    for tt=1:n
         b(tt)= (K(tt,tt)-2*K(tt,:)*S(:,tt)+S(:,tt)'*K*S(:,tt))^0.5;
         B(tt) = 1/(b(tt)+0.001);
%          B=abs(B);
    end
    
    S_Error = [S_Error;(norm(S-S_old)/norm(S_old))];
     iteration_No = [iteration_No;itt];
%     norm_S = [norm_S;norm(S-S_old)];
    if i>5 &((norm(S-S_old)/norm(S_old))<1e-04)
        break
    end
actual_ids= kmeans(F, c, 'emptyaction', 'singleton', 'replicates', 100, 'display', 'off');
Acc_MV = [Acc_MV; ClusteringMeasure( actual_ids,s)];
    
end
save('S_Error_robust.mat','S_Error');
save('Iteration_3rd.mat','iteration_No');
% save('norm_S.mat','norm_S');
save('Acc__MV_res_3rd.mat','Acc_MV');

actual_ids= kmeans(F, c, 'emptyaction', 'singleton', 'replicates', 100, 'display', 'off');
[result] = [ClusteringMeasure( actual_ids,s) itt];