function [result]=selfweightmkl_Robust(H,s,alpha,beta,gamma)
% s is the true class label.
warning off
[m,n,nn]=size(H);
S_Error=[];
iteration_No = [];
norm_S=[];
Acc=[];
sumK=zeros(n);
for i=1:nn
    sumK=sumK+H(:,:,i);
end
K=sumK/nn;
Z=eye(n);
for tt=1:n
  B(tt)= (K(tt,tt)-2*Z(tt,:)*K(tt,:)'+Z(tt,:)*K*Z(tt,:)'+1)^-0.5;
end

for ia=1:n
    E_all = zeros(n);
    E_all(ia,ia)= 1;
    E(:,:,ia)= E_all;
end
c=length(unique(s));
itt = 0;
j=0;
Accurc=[];
%options = optimset( 'Algorithm','interior-point-convex','Display','off');
for i=1:5
%     disp(j+1);
    itt = i;
    iteration=i
    Zold=Z;

    D = diag(sum(Z));
    L = D-Z;
    
    [F, temp, ev]=eig1(L, c, 0);
    sumH=zeros(n);
    sum1=0;
    for j=1:nn
        temp(j)=norm(H(:,:,j)-K,'fro');
        omega(j)=1/(2*temp(j));
        sumH=sumH+omega(j)*H(:,:,j);
        sum1=sum1+omega(j);
    end

    for ij=1:n
        for ji=1:n
            all(ji)=norm(F(ij,:)-F(ji,:));
        end
        ab=B(ij);

        Z(:,ij)=(4*K+gamma*eye(n))\(4*K(ij,:)'-alpha/4*all');
        % we use the free package to solve quadratic equation: http://sigpromu.org/quadprog/index.html
        %    [Z(:,ij),err,lm] = qpas(H,ff,[],[],ones(1,n),1,zeros(n,1),ones(n,1));
        % Z(:,ij)=quadprog(H,(beta/2*all'-(alpha-2)*K(:,ij))',[],[],ones(1,n),1,zeros(n,1),ones(n,1),Z(:,ij),options);
    end
    Z(find(Z<0))=0;
    Z= (Z+Z')/2;
    %      K=(Z'+beta*sumH)/(Z+beta*sum1*eye(n));
    
    sumB=zeros(n);
    for st=1:n
        Z_hat_all = zeros(n);
        Z_hat_all(st,:) = Z(st,:);
        Z_hat(:,:,st) = Z_hat_all;
        sumB = sumB + B(st)*(E(:,:,st)-2* Z_hat(:,:,st)+Z(st,:)'*Z(st,:));
    end
    
    K=(2*beta*sumH-sumB)/(2*beta*sum1);
    K(find(K<0))=0;
    K=(K'+K)/2;
    
    for tt=1:n
         B(tt)= (K(tt,tt)-2*Z(tt,:)*K(tt,:)'+Z(tt,:)*K*Z(tt,:)')^-0.5;
    end
    S_Error = [S_Error;(norm(Z-Zold)/norm(Zold))];
    iteration_No = [iteration_No;itt];
    %norm_S = [norm_S;(S-S_old)];
    if i>5 &((norm(Z-Zold)/norm(Zold))<1e-5)
        break
    end
% actual_ids= kmeans(F, c, 'emptyaction', 'singleton', 'replicates', 100, 'display', 'off');
% Acc = [Acc; ClusteringMeasure( actual_ids,s)];
end
save('S_Error.mat','S_Error');
% save('Iteration.mat','iteration_No');
%save('norm_S.mat','norm_S');
save('Acc_res.mat','Acc');

actual_ids= kmeans(F, c, 'emptyaction', 'singleton', 'replicates', 100, 'display', 'off');
[result] = [ClusteringMeasure( actual_ids,s) itt];