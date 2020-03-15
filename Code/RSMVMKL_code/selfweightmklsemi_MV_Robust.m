function [result]=selfweightmklsemi_MV_Robust(Kernel,s,labelind,alpha,beta,gamma)
% s is the true class label.
warning off
S_Error=[];
res=[];
[m,n,nn,mm]=size(Kernel);
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
ulabel= setdiff(1:n,labelind);%index of unlabeled data

ll=length(labelind);
Fl=zeros(ll,c);
for i=1:ll
    Fl(i,s(labelind(i)))=1;
end
Fu=zeros(n-ll,c);
for i=1:50
    iter=i
    Sold=S;
    
    for p = 1:mm
        for q = 1:nn
            inst(p,q)= norm((Kernel(:,:,q,p)-K),'fro');
            s1(p,q)= 1/(2*inst(p,q));
        end
    end
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
    
    sumB=zeros(n);
    for st=1:n
        S_hat_all = zeros(n);
        S_hat_all(:,st) = S(:,st);
        S_hat(:,:,st) = S_hat_all';
        sumB = sumB + B(st)*(EE(:,:,st)-2* S_hat(:,:,st)+S(:,st)'*S(:,st)*eye(n));
    end
    
    K=(2*beta*sumZH-sumB)/(2*beta*sumZ);
    K(find(K<0))=0;
    K=(K'+K)/2;
    F=[Fl;Fu];
    parfor ij=1:n
        %         for ji=1:n
        %             all(ji)=norm(F(ij,:)-F(ji,:));
        %         end
        [all]=veccomp2(ij,n,F);
        ab=B(ij);
        S(:,ij)=(ab*K+gamma*eye(n))\(ab*K(ij,:)'-alpha/4*all');
    end
    S(find(S<0))=0;
    S= (S+S')/2;
    
    for tt=1:n
         b(tt)= (K(tt,tt)-2*K(tt,:)*S(:,tt)+S(:,tt)'*K*S(:,tt))^0.5;
         B(tt) = 1/(b(tt)+0.001);
%          B=abs(B);
    end
    
    D = diag(sum(S));
%     L = D-S;
    L=eye(n)-sqrt(inv(D))*S*sqrt(inv(D));
    L=(L+L')/2;
    
    uu=zeros(n-ll,n-ll);
    for ii=1:(n-ll)
        for jj=1:(n-ll)
            uu(ii,jj)=L(ulabel(ii),ulabel(jj));
        end
    end
    ul=zeros(n-ll,ll);
    for ii=1:(n-ll)
        for jj=1:ll
            ul(ii,jj)=L(ulabel(ii),labelind(jj));
        end
    end
    Fu=-uu\(ul*Fl);
    
    S_Error = [S_Error;(norm(S-Sold)/norm(Sold))];
    
    if i>5 &((norm(S-Sold)/norm(Sold))<1e-3)
        break
    end
    
end
save('S_Error_semi_Robust.mat','S_Error');

[ur,uc]=size(ulabel);
[max_value,max_ind] = max(Fu,[],2);
cnt = 0;
for i = 1:uc
    if max_ind(i) == s(ulabel(i))
        cnt = cnt+1;
    end
end
result = cnt/uc;
res=[result iter];
result=res;
end

function [all]=veccomp2(ij,n,F)
warning off
for ji=1:n
    all(ji)=norm(F(ij,:)-F(ji,:));
end
end