
alpha=1e-3;
beta=1000;
mu=1;
res=[];

for i=1:50
    r=0.5;%rate of labeled data
    
    [m,n,rr,tt]=size(Kernel);
    c=length(unique(y)); % number of class
    numperc=floor(n/c); % number of data per class
    labelperc=floor(r*numperc); % number of labeled data per class
    labelindperc=sort(randperm(numperc,labelperc)); % index of labeled data selected
    labelind=[]; % labelind: index of known label
    for i=1:c
        labelind=[labelind labelindperc+(i-1)*numperc];
    end
    
    [result]=selfweightmklsemi_MV_Robust(Kernel,y,labelind,alpha,beta,mu)
    % plot=[plot;alpha gama result];
    res=[res;result];
end




