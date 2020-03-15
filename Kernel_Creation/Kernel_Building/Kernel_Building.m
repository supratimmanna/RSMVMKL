Data=X;
BuildKernels(Data);
load('kernel_file.mat')
ii=ii+1;
Kernel(:,:,:,ii)=K;
% k4=Kernel;
% load('KK.mat')
% i=i+1;
% Kernel(:,:,i)=K;
%Kernel(:,:,:,6)=K;

% for i=1:16
%     Data=X(:,i);
%     BuildKernels(Data);
%     load('kernel_file.mat')
%     Kernel(:,:,:,i)=K;
% end


