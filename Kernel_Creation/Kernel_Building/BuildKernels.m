function [f] = BuildKernels(dataset)
% global data_dir
% global kernel_dir
global KernelTypes KernelPostProcessTypes PolynomialDegrees PolyPlusDegrees GaussianDegrees 
% global kernel_path_part_1
%create a folder named by the name of dataset

% data_dir = fullfile(pwd, '..', 'data');
% 
% % kernel_dir = fullfile(data_dir, [dataset, '_kernel']);
% if exist(kernel_dir, 'dir') == 0
%     mkdir(kernel_dir);
% end
% 
% load(fullfile(data_dir, dataset));
% if exist('X', 'var')
%     % error([dataset, ' do not have variable X!']);
%     is_mv = 0;
% elseif exist('views', 'var')
%         % error([dataset, ' do not have variable views!']);
%         is_mv = 1;
%         for iView = 1:length(views)
%             if eval(sprintf('exist(''%s'', ''var'') == 0', ['X', num2str(iView)]))
%                 eval(sprintf('error([dataset, '' do not have variable %s!'']);', ['X', num2str(iView)]));
%             end
%         end
% else 
%     error([dataset, ' do not have variable X or views!']);
% end

KernelTypes = {'Gaussian', 'Linear', 'Polynomial', 'PolyPlus'};
KernelPostProcessTypes = {'Sample-Scale'};
PolynomialDegrees =  [2, 4];
PolyPlusDegrees = [2, 4];
GaussianDegrees = [0.01, 0.05, 0.1, 1, 10, 50, 100];


% if is_mv
%     for iView = 1:length(views)
%         kernel_path_part_1 = [dataset, '_kernel_view', num2str(iView) '_'];
%         clear X;
%         eval(sprintf('X = %s;', ['X', num2str(iView)]));
%         BuildSingleKernels(X);
%     end
% else
%     kernel_path_part_1 = [dataset, '_kernel_'];
%     BuildSingleKernels(X);
% end
% 

%X = Data;
BuildSingleKernels(dataset); %%%% ADded Later
 f = 1;

 end



function f2 = BuildSingleKernels(X)
%global kernel_dir
global KernelTypes KernelPostProcessTypes PolynomialDegrees PolyPlusDegrees GaussianDegrees 
% global kernel_path_part_1
j=0;

for kernel_type = KernelTypes
    kernel_option = [];
    switch lower(kernel_type{1})
        
        case lower('Gaussian')
            kernel_option.KernelType = 'Gaussian';
            for iKernelParam = GaussianDegrees
                kernel_option.t = iKernelParam;
                %kernel_path_part_2 = ['gaussian_', num2str(iKernelParam), '_'];
                for iPost = KernelPostProcessTypes
                    j=j+1;
                    D = EuDist2(X, [], 0);
                    max_D = max(D(:));
                    max_D = sqrt(max_D);
                    kernel_option.t = kernel_option.t * max_D;
                    K0 = constructKernel(X, [], kernel_option);
                    % K0 = exp(- D / (2 * kernel_option.t^2) );
                    K(:,:,j) = KernelNormalize(K0, iPost{1}); %ok
                    
                end
            end
            save('kernel_file_F.mat', 'K');
            
        case lower('Linear')
            kernel_option.KernelType = 'Linear';
            for iPost = KernelPostProcessTypes
                j=j+1;
                %kernel_path_part_2 = 'linear_';
                K0 = constructKernel(X, [], kernel_option);
                K(:,:,j) = KernelNormalize(K0, iPost{1}); %ok
                %save('kernel_file.mat', 'K');
            end
           
                
            
        case lower('Polynomial')
            kernel_option.KernelType = 'Polynomial';
            for iKernelParam = PolynomialDegrees
                kernel_option.d = iKernelParam;
                for iPost = KernelPostProcessTypes
                    j=j+1;
                    %kernel_path_part_2 = ['polynomial_', num2str(iKernelParam), '_'];
                    K0 = constructKernel(X, [], kernel_option);
                    K(:,:,j) = KernelNormalize(K0, iPost{1}); %ok
                    save('kernel_file_PLNML.mat', 'K');
                end
                               
            end
        case lower('PolyPlus')
            kernel_option.KernelType = 'PolyPlus';
            for iKernelParam = PolyPlusDegrees
                kernel_option.d = iKernelParam;
                for iPost = KernelPostProcessTypes
                    j=j+1;
                    %kernel_path_part_2 = ['polyplus_', num2str(iKernelParam), '_'];
                    K0 = constructKernel(X, [], kernel_option);
                    K(:,:,j) = KernelNormalize(K0, iPost{1}); %ok
                     
                end
            end
             save('kernel_file.mat', 'K');
        
        
    end
end
f2 = 1;
end