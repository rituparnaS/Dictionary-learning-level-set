
%output argument will also have x when required to display the dictionary
function [x, dict_Basis, sprse_code] = dict_learning(stack,show,discard,do_LOOC,param_dict)

[Nr,Nc,N] = size(stack);
% Stack all signals in rows
Data = zeros(Nr*Nc,N);
for ii = 1 : N
    dummy = stack(:,:,ii);
    Data(:,ii) = dummy(:);
end

if do_LOOC
    Data(:,discard) = [];
end

data_mean = mean(Data,2);
[n1 N1]= size(Data);
%[full_principal_components,~,pc_var] = pca(Data); % If needed to test with
%PCA
for jj=1:1:N1
    data_meanzeroed(:,jj) = Data(:,jj)-data_mean;
end


params.dictsize = param_dict.K;
param.InitializationMethod =  param_dict.InitializationMethod;
param.L = param_dict.L;   % number of elements in each linear combination.
param.K = param_dict.K; % number of dictionary elements
param.numIteration = param_dict.numIteration; % number of iteration to execute the K-SVD algorithm.
param.errorFlag = 0; 
param.preserveDCAtom = 0;

[Dksvd,g]=KSVD(data_meanzeroed,param);

dict_Basis = Dksvd;
sprse_code = g;

% If needed to display dictionary
if show
    n  = 2;
    m  = ceil(params.dictsize/n);
    D  = Dksvd;
    sz = [Nr Nc];
    x  = showdict(D,sz,n,m,'whitelines')  ;
end

    
    
    