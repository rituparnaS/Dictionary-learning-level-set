function [D] = PCADictionary(stack,energy_thresh,show,discard,do_LOOC)
%PCADICTIONARY Compute the PCA of the images and return the major principal
%vectors as column vectors of 'D'. By default, the first column of D is the mean signal 

[Nr,Nc,N] = size(stack);
% Stack all signals in rows
Data = zeros(N,Nr*Nc);
for ii = 1 : N
    dummy = stack(:,:,ii);
    Data(ii,:) = dummy(:)';
end

if do_LOOC
    Data(discard,:) = [];
end

data_mean = mean(Data);
[full_principal_components,~,pc_var] = pca(Data);

pc_var = pc_var/(sum(pc_var(:)));
t = 0;
flag = 0;
for ii = 1 : length(pc_var)
    t = t+pc_var(ii);
    if t <= energy_thresh
        flag = flag+1;
    else
       break; 
    end
end

% to_add = repmat(data_mean',1,flag);
D = [data_mean' full_principal_components(:,1:flag)];

if show
%    flag = 50; 
   k = ceil(sqrt(flag)); 
   figure(1); 
   for ii = 1 : flag
      data = D(:,ii+1);
      img = reshape(data,Nr,Nc);
      subplot(k,k,ii);imshow(img,[]);
   end
    
end

end

