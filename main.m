% PCA Chan-Vese for ultrasound images -- this is the main file
% Author:       Suvadip Mukherjee (sm5vp@virginia.edu)
% Co-Author:    Rituparna Sarkar  (rs9vj@virginia.edu)
% Affiliation:  VIVA Lab, Dept. of ECE, University of Virginia

% If you are using this code please cite the following paper
%
% "Sarkar, Rituparna, Suvadip Mukherjee, and Scott T. Acton. "Dictionary learning level set." 
% IEEE Signal Processing Letters 22.11 (2015): 2034-2038."
% -------------------------------------------------------------------------
% The KSVD toolbox is required to run this code. The toolbox is available
% at : http://www.cs.technion.ac.il/~elad/software/
clc,clear;
fname = './imCyl'; % directory containing the images
load(strcat(fname,'/gt_1.mat')); % ground truth for the images

% Add to path the KSVD toolbox folder
addpath('./ksvdbox');
addpath('./OMPbox');

imtype = 'png'; % change the image type 
resize_factor =  1;
img_stack       = readFolder(fname,imtype, resize_factor);
[Nr,Nc,N] = size(img_stack);

%% dictionay learning using k-svd_toolbox 
param_dict.InitializationMethod =  'DataElements';
param_dict.L = 3;   % number of elements in each linear combination.
param_dict.K = 8; % number of dictionary elements
param_dict.numIteration = 50; % number of iteration to execute the K-SVD algorithm.
param_dict.errorFlag = 0; % decompose signals until a certain error is reached. do not use fix number of coefficients.
%param.errorGoal = sigma;
display_dict        = 1;
param_dict.preserveDCAtom = 0;
do_LOOC = 1;
img_idx = 32;

[x,Dictionary,output]  = dict_learning(img_stack,display_dict,img_idx,do_LOOC,param_dict);

dict                = [ones(Nr*Nc,1),Dictionary];%/(Nr*Nc);  % First element is constant, the bias term
imshow(x)
%% Now use the DL-Chan vese

DL =1;
CV =0;
img_num =32;
if DL
    param.basis_vect        = dict(:,:);
elseif CV
    param.basis_vect        = dict(:,1);
end
param.Img               = img_stack(:,:,img_num);
param.n_regions         = 1;
param.mask_type         = 'multiball';
param.img_magnify       = 200;
if DL
    param.contour_color     = 'y';
elseif CV
    param.contour_color     = 'b';
end
param.display_intrvl    = 50;
param.num_iter          = 600;
param.convg_error       = 0.5;
param.length_term       = 0.4;   
param.lambda1           = 1;
param.lambda2           = 1; 
param.convg_count       = 10;
param.lambda_l2         = 0.3; 
param.evolution         = 0; 
param.evoldirectory     = './Evolution_results';
if DL
    param.fnameevolve   = strcat(param.evoldirectory ,'/',strtok(fname,'./'),'_DLevol',num2str(img_num));
elseif CV
    param.fnameevolve   = strcat(param.evoldirectory ,'/',strtok(fname,'./'),'_CVevol',num2str(img_num));
end

init_mask               = initRegions(param);
param.init_phi          = computeSDF(init_mask);
[phi,inside_image,outside_image,in_coef,out_coef,cn1,cn2] = ChanVeseDL(param);

%% calculate_DiCE coefficient 
final_segment = phi >= 0;
masked_final = param.Img.*final_segment;
% imshow(final_segment);
final_segment_rev = (final_segment);
grnd_truth = binary_stack(:,:,img_num);
image_intersec = and(final_segment_rev,grnd_truth);
dice_coeff = 2*sum(image_intersec(:))/(sum(final_segment_rev(:))+sum(grnd_truth(:)))
subplot(1,2,1);imshow(grnd_truth); 
subplot(1,2,2); imshow(final_segment_rev);

%%
figure(2); 
subplot(1,3,1); imshow(param.Img,'InitialMagnification',param.img_magnify);title('Initial');
subplot(1,3,2);imshow(inside_image,'InitialMagnification',param.img_magnify);title('p_{in}');
subplot(1,3,3);imshow(outside_image,'InitialMagnification',param.img_magnify);title('p_{out}');

