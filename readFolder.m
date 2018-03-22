function [ imstack ] = readFolder(fname,imtype,f)
%READFOLDER read all the images in the folder 'fname' and type 'imtype'.
%Resize by the factor 'f';

file_type       = strcat('*.',imtype);
filenames       = dir(fullfile(fname, file_type));  % read all images with specified filetype
total_files     = numel(filenames);

% Extract a single file to get dimensions

file_name = fullfile(fname,filenames(1).name);
needs_conversion = 0;
img = imread(file_name);
if length(size(img)) > 2
   needs_conversion = 1; 
   img = rgb2gray(img);
   img = imresize(img,f);
   [Nr,Nc,junk] = size(img);
   
else
   img = imresize(img,f); 
   [Nr,Nc] = size(img);
end

imstack = zeros(Nr,Nc,total_files);

for ii = 1 : total_files
    
    file_name = fullfile(fname,filenames(ii).name);
    img = imread(file_name);
    if needs_conversion
        img = rgb2gray(img);
    end
    img = imresize(img,f);
    img = im2double(img);
    imstack(:,:,ii) = img;
    
end

