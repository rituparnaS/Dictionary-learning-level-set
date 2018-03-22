I = imread('28.png');
I2 = imresize(I,[60 100]);
imshow(I2)
imwrite(I2,'28_resz.png');

%%