
IMAGE_PATH="Lenna.png";

%% 1) Import image and make sure it's gray scale and square
img_original = imread(IMAGE_PATH);
img_cropped = imcrop(img_original, [0 0 min(size(img_original(:,:,1))) min(size(img_original(:,:,1))) ]);
img_gs = rgb2gray(img_cropped);

img = double(img_gs);

imwrite(uint8(img), "images/grayscale.png");

%% 2) use svd to extract singular values of the image

[u,s,v] = svd(img);

%% Calculate the summation of the singular values

sv_sum = trace(s)

%% Sum of the largest 10 elements vs total sum
% Since 'svd' returns the eigenvalues sorted in decending order we don't need to sort

sum_of_10 = sum(diag(s(1:10,1:10)))

ratio_for_10sv = sum_of_10 / sv_sum

%% Image reconstruction using 10 largest singular values

img_10sv = u(:,1:10)*s(1:10,1:10)*v(:,1:10)';
imwrite(uint8(img_10sv), "images/img_reconstructed_10sv.png");

%% Repeat with 50 singular values
sum_of_50 = sum(diag(s(1:50,1:50)))

ratio_for_50sv = sum_of_50 / sv_sum

img_50sv = u(:,1:50)*s(1:50,1:50)*v(:,1:50)';
imwrite(uint8(img_50sv), "images/img_reconstructed_50sv.png");
