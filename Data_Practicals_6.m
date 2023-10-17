% Landsat images
clear, close all, clc
% Read separate files for each spectral band
% We use band 4 (Red, 640–670 nm), band 3
% Green, 530–590 nm) and band 2 (Blue, 450–510 nm), each of which has a 30
% meters resolution.
datafolder = 'C:\Users\alexamin\Dropbox (UiO)\practical-data\workshop-sandbox\mansour\LC08_L2SP_015009_20200830_20200906_02_T1\'

I1 = imread([datafolder,'LC08_L2SP_015009_20200830_20200906_02_T1_SR_B4.TIF']);
I2 = imread([datafolder,'LC08_L2SP_015009_20200830_20200906_02_T1_SR_B3.TIF']);
I3 = imread([datafolder,'LC08_L2SP_015009_20200830_20200906_02_T1_SR_B2.TIF']);

% Adaptive histogram equilization of the image
I1 = adapthisteq(I1,'ClipLimit',0.1,'Distribution','Rayleigh');
I2 = adapthisteq(I2,'ClipLimit',0.1,'Distribution','Rayleigh');
I3 = adapthisteq(I3,'ClipLimit',0.1,'Distribution','Rayleigh');
% 'ClipLimit' limits the contrast enhancement,
% 'Distribution' The histogram shape, such as Uniform, Rayleigh, or Exponential. 

% display image format info
whos

% The three bands are concatenated to a 24-bit RGB images
I = cat(3,I1,I2,I3);
figure, imshow(I)
%%
figure(1),clf
imshow(I,'InitialMagnification',10)
axis on
%%
I1 = I(6500:7000,3400:3900,:);
figure, imshow(I1), axis on
%%
I2 = imsharpen(I1);
% I2(:,:,1) = adapthisteq(I1(:,:,1),'ClipLimit',0.01,'Distribution','exponential');
% I2(:,:,2) = adapthisteq(I1(:,:,2),'ClipLimit',0.01,'Distribution','exponential');
% I2(:,:,3) = adapthisteq(I1(:,:,3),'ClipLimit',0.01,'Distribution','exponential');

%
figure
subplot(121),imshow(I1)
axis on
subplot(122),imshow(I2)
%%
I3 = rgb2gray(I2);

%%
I3 = adapthisteq(I3,'ClipLimit',0.01,'Distribution','exponential');
figure,imshow(I3)
axis on
%%
% I4 = imbinarize(I3, 'adaptive');
% figure(10),clf
% imshow(I4), hold on
% axis on
%
figure
I4 = edge(I3,'canny',0.3,3);
imshow(I4), hold on
axis on
%%
[H,theta,rho] = hough(I4,'theta',-89:0.1:89);
%
%H = medfilt2(medfilt2(H));
P  = houghpeaks(H,150,'Threshold',30);
figure, 
imagesc(theta,rho,H), colormap(hot)
xlabel('\theta'), ylabel('\rho');
hold on,
x = theta(P(:,2)); y = rho(P(:,1));
plot(x,y,'s','color','k');
%%
lines1 = houghlines(I4,theta,rho,P,'FillGap',3,'MinLength',7);
figure
subplot(121)
imshow(I3)
hold on,
for k = 1:length(lines1)
    xy = [lines1(k).point1; lines1(k).point2];
    plot(xy(:,1),xy(:,2),...
        'LineWidth',1,...
        'Color',[1 0 0]);
end
subplot(122)
imshow(I3)

%% testing Hough transform 
% Load and display the image
% I = imread('example_image.jpg');
% if size(I,3) == 3
%     I = rgb2gray(I);
% end
% figure,
% I =false(100);
% 
% I(50,20:80)=true;
% I(20:80,50)=true;
% imshow(I);
%% testing
% Define the edge detection filter (Sobel, for instance)
% sobel_x = [-1 0 1; -2 0 2; -1 0 1];
% sobel_y = sobel_x';
% 
% % Apply edge detection
% Gx = conv2(double(I), sobel_x, 'same');
% Gy = conv2(double(I), sobel_y, 'same');
% edges = sqrt(Gx.^2 + Gy.^2);
% edges = edges > 50; % threshold
% 
% % Hough transform parameters
% theta_resolution = 0.01;
% rho_resolution = 1;
% [r, c] = size(I);
% diagonal = sqrt(r^2 + c^2);
% rhos = -diagonal:rho_resolution:diagonal;
% thetas = -pi/2:theta_resolution:pi/2;
% 
% % Initialize Hough accumulator
% H = zeros(length(rhos), length(thetas));
% 
% % Fill the Hough accumulator
% for x = 1:c
%     for y = 1:r
%         if edges(y,x)
%             for theta_index = 1:length(thetas)
%                 theta = thetas(theta_index);
%                 rho = x * cos(theta) + y * sin(theta);
%                 rho_index = round(rho/diagonal * numel(rhos)/2 + numel(rhos)/2);
%                 H(rho_index, theta_index) = H(rho_index, theta_index) + 1;
%             end
%         end
%     end
% end
% 
% % Display the Hough accumulator
% figure;
% imshow(mat2gray(H), 'XData', thetas, 'YData', rhos);
% xlabel('\theta (radians)');
% ylabel('\rho');
% title('Hough Transform');
% axis on;
% colormap(gca, hot);


