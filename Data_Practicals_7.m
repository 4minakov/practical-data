%% clear workspace
clear all
close all
%% dowload and read image
websave('moon.jpg','https://www.solarsystemscope.com/textures/download/8k_moon.jpg')
I1 = imread('moon.jpg');
%% convert to grayscale intensity 
I2 = rgb2gray(I1);
%% plot image
figure,
imshow(I2,'XData',[-180 180],'YData',[90 -90]), axis on,
set(gca,'Ydir','normal'), xlabel('longitude (^{o})'), ylabel('latitude (^{o})')
title('Surface of the Moon'),
%% planet parameters
a = 1738; % equatorial radius of the Moon (km)
f = 0.0012; % flattening
b = a*(1-f); % polar axis of the Moon
%% 3D view
figure
[x,y,z] = ellipsoid(0,0,0,a,a,b,400); % Moon's ellipsoid 
surface(x,y,z,'FaceColor','texture','EdgeColor','none',...
    'CData',flipud(I1),'DiffuseStrength',1,'SpecularStrength',0,'FaceAlpha',1);
axis equal tight off, colormap(gray)    
ax = gca; ax.Clipping = 'off'; view(3);
%% select an area
I3 = I2(1200:1456,500:756);
figure, imshow(I3), title('Selected area')
%% Sharpening of the image (high-pass gaussian filter)
I4 = imsharpen(I3);
figure, imshow(I4), title('Sharpened image')
%% adaptive histogram equilization (CLAHE)
I4 = adapthisteq(I3,'ClipLimit',0.1,'Distribution','Rayleigh');
figure, imshow(I4), title('Enhanced contrast')
%% median filter to remove speck noise
I5 = medfilt2(I4);
figure, imshow(I5), title('Denoised image')
%% Find edges in the image
% The Canny method finds edges by looking for local maxima of the gradient of image. 
% Parameter 1: upper threshold for detection
% Parameter 2: the standard deviation of the Gaussian filter to compute
% gradient magnitude
I6 = edge(I5,'canny',0.4, 3);
figure, imshow(I6), title('Edges detected')
%% circle detection using Hough transform
[centers, radii] = imfindcircles(I6, [5 30], 'Sensitivity', 0.85);
% The routine finds circles with radii in the search range [MIN_RADIUS MAX_RADIUS]
% Sensitivity factor for finding circles in the range [0 1].
% A high sensitivity detecting more circles at the risk of a higher false detection rate.
%% Draw the detected circles
figure,
subplot(131),imshow(I6);
hold on, axis on;
plot(centers(:,1), centers(:,2), '+b');
subplot(132),imshow(I5);
viscircles(centers, radii, 'EdgeColor', 'b');
hold on, axis on
subplot(133),imshow(I5);
hold on, plot(centers(:,1), centers(:,2), '+b'), axis on
%% crater statistics
% approximate conversion of pixels to km
dx = pi*a/size(I1,1);
D_km = 2*radii*dx; %crater diameter in km
% Fitting power-law distribution
% Create a histogram with specified bin edges
bin_edges = 15:6:60; % Adjust bin edges as needed
[N, edges] = histcounts(D_km, bin_edges);
% Calculate bin centers
bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
% Plot the results
figure(20),clf
hold on
plot(bin_centers, N,'s'), 
bar(bin_centers, N, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k'); % Observed data
logN = log(N);
logX = log(bin_centers);
k_fit = logN/[logX; logX*0+1];
fitted_N = exp(k_fit(2))*bin_centers.^(k_fit(1));
plot(bin_centers, fitted_N, 'r', 'LineWidth', 2)
hold off
xlabel('log(Crater Diameter) (km)');
ylabel('log(Crater Frequency)');
title(['A = ', num2str(k_fit(2)), ', k = ', num2str(k_fit(1))]);
set(gca, 'XScale', 'log', 'YScale', 'log');

%%  testing 
% % Create a black image
% image = zeros(500, 500, 3, 'uint8');
% % Draw circle 1 on the image
% center = [250, 250];
% radius = 50;
% image = insertShape(image, 'FilledCircle', [center, radius], 'Color', [0, 255, 0], 'Opacity', 1);
% % Draw  circle 2 on the image
% center = [150, 150];
% radius = 20;
% image = insertShape(image, 'FilledCircle', [center, radius], 'Color', [0, 255, 0], 'Opacity', 1);
% % Draw  circle 3 on the image
% center = [170, 130];
% radius = 30;
% image = insertShape(image, 'FilledCircle', [center, radius], 'Color', [0, 255, 0], 'Opacity', 1);
% % Draw intersecting rectangle and line
% %image = insertShape(image, 'Rectangle', [200, 200, 100, 100], 'Color', [255, 255, 0], 'LineWidth', 5);
% %image = insertShape(image, 'Line', [200, 250, 300, 250], 'Color', [255, 255, 255], 'LineWidth', 5);
% % Add Gaussian noise
% image_gray = rgb2gray(image);
% image_noisy = imnoise(image_gray, 'gaussian', 0, 0.02);
% image_edges = edge(image_noisy,'canny',0.3, 2);
% figure,imshow(image_edges)
% % Show the noisy image with intersecting features
% figure,imshow(image_noisy);
% title('Noisy Image with Intersecting Features');
% % Detect circles using Hough Transform
% [centers, radii] = imfindcircles(image_edges, [10 100], 'Sensitivity', 0.85);
% % Draw the detected circles
% figure,imshow(image_noisy),
% title('Detected Circles');
% viscircles(centers, radii, 'EdgeColor', 'b');


