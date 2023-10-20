%%


websave('planet_img.jpg','https://www.solarsystemscope.com/textures/download/8k_moon.jpg')
I1 = imread('planet_img.jpg');
%%

figure,
imshow(I1,'XData',[-180 180],'YData',[90 -90]), axis on,
set(gca,'Ydir','normal'), xlabel('longitude (^{o})'), ylabel('latitude (^{o})')
title('Surface of the Moon')

%%
I2 = rgb2gray(I1);

figure
a = 1737;
f = 0.0012;
b = a*(1-f);
[x,y,z] = ellipsoid(0,0,0,a,a,b,400);
planet_surface = surface(x,y,z,'FaceColor','texture',...
            'EdgeColor','none','CData',flipud(I2),'DiffuseStrength',...
            1,'SpecularStrength',0,'FaceAlpha',1);
axis equal tight off, colormap(gray)    
    
% set axis clipping
ax = gca;
ax.Clipping = 'off';
% equal data unit lengths along each axis
axis equal;
% 3D view
view(3);
% Create a black image
%image = zeros(500, 500, 3, 'uint8');

% Draw a circle on the image
% center = [250, 250];
% radius = 50;
% image = insertShape(image, 'FilledCircle', [center, radius], 'Color', [0, 255, 0], 'Opacity', 1);

% Draw intersecting rectangle and line
% image = insertShape(image, 'Rectangle', [200, 200, 100, 100], 'Color', [255, 255, 0], 'LineWidth', 5);
% image = insertShape(image, 'Line', [200, 250, 300, 250], 'Color', [255, 255, 255], 'LineWidth', 5);
% Add Gaussian noise
%%
I2 = I2(1200:1456,500:756);

figure, 
imshow(I2)
%%
I3 = imsharpen(I2);
figure, 
imshow(I3)
%%
I4 = adapthisteq(I3,'ClipLimit',0.1,'Distribution','Rayleigh');
figure, 
imshow(I4)
%image_noisy = imnoise(image_gray, 'gaussian', 0, 0.02);

%%
I5 = medfilt2(I4);
figure, 
imshow(I5)
%%
% figure
% Show the noisy image with intersecting features
% imshow(image_noisy);
% title('Noisy Image with Intersecting Features');
% hold on
% Detect circles using Hough Transform

I6 = edge(I5,'canny',0.4, 3);
figure, 
imshow(I6)
%%
[centers, radii] = imfindcircles(I6, [5 30], 'Sensitivity', 0.85);

% Draw the detected circles
figure,
subplot(131),imshow(I6);
hold on, axis on
%viscircles(centers, radii, 'EdgeColor', 'b');
plot(centers(:,1), centers(:,2), '+b');
subplot(132),imshow(I5);
viscircles(centers, radii, 'EdgeColor', 'b');
hold on, axis on
subplot(133),imshow(I5);
hold on,
plot(centers(:,1), centers(:,2), '+b');
 axis on
 
%% crater statistics
% approximate conversion of pixels to km
dx = pi*a/size(I1,1);
D_km = 2*radii*dx; %crater diameter in km
% figure, 
% histogram(D_km,6), xlabel('Diameter (km)'), ylabel('Crater frequency (km)')
% title(['Number of craters ',num2str(length(radii))])
% hold on,
%
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
%
%plot(bin_centers, fitted_N, 'r', 'LineWidth', 2); % Fitted curve
logN = log(N);
logX = log(bin_centers);
k_fit = logN/[logX; logX*0+1];
fitted_N = exp(k_fit(2))*bin_centers.^(k_fit(1));
plot(bin_centers, fitted_N, 'r', 'LineWidth', 2)
hold off
xlabel('log(Crater Diameter) (km)');
ylabel('log(Crater Frequency)');
%legend('Observed', 'Fitted');
title(['A = ', num2str(k_fit(2)), ', k = ', num2str(k_fit(1))]);
set(gca, 'XScale', 'log', 'YScale', 'log');

%% testing
% function my_imfindcircles(image, radius_range)
%     % Convert image to grayscale if it's not
%     if size(image, 3) == 3
%         image = rgb2gray(image);
%     end
%     
%     % Blur the image
%     h = fspecial('gaussian', [5 5], 2);
%     image = imfilter(image, h);
%     
%     % Use edge detection
%     edges = edge(image, 'canny');
%     
%     [rows, cols] = size(image);
%     accumulator = zeros(rows, cols);
%     
%     % Iterate through the image
%     [y, x] = find(edges);
%     for i = 1:length(x)
%         for radius = radius_range(1):radius_range(2)
%             for theta = 0:360
%                 a = round(x(i) - radius * cosd(theta));
%                 b = round(y(i) - radius * sind(theta));
%                 
%                 if a > 0 && a <= rows && b > 0 && b <= cols
%                     accumulator(a, b) = accumulator(a, b) + 1;
%                 end
%             end
%         end
%     end
%     
%     % Find the circles
%     max_accumulator = max(accumulator(:));
%     [center_y, center_x] = find(accumulator >= max_accumulator * 0.9); % Threshold
%     
%     % Draw circles
%     figure, imshow(image, []), hold on;
%     for i = 1:length(center_x)
%         viscircles([center_x(i), center_y(i)], radius_range(1):1:radius_range(2), 'EdgeColor', 'b');
%     end
%     hold off;
% end


