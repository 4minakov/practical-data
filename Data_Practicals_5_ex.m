clear all
close all
%% download and save image 
websave('ripples.jpg','https://www.uahirise.org/images/2009/details/cut/ESP_011765_1780_cut.jpg')
I = imread('ripples.jpg');
%%
% convert to grayscale
I1 = double(rgb2gray(I));
% subtract mean
meanI = mean(I1(:));
I1 = I1 - meanI;
% make coordinate vectors
Lx = 400; % width of the image is 400 m
[ny,nx] = size(I1); %number of pixels in x and y direction
dx = Lx/nx;% size of the pixel
x = dx/2:dx:dx*nx-dx/2; % x-vector
y = dx/2:dx:dx*ny-dx/2; % y-vector
[X,Y]=meshgrid(x,y);
%% Make wavenumbers
%Nyquist wavenumbers
knx = 1/2/dx;
kny = 1/2/dx;
%wavenumber intervals
dkx = 1/dx/nx; 
dky = 1/dx/ny; 
%wavenumber vectors
kx = -knx:dkx:knx-dkx; 
ky = -kny:dky:kny-dky; 
[KX,KY]=meshgrid(kx,ky);
K = sqrt(KX.^2+KY.^2);
%% Fourier transformed image (complex numbers)
Sv = fftshift(fft2(I1));
%% Isotropic Gaussian filtering
ax = 0.1; %smoothing along X-direction
ay = ax; %smoothing along Y-direction
FLT = exp(-(KX.^2/ax^2 + KY.^2/ay^2)); %gaussian filter
Sv_flt = Sv.*FLT; %convolution
I1_flt = real(ifft2(ifftshift(Sv_flt))); %convert to spatial domain via inverse FFT
%% Residual image and spectrum
I1_res = I1 - I1_flt;
Sv_res = Sv-Sv_flt;
%% Plotting
figure(1),
subplot(231),
imagesc(x,y,I1), xlabel('x [m]'), ylabel('y [m]'), axis equal tight
colormap(gray), title('Meridiani Planum')
subplot(234),
imagesc(kx,ky,abs(Sv)), shading flat, colormap(gray(10)),
axis equal tight, set(gca,'Ydir','normal')
xlabel('k_x [m^{-1}]'), ylabel('k_y [m^{-1}]'), caxis([0 1]*1e6)
title('Spectrum |S|')
xlim([-1 1]/2), ylim([-1 1]/2)
subplot(232)
imagesc(x,y,I1_flt), shading flat, colormap(gray(10)),
axis equal tight,  
xlabel('x [m]'), ylabel('y [m]'), axis equal tight
title('Filtered Image')
subplot(235)
imagesc(kx,ky,abs(Sv_flt)), shading flat, colormap(gray(10)),
axis equal tight, set(gca,'Ydir','normal'), caxis([0 1]*1e6)
xlabel('k_x [m^{-1}]'), ylabel('k_y [m^{-1}]'), 
title('Filtered spectrum')
xlim([-1 1]/2), ylim([-1 1]/2)
subplot(233)
imagesc(x,y,I1_res), shading flat, colormap(gray(10)),
axis equal tight,  caxis([-1 1]*2*std(I1_res(:)))
xlabel('x [m]'), ylabel('y [m]'), axis equal tight
title('Residual Image ')
subplot(236)
imagesc(kx,ky,abs(Sv_res)), shading flat, colormap(gray(10)),
axis equal tight, set(gca,'Ydir','normal'), caxis([0 1]*1e6)
xlabel('k_x [m^{-1}]'), ylabel('k_y [m^{-1}]'), 
title('Residual spectrum')
xlim([-1 1]/2), ylim([-1 1]/2)
%% Directional Gaussian filtering
theta = -5 * pi/180; %rotation angle
% Make 2D rotation matrix
R = [cos(theta) -sin(theta);
     sin(theta)  cos(theta)];
% Apply rotation to wavenumbers
KX1 = R(1,1)*KX + R(1,2)*KY;
KY1 = R(2,1)*KX + R(2,2)*KY;
%% Directional Gaussian filtering
ax = 2; %Smoothing along x-direction
ay = 0.1;%Smoothing along y-direction
FLT1 = exp(-(KX1.^2/ax^2 + KY1.^2/ay^2)); %gaussian filter
Sv_flt1 = Sv.*FLT1.*(1-FLT);%Convolution of image with combination of anisotropic 
% low-pass filter and isotropic low-cut filter
I1_flt1 = real(ifft2(ifftshift(Sv_flt1))); %convert to spatial domain via inverse FFT
%% Residual image and spectrum
I1_res1 = I1 - I1_flt1;
Sv_res1 = Sv-Sv_flt1;
%% Plotting
figure(2),
subplot(231),
imagesc(x,y,I1), xlabel('x [m]'), ylabel('y [m]'), axis equal tight
colormap(gray), title('Meridiani Planum')
subplot(234),
imagesc(kx,ky,abs(Sv)), shading flat, colormap(gray(10)),
axis equal tight, set(gca,'Ydir','normal')
xlabel('k_x [m^{-1}]'), ylabel('k_y [m^{-1}]'), caxis([0 1]*1e6)
title('Spectrum |S|')
xlim([-1 1]/2), ylim([-1 1]/2)
subplot(232)
imagesc(x,y,I1_flt1), shading flat, colormap(gray(10)),
axis equal tight,  caxis([-1 1]*2*std(I1_flt1(:)))
xlabel('x [m]'), ylabel('y [m]'), axis equal tight
title('Filtered Image ')
subplot(235)
imagesc(kx,ky,abs(Sv_flt1)), shading flat, colormap(gray(10)),
axis equal tight, set(gca,'Ydir','normal'), caxis([0 1]*1e6)
xlabel('k_x [m^{-1}]'), ylabel('k_y [m^{-1}]'), 
title('Filtered spectrum')
xlim([-1 1]/2), ylim([-1 1]/2)
subplot(233)
imagesc(x,y,I1_res1), shading flat, colormap(gray(10)),
axis equal tight,  caxis([-1 1]*2*std(I1_res(:)))
xlabel('x [m]'), ylabel('y [m]'), axis equal tight
title('Residual Image ')
subplot(236)
imagesc(kx,ky,abs(Sv_res1)), shading flat, colormap(gray(10)),
axis equal tight, set(gca,'Ydir','normal'), caxis([0 1]*1e6)
xlabel('k_x [m^{-1}]'), ylabel('k_y [m^{-1}]'), 
title('Residual spectrum')
xlim([-1 1]/2), ylim([-1 1]/2)
%% Find dominant wavenumber and wavelength
% smoothing using 3x3 cells running average of the 2D amplitude spectrum
Sv_med = medfilt2(medfilt2(medfilt2(abs(Sv_flt1)))); 
% Find index corresponding to the maximum value of the spectrum
[~,ii] = max(Sv_med(:));
%% Plotting
figure(3),clf
pcolor(KX1,KY1,Sv_med), shading flat, colormap(gray(10)),
hold on, plot(KX1(ii),KY1(ii),'^','MarkerFaceColor','r')
text(KX1(ii),KY1(ii),['\lambda_x = ',num2str(abs(round(1./KX1(ii),1))),' m'],...
    'color','r','FontSize',14)
axis equal tight, set(gca,'Ydir','normal'), caxis([0 1]*1e6)
xlabel('k_x [m^{-1}]'), ylabel('k_y [m^{-1}]'), 
title('Filtered spectrum (Running av)')
xlim([-1 1]/2), ylim([-1 1]/2)
%% Testing on a simple model
% nx=2^7;
% L=3*pi;
% [x1,y1] = meshgrid(linspace(-L/2,L/2,nx),linspace(-L/2,L/2,nx));
% 
% theta = pi/4;
% 
% zz = sin(x1)+cos(3*y1);
% zz = zz - mean(zz(:));
% 
% xr = x1*cos(theta)-y1*sin(theta);
% yr = x1*sin(theta)+y1*cos(theta);
% 
% zzr = sin(xr)+cos(3*yr);
% zzr = zzr - mean(zzr(:));
% 
% dx = x1(1,2)-x1(1,1);
% %Nyquist wavenumbers
% knx = 1/2/dx;
% kny = 1/2/dx;
% %wavenumber intervals
% dkx = 1/dx/nx; 
% dky = 1/dx/nx; 
% %wavenumber vectors
% kx = -knx:dkx:knx-dkx; 
% ky = -kny:dky:kny-dky; 
% [kx,ky]=meshgrid(kx,ky);
% Fzz = fftshift(fft2(zz));
% Fzzr = fftshift(fft2(zzr));
% 
% figure(4), 
% subplot(221),contourf(x1,y1,zz), axis equal tight
% subplot(222),contourf(kx,ky,abs(Fzz)), axis equal tight,
% xlim([-1 1]/L*10),ylim([-1 1]/L*10),
% subplot(223),contourf(x1,y1,zzr), axis equal tight
% subplot(224),contourf(kx,ky,abs(Fzzr)), axis equal tight,
% xlim([-1 1]/L*10),ylim([-1 1]/L*10),




