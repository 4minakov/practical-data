%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
close all
% Parameters

Ix = 300;
Iy = 100;
Iz = 100;

nx = 501;               % Must be odd
ny = 251;
nz = 1;%251;

nx_inn = 301;           % Must be odd
ny_inn = 151;
nz_inn = 1;%151;

dx = 10;
dy = 10;
dz = 10;

sigma = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grid

% Extendd domain

x = 0:dx:dx*(nx-1);
y = 0:dy:dy*(ny-1);
z = 0:dz:dz*(nz-1);

[X,Y,Z] = ndgrid(x,y,z);

% Inner subdomain

dx_inn = 0.5*(nx-nx_inn);
dy_inn = 0.5*(ny-ny_inn);
dz_inn = 0.5*(nz-nz_inn);

v_inn=zeros(nx_inn,ny_inn);  

x_inn = 0:dx:dx*(nx_inn-1);
y_inn = 0:dy:dy*(ny_inn-1);
z_inn = 0:dz:dz*(nz_inn-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Random field generation

% m=rand(nx,ny,nz) - 0.5;%array of random numbers
m=randn(nx,ny,nz); %array of standard normal random numbers

fm=fftn(m);

cx=(nx-1)*dx/2;
cy=(ny-1)*dy/2;
cz=(nz-1)*dz/2;
    
for ix=1:nx
    for iy=1:ny
        for iz=1:nz
             
             Nr(ix,iy,iz) = exp(-sqrt((x(ix)-cx)^2/Ix^2 + (y(iy)-cy)^2/Iy^2 + (z(iz)-cz)^2/Iz^2));

        end
    end
end


% Define the angle of rotation in radians (e.g., 45 degrees)
theta = pi/4; % 45 degrees

% Create the 2D rotation matrix
R = [cos(theta) -sin(theta);
     sin(theta)  cos(theta)];

% Display the rotation matrix
disp('Rotation Matrix:');
disp(R);

for ix=1:nx
    for iy=1:ny
        for iz=1:nz
             
            x1 = R(1,1)*(x(ix)-cx) + R(1,2)*(y(iy)-cy);
            y1 = R(2,1)*(x(ix)-cx) + R(2,2)*(y(iy)-cy);
            
            %Nr(ix,iy,iz) = exp(-sqrt(x1^2/Ix^2 + y1^2/Iy^2 + z(iz)^2/Iz^2)); %exponential
            Nr(ix,iy,iz) = exp(-(x1^2/Ix^2 + y1^2/Iy^2 + z(iz)^2/Iz^2)); %gaussian

        end
    end
end




fNr=fftn(Nr) ;
fNr=sqrt(fNr) ;

fm_flt = fm.*fNr ;
dv=sigma*real(ifftn(fm_flt));
% dv = dv/max(dv(:));

v(:,:)=dv(:,:,1);

v_inn=v([dx_inn+1:dx_inn+nx_inn],[dy_inn+1:dy_inn+ny_inn]);

figure,imagesc(v_inn'),colormap(flipud(jet)), axis square, colorbar

%% 
%Nyquist wavenumbers
knx = 1/2/dx; kny= 1/2/dy; knz = 1/2/dz;
%wavenumber intervals
dkx = 1/dx/nx; dky=1/dy/ny; dkz = 1/dz/nz;
%wavenumber vectors
ksx = -knx:dkx:knx-dkx; 
ksy = -kny:dky:kny-dky; 
ksz = -knz:dkz:knz-dkz;
%circular wavenumbers
ksx = ksx*2*pi; 
ksy = ksy*2*pi; 
ksz = ksz*2*pi;


theta = -pi/4; % 45 degrees

% Create the 2D rotation matrix
R2 = [cos(theta) -sin(theta);
     sin(theta)  cos(theta)];

for ix=1:nx
    for iy=1:ny
        for iz=1:nz
             
            ksx1 = R(1,1)*ksx(ix) + R(1,2)*ksy(iy);
            ksy1 = R(2,1)*ksx(ix) + R(2,2)*ksy(iy);
            
             ff(ix,iy,iz) = Ix*Iy/(8*pi^(3/2))*...
                 ( exp(- 1/4 *( ksx1.^2 * Ix^2 + ksy1.^2 * Iy^2 ) ) ) ;
        end
    end
end

fm_flt = fftshift(fm).*sqrt(ff) ;
dv=sigma*real(ifft2(ifftshift(fm_flt)));
v(:,:)=dv(:,:,1);

v_inn=v([dx_inn+1:dx_inn+nx_inn],[dy_inn+1:dy_inn+ny_inn]);

figure,imagesc(v_inn'),colormap(flipud(jet)), axis square, colorbar

%% anisotropic gaussian filter

[KX,KY] = ndgrid(ksx,ksy);
KX1 = R(1,1)*KX + R(1,2)*KY;
KY1 = R(2,1)*KX + R(2,2)*KY;
ax = 0.002;
ay = ax*5;
FLT = exp(-(KX1.^2/ax^2 + KY1.^2/ay^2)); %gaussian filter

Sv = fft2(v);

figure, imagesc(ksx,ksy,abs(fftshift(Sv))), 
axis square, axis([-0.1 0.1 -0.1 0.1])
%
v_flt = real(ifft2(ifftshift(FLT).*Sv));
figure, imagesc(ksx,ksy,FLT), axis square, axis([-0.1 0.1 -0.1 0.1])

v_inn_flt=v_flt([dx_inn+1:dx_inn+nx_inn],[dy_inn+1:dy_inn+ny_inn]);
figure, 
subplot(121),imagesc(x,y,v_inn), axis square
subplot(122),imagesc(x,y,v_inn_flt), axis square


