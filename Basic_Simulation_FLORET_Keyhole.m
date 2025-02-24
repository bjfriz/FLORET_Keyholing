%% Clean up workspace
clear;clc;close all;fclose all;
%% Find the path that this code is on to make things self-contained
parent_path = which('Basic_Simulation_FLORET_Keyhole');
idcs = strfind(parent_path,filesep);%determine location of file separators
parent_path = parent_path(1:idcs(end)-1);%remove file

%% Need a phantom to work with
ImageSize = 96; %Size of ventilation images that we collect
PhantSize = ImageSize; %Size of initial phantom (for better simulation)
%Create a phantom much larger than the image size
MyPhantom = Keyhole_Tools.phantom3d(PhantSize);
%This phantom has lots of areas of low signal intensity. Adjust a bit here:
MyPhantom(MyPhantom > 0.19 & MyPhantom < 0.22) = 0.8;
MyPhantom(MyPhantom < 0) = 0.5;

%% Now, create image signal decay - RF induced Decay
FA = 1;
%Uniform Flip Angle
%anglemap = ones(PhantSize,PhantSize,PhantSize) *FA;
% % % % % % %Sharply Varying flip angles
% BoxSize = PhantSize/2;
% Cube = ones(BoxSize,BoxSize,BoxSize);
% Oct1 = 0.8*Cube;
% Oct2 = 0.9*Cube;
% Oct3 = 1.0*Cube;
% Oct4 = 1.1*Cube;
% anglemap = cat(1,cat(2,cat(3,Oct1,Oct2),cat(3,Oct3,Oct4)),cat(2,cat(3,Oct3,Oct4),cat(3,Oct1,Oct2)));
% % % % % % %Smoothly Varying flip angles
 x = 1:PhantSize;
 filt = (FA^(1/3))*exp(-((x-PhantSize/2).^2)/2*(0.03).^2);
 [X,Y,Z] = meshgrid(filt,filt,filt);
 anglemap = X.*Y.*Z;

%% Now, create image signal decay - T1 induced Decay
T1 = 30;
%Uniform Flip Angle
%T1map = ones(PhantSize,PhantSize,PhantSize) * T1;
% % % % % %Sharply Varying T1
% BoxSize = PhantSize/2;
% Cube = ones(BoxSize,BoxSize,BoxSize);
% Oct1 = 22*Cube;
% Oct2 = 26*Cube;
% Oct3 = 30*Cube;
% Oct4 = 34*Cube;
% T1map = cat(1,cat(2,cat(3,Oct1,Oct2),cat(3,Oct3,Oct4)),cat(2,cat(3,Oct3,Oct4),cat(3,Oct1,Oct2)));
% % % % % %Smoothly Varying T1
 x = 1:PhantSize;
 filt = (T1^(1/3))*exp(-((x-PhantSize/2).^2)/2*(0.03).^2);
 [X,Y,Z] = meshgrid(filt,filt,filt);
 T1map = X.*Y.*Z;
%% Create Signal Decay Matrix
TR = 0.008; %8 ms between gaseous excitations
SigDecay = cosd(anglemap).*exp(-TR./T1map);
%% Simulate Decay and Sample Image - First need trajectories
% Use the trajectories that are actually used in imaging
traj_file = fullfile(parent_path,'Traj_files','Vent_GasExchange_20210819_Traj.dat');
traj_twix = Keyhole_Tools.mapVBVD(traj_file);
Im_file = fullfile(parent_path,'Traj_files','Noise_File.dat');
Im_twix = Keyhole_Tools.mapVBVD(Im_file);
traj = Keyhole_Tools.spiral_coords_from_dat(traj_twix,Im_twix);

%% Sampling data is a problem when datapoints aren't on a grid. 
%Need to change trajectory points to pixel units and round to nearest grid
%point
pix_traj = traj*ImageSize + PhantSize/2 + 1; 
pix_traj = round(pix_traj);
pix_traj(pix_traj>PhantSize) = PhantSize;
%since trajectories were rounded to nearest point, need to go back to
%unitless trajectories with the rounded points:
reco_traj = (pix_traj - PhantSize/2 -1)/ImageSize;
%if interested, can visualize trajectories:
%unaltered
Keyhole_Tools.disp_traj(traj,false,false,100);
%Rounded
Keyhole_Tools.disp_traj(reco_traj,false,false,100);
%% Decay and sample data
%Preallocate memory
raw = zeros(size(traj,2),size(traj,3));
%Loop through all spiral arms
for i = 1:size(traj,3)
   % tic
    %Decay data using pre-derived flip angle and T1 matrices
    MyPhantomDecay = MyPhantom.*sind(anglemap).*SigDecay.^(i-1);
    %Convert to k-space for image sampling
    MyPhantomKSpace = fftshift(ifftn(ifftshift(MyPhantomDecay)));
    %Sample data using pre-made pixel trajectory points
    my_indices = sub2ind(size(MyPhantomKSpace),squeeze(pix_traj(1,:,i))',squeeze(pix_traj(2,:,i))',squeeze(pix_traj(3,:,i))');
    raw(:,i) = MyPhantomKSpace(my_indices);
   % toc
end
%% Sanity check of Raw data
figure('Name','Raw Data Sanity Check')
imagesc(abs(raw));
colormap(gray);

%% Reconstruct Image
%Image recon requires data in columns - reshape here
reco_raw = reshape(raw,1,[])';
reco_traj_C = Keyhole_Tools.column_traj(reco_traj);

%Make sure no points are larger than 0.5 - otherwise, recon will hang
rad = sqrt(reco_traj_C(:,1).^2+reco_traj_C(:,2).^2+reco_traj_C(:,3).^2);
reco_traj_C(rad>0.5,:) = [];
reco_raw(rad>0.5) = [];

%Call Image reconstruction - Want to make sure that this works before
%starting the keyhole process.
Full_Im = AllinOne_Recon.base_floret_recon(ImageSize,reco_raw,reco_traj_C);
imslice(abs(Full_Im),'Original Image Recon');

%% Keyhole Data
%Number of keys, 2
NumKeys = 2;
maxRad = .025;

%Number of Projections(spirals), and how many in each key
NumSpirals = 1238;
SpiralsInKey = NumSpirals/NumKeys; %619 spirals in each key

%calculate the radius of points in projection rad = sqrt(...)
rad = squeeze(sqrt(reco_traj(1,:,:).^2+reco_traj(2,:,:).^2+reco_traj(3,:,:).^2));

%seperate keys from key hole
keyhole = raw;
keyhole(rad<maxRad) = nan;

%Create the keys
key1 = keyhole;
key2 = keyhole;
for i = 1:(NumSpirals/2)
    key1(rad(:,i)<maxRad,1:(NumSpirals/2)) = raw(rad(:,i)<maxRad,1:(NumSpirals/2));
    key2(rad(:,i+NumSpirals/2)<maxRad,(NumSpirals/2+1):NumSpirals) = raw(rad(:,i+NumSpirals/2)<maxRad,(NumSpirals/2+1):NumSpirals);
end

%Average k0 of key
key1av = mean(abs(raw(1,1:(NumSpirals/2))));
key2av = mean(abs(raw(1,(SpiralsInKey+1):NumSpirals)));

%Average keyhole (high frequency data) with each keys k0
%keyholepoint * (average k0/k0 for that column)
avdKey1 = key1;
avdKey2 = key2;
for i = 1:NumSpirals
    avdKey1(rad(:,i)>maxRad,i) = keyhole(rad(:,i)>maxRad,i)*(key1av/raw(1,i));
    avdKey2(rad(:,i)>maxRad,i) = keyhole(rad(:,i)>maxRad,i)*(key2av/raw(1,i));
end

%reconstruct the keys
%Image recon requires data in columns - reshape here
reco_Key1 = reshape(avdKey1,1,[])';
reco_Key2 = reshape(avdKey2,1,[])';
reco_traj_C1 = Keyhole_Tools.column_traj(reco_traj);
reco_traj_C2 = Keyhole_Tools.column_traj(reco_traj);

%Make sure there are no NaN points
reco_traj_C1(isnan(avdKey1),:) = [];
reco_Key1(isnan(avdKey1)) = [];
reco_traj_C2(isnan(avdKey2),:) = [];
reco_Key2(isnan(avdKey2)) = [];

Full_Key1 = AllinOne_Recon.base_floret_recon(ImageSize,reco_Key1,reco_traj_C1);
imslice(abs(Full_Key1),'Key 1 image recon');
clim([0 max(abs(Full_Key1(:)))]);

Full_Key2 = AllinOne_Recon.base_floret_recon(ImageSize,reco_Key2,reco_traj_C2);
imslice(abs(Full_Key2),'Key 2 image recon');
clim([0 max(abs(Full_Key1(:)))]);
%figure(5);montage(abs(Full_Key2),Indices=25:75);clim([0 max(abs(Full_Key1(:)))]);

%% Finding c1
c1 = (Full_Key2./Full_Key1).^(1/(NumSpirals/2));

figure('Name','Calculated decay');montage(real(c1));clim([.999 max(abs(SigDecay(:)))]);
colormap("jet")
imslice(abs(c1),'Calculated decay');
clim([.999 max(abs(SigDecay(:)))]);
colormap("jet")

mask = MyPhantom;mask(MyPhantom>0) = 1;
SigDecay = SigDecay.*mask;
figure('Name','Expected Decay');montage(real(SigDecay));clim([0.999 max(abs(SigDecay(:)))])
colormap("jet")
imslice(abs(SigDecay),'Expected Decay');
clim([.999 max(abs(SigDecay(:)))]);
colormap("jet")

%% Finding the fractional signal Attenuation
attenuation = 1 - (1/NumSpirals)*((1-(c1.^NumSpirals))./(1-c1));
figure(31);montage(real(attenuation));
%% reconstruct unattenuated image (recovered from original image and key images)
% image = original image/ (1-attenuation)
% image = original image * ((number spirals *(1-c1))/(1-(c1^number spirals)))
unattImage1 = Full_Im./(1-attenuation);
imslice(abs(unattImage1),'Corrected Image');
clim([0 max(abs(Full_Im(:)))]);
