addpath('Urban')
addpath('Histology')
%%
N=4;                % Number of End-members
Nsamples=64;       % Size of the Squared Image Nsamples x Nsamples 
SNR=0;             % Additive Gaussian noise
density=0.00;       % Density of sparse noise component
ModelType=0;        % 0 --> Linear Mixing Model and 5 --> Multilinear Model
EndMembersSynth=1;  % 0--> USGS Spectral Library Version & 1 --> Craneotomy 

%%
[Zlmm,~,~,~,~]=MatternGaussian_Sparse_Synth(SNR,density,0);
nRow=Nsamples;
nCol=Nsamples;
x=linspace(450,700,281);
%%
r=225;
g=113;
b=1;
HSI_lmm = reshape(Zlmm',nCol,nRow,281);
% Assuming 'HSI' is your hyperspectral image [height x width x 281]
blue_channel = HSI_lmm(:, :, b);   % Approx. 472 nm
green_channel = HSI_lmm(:, :, g);  % Approx. 532 nm
red_channel = HSI_lmm(:, :, r);   % Approx. 665 nm

% Normalize channels to [0, 1] for display
blue_channel = mat2gray(blue_channel);
green_channel = mat2gray(green_channel);
red_channel = mat2gray(red_channel);

% Combine into RGB
Synth_lmm = cat(3, red_channel, green_channel, blue_channel);
%%

Hist = imread('0024322108_x20_C12.png');
Urban = imread('UrbanRGB.png');
%%
Fs=14;
figure(1)
clf
subplot(221)
imshow(Synth_lmm);
title('(a)','FontSize', Fs)
axis('off')

subplot(222)
imshow(Urban);
title('(b)','FontSize', Fs)
axis 'off';

subplot(212)
imshow(Hist)
title('(c)','FontSize', Fs)
axis 'off'

