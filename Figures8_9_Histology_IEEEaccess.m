%%%%%%% Unmix de las imagenes de histología de la carpeta 00_HS_Preprocessed
%%%%%%%  Ya se cuentan con las mascaras y se van almacenado los end-members
%%%%%%%  extra, estimado los mapas de abundancias.
%%%%%%%  HSI transform
addpath('EBEAE')
addpath('Histology')
load('pure_stains_128bands.mat');
wl=wavelength;
Error=zeros(3,1);
%% Commun Parameters
initcond=5;
rho=0.01;
lambda=0.2;
epsilon=1e-4;
maxiter=30;
downsampling=0;
parallel=0;
display=0;
lm=0.1;
parameters=[initcond,rho,lambda,lm,epsilon,maxiter,downsampling,parallel,display];
%% 
 filename = 'preProcessed_0024322108_x20_C12.mat';
 load(filename);
 %%
 [Ny,Nx,Nz]=size(preProcessedImage);
  Z= replace_zeros_with_mean(reshape(preProcessedImage,Nx*Ny,Nz)'); %vectorizamos
  Zn = Z./repmat(sum(Z,1),Nz,1); %normalizar suma a 1
  load(strcat(filename(1,14:end-4),'.mat' ));
  filergb=strcat(filename(1,14:end-4),'.png' );
  rgb = imread(filergb);
  nwindex = reshape(BW,Nx*Ny,1)';
  Zsub = Zn(:,nwindex);
  Po2 = [Eosin, hematoxylin];

   [L, n]=size(Po2); %numero de end-members
   Po2n = (Po2./repmat(sum(Po2),Nz,1) );
%%
    tic
    N=3;
    [P3,A3,S3,Yh3]=ESSEAE(Zsub,N,parameters,Po2n);
    t=toc;
    Error(1) = norm(Yh3-Zsub,'fro')/norm(Zsub,'fro');
    disp(filename);
    disp(['Error with 3 EM, Hybrid 2 fixed 1 Unknown = ' num2str(Error(1))]);
    disp(['Computational Time (s) =' num2str(t)]);
    Afull = zeros(N,Nx*Ny);
    Afull(:,nwindex) = A3;
    a3_1 = reshape(Afull(1,:),Ny,Nx);
    a3_2 = reshape(Afull(2,:),Ny,Nx);
    a3_3 = reshape(Afull(3,:),Ny,Nx);
    maskedRGBImage(:,:,1)=single(rgb(:,:,1)).*single(BW);
    maskedRGBImage(:,:,2)=single(rgb(:,:,2)).*single(BW);
    maskedRGBImage(:,:,3)=single(rgb(:,:,3)).*single(BW);
    %%
    tic
    N=4;
    [P4,A4,S4,Yh4]=ESSEAE(Zsub,N,parameters,Po2n);
    t=toc;
    Error(2) = norm(Yh4-Zsub,'fro')/norm(Zsub,'fro');
    disp(filename);
    disp(['Error with 4 EM, Hybrid 2 fixed 2 Unknown = ' num2str(Error(2))]);
    disp(['Computational Time (s) =' num2str(t)]);
    Afull = zeros(N,Nx*Ny);
    Afull(:,nwindex) = A4;
    a4_1 = reshape(Afull(1,:),Ny,Nx);
    a4_2 = reshape(Afull(2,:),Ny,Nx);
    a4_3 = reshape(Afull(3,:),Ny,Nx);
    a4_4 = reshape(Afull(4,:),Ny,Nx); 
    maskedRGBImage(:,:,1)=single(rgb(:,:,1)).*single(BW);
    maskedRGBImage(:,:,2)=single(rgb(:,:,2)).*single(BW);
    maskedRGBImage(:,:,3)=single(rgb(:,:,3)).*single(BW);
%%
tic
    N=5;
    [P5,A5,S5,Yh5]=ESSEAE(Zsub,N,parameters,Po2n);
    t=toc;
    Error(3) = norm(Yh5-Zsub,'fro')/norm(Zsub,'fro');
    disp(filename);
    disp(['Error with 5 EM, Hybrid 2 fixed 3 Unknown = ' num2str(Error(3))]);
    disp(['Computational Time (s) =' num2str(t)]);
    Afull = zeros(N,Nx*Ny);
    Afull(:,nwindex) = A5;
    a5_1 = reshape(Afull(1,:),Ny,Nx);
    a5_2 = reshape(Afull(2,:),Ny,Nx);
    a5_3 = reshape(Afull(3,:),Ny,Nx);
    a5_4 = reshape(Afull(4,:),Ny,Nx);
    a5_5 = reshape(Afull(5,:),Ny,Nx); 
    maskedRGBImage(:,:,1)=single(rgb(:,:,1)).*single(BW);
    maskedRGBImage(:,:,2)=single(rgb(:,:,2)).*single(BW);
    maskedRGBImage(:,:,3)=single(rgb(:,:,3)).*single(BW);

%%
figure(2)
clf()
    subplot(221)
    hold on;
    grid on;
    plot(wl, Po2n(:,1), 'LineWidth', 2, 'color', 'r');
    plot(wl, Po2n(:,2), 'LineWidth', 2, 'color', 'b');
    axis tight;
    legend({'Eosin', 'Hematoxylin'} );
    title('Initial End-members');
    subplot(222)
    hold on;
    grid on;
    plot(wl, P3(:,1), 'LineWidth', 2, 'color', 'r');
    plot(wl, P3(:,2), 'LineWidth', 2, 'color', 'b');
    plot(wl, P3(:,3), '-.', 'LineWidth', 2, 'color', 'g');
    axis tight;
    legend({'Eosin', 'Hematoxylin', 'Unknown 1'} );
    title('Resulting End-members  N=3');
subplot(223)
    hold on;
    grid on;
    plot(wl, P4(:,1), 'LineWidth', 2, 'color', 'r');
    plot(wl, P4(:,2), 'LineWidth', 2, 'color', 'b');
    plot(wl, P4(:,3), '-.', 'LineWidth', 2, 'color', 'g');
    plot(wl, P4(:,4), '-.', 'LineWidth', 2, 'color', 'k');
    axis tight;
    legend({'Eosin', 'Hematoxylin', 'Unknown 1', 'Unknown 2'} );
    title('Resulting End-members N=4');
    subplot(224)
    hold on;
    grid on;
    plot(wl, P5(:,1), 'LineWidth', 2, 'color', 'r');
    plot(wl, P5(:,2), 'LineWidth', 2, 'color', 'b');
    plot(wl, P5(:,3), '-.', 'LineWidth', 2, 'color', 'g');
    plot(wl, P5(:,4), '-.', 'LineWidth', 2, 'color', 'k');
    plot(wl, P5(:,5), '-.', 'LineWidth', 2, 'color', 'm');
    axis tight;
    legend({'Eosin', 'Hematoxylin', 'Unknown 1', 'Unknown 2','Unknown 3'} );
    title('Resulting End-members N=5');

%%
figure(1)
clf()

% Create a 4x5 tiled layout for uniform sizing and control over spacing
t = tiledlayout(4, 5, 'TileSpacing', 'compact', 'Padding', 'compact');

% First row (1x2 layout, center the images)
nexttile(1); imshow(rgb); title('RGB'); axis 'off';  % Spans 2 tiles
nexttile(2); imshow(uint8(maskedRGBImage)); title('ROI'); axis 'off';  % Spans 2 tiles

% Second row (N=3) - center the images by manually spanning columns
nexttile(6); imshow(a3_1); title('Eosin abundance map'); axis 'off';
text(-0.1, 0.5, 'N=3', 'Units', 'normalized', 'FontSize', 14, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' );
nexttile(7); imshow(a3_2); title('Hematoxylin abundance map'); axis 'off';
nexttile(8); imshow(a3_3); title('Unknown 1 abundance map'); axis 'off';

% Third row (N=4)
nexttile(11); imshow(a4_1); title('Eosin abundance map'); axis 'off';
text(-0.1, 0.5, 'N=4', 'Units', 'normalized', 'FontSize', 14, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' );
nexttile(12); imshow(a4_2); title('Hematoxylin abundance map'); axis 'off';
nexttile(13); imshow(a4_3); title('Unknown 1 abundance map'); axis 'off';
nexttile(14); imshow(a4_4); title('Unknown 2 abundance map'); axis 'off';

% Fourth row (N=5)
nexttile(16); imshow(a5_1); title('Eosin abundance map'); axis 'off';
text(-0.1, 0.5, 'N=5', 'Units', 'normalized', 'FontSize', 14, 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' );
nexttile(17); imshow(a5_2); title('Hematoxylin abundance map'); axis 'off';
nexttile(18); imshow(a5_3); title('Unknown 1 abundance map'); axis 'off';
nexttile(19); imshow(a5_4); title('Unknown 2 abundance map'); axis 'off';
nexttile(20); imshow(a5_5); title('Unknown 3 abundance map'); axis 'off';

% Create a general colorbar that spans from the second to the fourth rows
cb = colorbar('Position', [0.92 0.35 0.02 0.45]);  % Manually position the colorbar
cb.Label.FontSize = 14;
colormap parula;


