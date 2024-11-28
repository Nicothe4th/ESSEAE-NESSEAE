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
par  =[initcond,rho,lambda,epsilon,maxiter,downsampling,parallel,display];

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
  maskedRGBImage(:,:,1)=single(rgb(:,:,1)).*single(BW);
  maskedRGBImage(:,:,2)=single(rgb(:,:,2)).*single(BW);
  maskedRGBImage(:,:,3)=single(rgb(:,:,3)).*single(BW);


   [L, n]=size(Po2); %numero de end-members
   Po2n = (Po2./repmat(sum(Po2),Nz,1) );

%%
    tic
    N=2;
    [P2,A2,S2,Yh2,J2]=EBEAE(Zsub,N,par,Po2n,1);
    t=toc;
    Error(4)=norm(Yh2-Zsub,'fro')/norm(Zsub,'fro');
    disp(filename);
    disp(['Error with 2 EM, Supervised 2 fixed = ' num2str(Error(4))]);
    disp(['Computational Time (s) =' num2str(t)]);
    Afull = zeros(N,Nx*Ny);
    Afull(:,nwindex) = A2;
    a2_1 = reshape(Afull(1,:),Ny,Nx);
    a2_2 = reshape(Afull(2,:),Ny,Nx);
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

%%
figure(1)
clf()

% Create a 4x5 tiled layout for uniform sizing and control over spacing
t = tiledlayout(4, 6);

% First row (1x2 layout, center the images)
nexttile(1); 
hold on;
grid on;
plot(wl, Po2n(:,1), 'LineWidth', 2, 'color', 'r');
plot(wl, Po2n(:,2), 'LineWidth', 2, 'color', 'b');
axis tight;
xlabel('Wavelength (nm)')
ylabel('Normalized intensity')
legend({'Eosin', 'Hematoxylin'} );
text(0.5, 1.11, '(a)', 'Units', 'normalized', 'fontsize', 14, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' );

% Second row (N=2)
nexttile(2); imshow(a2_1); title('Eosin'); axis 'off';
text(1.1, 1.2, '(b)', 'Units', 'normalized', 'fontsize', 14, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' );
nexttile(3); imshow(a2_2); title('Hematoxylin'); axis 'off';


% Second row (N=3) - center the images by manually spanning columns
nexttile(7);
hold on;
grid on;
plot(wl, P3(:,1), 'LineWidth', 2, 'color', 'r');
plot(wl, P3(:,2), 'LineWidth', 2, 'color', 'b');
plot(wl, P3(:,3), '-.', 'LineWidth', 2, 'color', 'g');
axis tight;
legend({'Eosin', 'Hematoxylin', 'Unknown-1'} );
xlabel('Wavelength (nm)')
ylabel('Normalized intensity')
text(0.5, 1.1, '(c)', 'Units', 'normalized', 'FontSize', 14,  'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' );

nexttile(8); imshow(a3_1);  title('Eosin'); axis 'off';
nexttile(9); imshow(a3_2);  title('Hematoxylin'); axis 'off';
text(.5, 1.2, '(d)', 'Units', 'normalized', 'fontsize', 14, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' );
nexttile(10); imshow(a3_3);  title('Unknown-1');axis 'off';
 
% Third row (N=4)
nexttile(13);
hold on
grid on
plot(wl, P4(:,1), 'LineWidth', 2, 'color', 'r');
plot(wl, P4(:,2), 'LineWidth', 2, 'color', 'b');
plot(wl, P4(:,3), '-.', 'LineWidth', 2, 'color', 'g');
plot(wl, P4(:,4), '-.', 'LineWidth', 2, 'color', 'k');
hold off
axis tight;
legend({'Eosin', 'Hematoxylin', 'Unknown-1', 'Unknown-2'} );
xlabel('Wavelength (nm)')
ylabel('Normalized intensity')
text(0.5, 1.11, '(e)', 'Units', 'normalized', 'FontSize', 14, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' );

nexttile(14); imshow(a4_2);  title('Eosin'); axis 'off';
nexttile(15); imshow(a4_1); title('Hematoxylin'); axis 'off';
text(1.1, 1.2, '(f)', 'Units', 'normalized', 'fontsize', 14, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' );
nexttile(16); imshow(a4_3); title('Unknown-1'); axis 'off';
nexttile(17); imshow(a4_4); title('Unknown-2'); axis 'off'; 

% Fourth row (N=5)
nexttile(19)
hold on;
grid on;
plot(wl, P5(:,1), 'LineWidth', 2, 'color', 'r');
plot(wl, P5(:,2), 'LineWidth', 2, 'color', 'b');
plot(wl, P5(:,3), '-.', 'LineWidth', 2, 'color', 'g');
plot(wl, P5(:,4), '-.', 'LineWidth', 2, 'color', 'k');
plot(wl, P5(:,5), '-.', 'LineWidth', 2, 'color', 'm');
axis tight;
legend({'Eosin', 'Hematoxylin', 'Unknown-1', 'Unknown-2','Unknown-3'} );
xlabel('Wavelength (nm)')
ylabel('Normalized intensity')
text(0.5, 1.11, '(g)', 'Units', 'normalized', 'FontSize', 14,  'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' );



nexttile(20); imshow(a5_1); title('Eosin'); axis 'off';
nexttile(21); imshow(a5_2); title('Hematoxylin');  axis 'off';
nexttile(22); imshow(a5_3); title('Unknown-1');axis 'off';
text(.5, 1.2, '(h)', 'Units', 'normalized', 'fontsize', 14, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' );
nexttile(23); imshow(a5_4); title('Unknown-2'); axis 'off';
nexttile(24); imshow(a5_5); title('Unknown-3'); axis 'off';

% Create a general colorbar that spans from the second to the fourth rows
 cb = colorbar('Position', [0.96 0.1 0.01 0.8]);  % Manually position the colorbar
cb.Label.FontSize = 14;
colormap parula;
