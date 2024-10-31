%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NESSEAE for in-vivo brain VNIR Datasets
%
% Daniel Ulises Campos Delgado & Nicolas Mendoza-Chavarria
% FC-UASLP & ULPGC
% Version: Oct/2024
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 
addpath('EBEAE');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load in-vivo brain VNIR dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=10;
analysisType=1;    %% 0 --> Reflectance & 1 --> Absorbance 
plotJournal=0;     %% 0 --> Include titles & 1 --> No titles

     
file1='15C1';
titleL='(a) Op15C1';
file2='P3GoldenReference.mat';


load(['Craneotomy/VNIRimagesOp' file1 '.mat']);
Image=preProcessedImage; 
fileImage='preProcessedImage';
load(['Craneotomy/' file2]);
load('DatasetSynth/EMHbHbO2_3.mat');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust the number of spectral channels for initial datasets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Orignal wavelength range: 448 - 910

lmin=440;
lmax=910;
if lmin<440, lmin=440; end
if lmax>910, lmax=910; end
Is=find(wave2>lmin & wave2 <lmax);

IndexY=1:size(Image,1);
IndexX=1:size(Image,2);
IndexS=Is;
[Ny,Nx,Nz]=size(Image(IndexY,IndexX,IndexS));
Zz=reshape(Image(IndexY,IndexX,IndexS),Nx*Ny,Nz)';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transform from reflectance to absorbance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[L,K]=size(Zz);
if analysisType==0
    Z=Zz;
else
    if min(Zz(:)) < 1e-4
        NonZeroThreshold=1e-2;
    else
        NonZeroThreshold=0;
    end
    Z=-log10(Zz+NonZeroThreshold);
end

Z=(Z-min(Z(:)))/(max(Z(:))-min(Z(:)));
Z=Z./repmat(sum(Z,1),L,1);
wave2=wave2(IndexS);
Nz=size(Z,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load molar extinction coefficients of Hb & HbO2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I1=find(goldenStandardMap(:)==1);
I2=find(goldenStandardMap(:)==2);
I3=find(goldenStandardMap(:)==3);
I4=find(goldenStandardMap(:)==4);

P011=NFINDR(Z(:,I1),2);
P012=NFINDR(Z(:,I2),2);
P013=NFINDR(Z(:,I3),2);
P014=NFINDR(Z(:,I4),4);
F=size(P011,2)+size(P012,2)+size(P013,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Estimated End-members
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
subplot(1,3,1);
plot(wave2,Z(:,I1),'g:', 'LineWidth',0.25); 
hold on; grid on;
plot(wave2,P011,'k','linewidth',2);  axis tight;
title('(a)');
xlabel('wavelength (nm)'); ylabel('Intensity')
subplot(1,3,2);
plot(wave2,Z(:,I2),'r:', 'LineWidth',0.25);
hold on; grid on;
plot(wave2,P012,'k','LineWidth',2); axis tight;
xlabel('wavelength (nm)'); ylabel('Intensity')
title('(b)');
subplot(1,3,3);
plot(wave2,Z(:,I3),'b:', 'LineWidth',0.25);
hold on; grid on;
xlabel('wavelength (nm)'); ylabel('Intensity')
plot(wave2,P013,'k','LineWidth',2); axis tight;
title('(c)');


P01=[P011 P012 P013];
P001=[P011 P012 P013 P014];

if analysisType==0
    TypeD='Reflectance dataset';
else 
    TypeD='Absorbance dataset';
end

P0=P01./repmat(sum(P01,1),L,1);
P00=P001./repmat(sum(P001,1),L,1);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EBEAE Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                                                                                                                                ;             % Weight on Entropy of Abundances \in [0,1)
epsilon=1e-3;            % Threshold for convergence
maxiter=20;              % Maximum number of iteration in alternated least squares approach
parallel=1;              % Parallelization in abundance estimation process
display_iter=0;          % Display results of iterative optimization
initcond=5;
downsampling=0;
rho=0.25;                   % Weight on Regularization Term  >= 0
lambda=0.1; 
lm=0.05;
paramvecSN=[initcond,rho,lambda,lm,epsilon,maxiter,downsampling,parallel,display_iter];


tic;
[P1,A1,D1,S1,Zh1,V1,J1]=NESSEAE(Z,N,paramvecSN,P0);
Thebeae=toc;

rho=0.1;
lambda=0.2;
paramvecSN=[initcond,rho,lambda,epsilon,maxiter,downsampling,parallel,display_iter];

tic;
[P2,A2,D2,S2,Zh2,J2]=NEBEAE3(Z,size(P00,2),paramvecSN,P00,1);
Tnebeae=toc;
%%

analysisEBEAE='NESSEAE';
analysisEBEAE2='NEBEAE/Supervised';


disp('%%%%%%%%%%%%%%%%%%%%');
disp(['in-vivo Brain Dataset ' num2str(dataset)]);
disp(TypeD);
disp(titleL);
disp(['Number of fixed end-members=' num2str(size(P0,2))]);
disp(['Number of estimated end-members=' num2str(N-size(P0,2))]);
disp('%%%%%%%%%%%%%%%%%%%%');
disp(analysisEBEAE);

disp(['Estimation Error = ' num2str(norm(Zh1-Z,'fro')/norm(Z,'fro'))]);
disp(['Computation time = ' num2str(Thebeae)]);
disp(['Weight of Sparse Noise Energy =' num2str(norm(V1,'fro')/norm(Z,'fro'))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Abundance Maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Aa=[];
Ac=zeros(4,K);
legendP={'NT 1','NT 2','TT 1', 'TT 2','HT 1','HT 2'};
legendC='Normal, Tumor, Blood Vessels, Background';
index=1;
F=size(P011,2)+size(P012,2)+size(P013,2);
indexC=[size(P011,2),size(P012,2),size(P013,2),N-F];
for i=1:4
    if i>=1
        Ac(i,:)= sum(A1(index:index+indexC(i)-1,:));
    else
        Ac(i,:)=A1(index:index+indexC(i)-1,:);
    end
    Aa=[Aa ones(Ny,5) reshape(Ac(i,:),Ny,Nx)];
    if i>3
        for j=1:indexC(i)
            str=[' Estimation ' num2str(j) ','];
            legendP(sum(indexC(1:3))+j)={str};
        end
    end
    index=index+indexC(i);
end

[tt,maxC]=max(Ac,[],1);

figure;
imagesc(2*rgbCropped)
if plotJournal==0, title('(a) RGB image'); else title('(a)'); end
figure;
imagesc(reshape(goldenStandardMap, Ny,Nx))
mycolors0=[1 1 1;0 1 0; 1 0 0; 0 0 1; 0 0 0];
colormap(mycolors0)
if plotJournal==0, title('(b) Ground-truth pixels'); else title('(b)'); end
figure;
imagesc(reshape(maxC,Ny,Nx))
mycolors=[0 1 0; 1 0 0; 0 0 1; 0 0 0];
colormap(mycolors)
if plotJournal==0, title('(c) Classified image (NESSEAE)'); else title('(c)'); end
figure;
Q=quantile([D1;D2],[0.25 0.75]);
imagesc(reshape(D1,Ny,Nx),[Q(1)-1.5*(Q(2)-Q(1)) Q(2)+1.5*(Q(2)-Q(1))]); 
colorbar('eastoutside');
if plotJournal==0, title('(d) Nonlinear interaction levels (NESSEAE)'); else title('(d)'); end


%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%

figure;
subplot(2,2,1);
imshow(2*rgbCropped); title({analysisEBEAE, titleL})
subplot(2,1,2);
imagesc(Aa,[0,1]); colorbar('eastoutside');
title({'(c) Abundance Maps', legendC});
subplot(2,2,2)
plot(wave2,P1,'linewidth',2); axis tight; legend(legendP); 
xlabel('Wavelength (nm)'); grid on;
title({TypeD, '(b) End-members'}); 

figure;
subplot(1,2,1);
imagesc(V1);
xlabel('Spatial location'); ylabel('Spectral channel'); colormap('pink'); colorbar('southoutside')
if plotJournal==0, title('(a) Sparse noise matrix'); else title('(a)'); end
subplot(1,2,2)
bar(wave2,sum(V1,2)); xlabel('Wavelength (nm)'); ylabel('Accumulated noise intensities');
axis tight
if plotJournal==0, title('(b) Spectral noise'); else title('(b)'); end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('%%%%%%%%%%%%%%%%%%%%');
disp(analysisEBEAE2);
disp(['Estimation Error = ' num2str(norm(Zh2-Z,'fro')/norm(Z,'fro'))]);
disp(['Computation time = ' num2str(Tnebeae)]);

Ac2=zeros(4,K);
Aa2=[];
index=1;
indexC2=[size(P011,2) size(P012,2) size(P013,2) size(P014,2)];
for i=1:4
    if i>=1
        Ac2(i,:)= sum(A2(index:index+indexC2(i)-1,:));
    else
        Ac2(i,:)=A2(index:index+indexC2(i)-1,:);
    end
    Aa2=[Aa2 ones(Ny,5) reshape(Ac2(i,:),Ny,Nx)];
    index=index+indexC2(i);
end

[tt,maxC2]=max(Ac2,[],1);

figure
imagesc(reshape(maxC2,Ny,Nx))
mycolors=[0 1 0; 1 0 0; 0 0 1; 0 0 0];
colormap(mycolors)
if plotJournal==0, title('(e) Classified image (NEBEAE-SN/Supervised)'); else title('(e)'); end
figure
imagesc(reshape(D2,Ny,Nx),[Q(1)-1.5*(Q(2)-Q(1)) Q(2)+1.5*(Q(2)-Q(1))]); 
colorbar('eastoutside');
if plotJournal==0, title('(f) Nonlinear interaction levels (NEBEAE-SN/Supervised)'); else title('(f)'); end


%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%

figure;
subplot(2,2,1);
imshow(2*rgbCropped); title({analysisEBEAE2, titleL})
subplot(2,1,2);
imagesc(Aa2,[0,1]); colorbar('eastoutside');
title({'(c) Abundance Maps', legendC});
subplot(2,2,2)
plot(wave2,P2,'linewidth',2); axis tight; legend(legendP); 
xlabel('Wavelength (nm)'); grid on;
title({TypeD, '(b) End-members'}); 

P014h=P1(:,F+1:N);
figure; 
plot(wave2,P014,':','LineWidth',2);
hold on; grid on;
plot(wave2, P014h, 'LineWidth',2 )
axis tight
xlabel('Wavelength (nm)');
ylabel('Intensity');
legend('Background-1','Background-2','Background-3','Background-4','Unknown-1','Unknown-2','Unknown-3','Unknown-4')


Ii=find(not(goldenStandardMap==0));
T1=goldenStandardMap(Ii);
O1=maxC(Ii)';
O2=maxC2(Ii)';
disp('Confusion Matrix (NESSEAE)')
CM1=confusionmat(T1,O1)
disp('Confusion Matrix (NEBEAE-SN/Supervised)')
CM2=confusionmat(T1,O2)
disp(['Precision (NESSEAE) = ' num2str(sum((T1-O1)==0)/length(T1))]);
display(['Precision (NEBEAE-SN/Supervised) =' num2str(sum((T1-O2)==0)/length(T1))]);
display(['Error in estimated end-members = ' num2str(errorendmembers(P014,P014h))]);
display(['SAM in estimated end-members = ' num2str(errorSAM(P014,P014h))]);