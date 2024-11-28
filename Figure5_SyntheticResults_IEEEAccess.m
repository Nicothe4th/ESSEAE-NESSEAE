%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Experimental Validation of ESSEAE & NESSEAE with Synthetic Dataset
% 
%
% OCT/2024
% JNMC-DUCD-UASLP
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 
addpath('EBEAE');
SNR=30;   
density=0.01;

N=4;                % Number of End-members
Nsamples=64;
nCol=Nsamples;
nRow=Nsamples;
ModelType=0;        % 0 --> Linear Mixing Model and 5 --> Multilinear Model
EndMembersSynth=1;  % 0--> USGS Spectral Library Version & 1 --> Craneotomy 

initcond=6;             % Initial condition of end-members matrix: 6 (VCA) and 8 (SISAL).
rho=0.1;               % Similarity weight in end-members estimation
lambda=0.1;             % Etropy weight for abundance estimation
epsilon=1e-3;
maxiter=20;
parallel=1;
downsampling=0.0;       % Downsampling in end-members estimation
display_iter=0;            % Display partial performance in BEAE
lm=0.1;
paramvec=[initcond,rho,lambda, lm, epsilon,maxiter,downsampling,parallel,display_iter];
n=3;                    % Three unknown end-members in ESSEAE & NESSEAE
rNs =  sort(randperm(4, 4-n));  % Random selection of unknown end-members

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('Synthetic Datasets');
disp('Absorbance HSI');
disp('ESSEAE');
disp(['SNR =' num2str(SNR) ' dB']);
disp(['density =' num2str(density) ]);
disp(['Number of unknown end-members=' num2str(n)]);
    
[Z,P0,A0,V0,D0]=MatternGaussian_Sparse_Synth(SNR,density,ModelType);
titleL='Synthetic Dataset with Absorbance Spectral Response of Hb, HbO2, Fat and Water';
[L,K]=size(Z);
        
%%
Pu=P0(:,rNs);
tic;
[P1,A1,S1,Zh1,V1,J1]=ESSEAE(Z,N,paramvec,Pu);
tesseae=toc;
disp(['E_Z =' num2str(norm(Zh1-Z,'fro')/norm(Z,'fro'))])
disp(['E_a =' num2str(errorabundances(A0,A1))]);
disp(['E_p =' num2str(errorendmembers(P0,P1))]);
disp(['E_{SAM} =' num2str(errorSAM(P0,P1))]);
disp(['Computational time =' num2str(tesseae)]);

wavelength=linspace(450,700,size(P0,1));
figure;
subplot(3,2,3)
semilogy(wavelength,P1,'linewidth',2); grid on;
legend('Fixed','Unknown-1','Unknown-2','Unknown-3','FontSize',12);
xlabel('wavelength (nm)','FontSize',12); axis([min(wavelength) max(wavelength) min(P1(:)) max(P1(:))]);
ylabel('Normalized absorbance','FontSize',12)
title('(c)','FontSize',14,'FontWeight','normal')

for i=1:4
    subplot(3,8,12+i);
    imagesc(reshape(A1(i,:),nRow,nCol),[0,1]); axis off;
    xlabel('pixels'); ylabel('pixels'); 
    if i>1 
        title(['Unknown-' num2str(i-1)],'FontWeight','normal');
    else
        title('Fixed','FontWeight','normal')
    end
end

annotation('textbox', [0.7, 0.6 0.1, 0.1], 'string', '(d)', 'FontSize', 14, 'FontWeight', 'normal', 'EdgeColor', 'none');


%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%
ModelType=5;        % 0 --> Linear Mixing Model and 5 --> Multilinear Model
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('Synthetic Datasets');
disp('Absorbance HSI');
disp('NESSEA');
disp(['SNR =' num2str(SNR) ' dB']);
disp(['density =' num2str(density) ]);
disp(['Number of unknown end-members=' num2str(n)]);
    
[Z,P0,A0,V0,D0]=MatternGaussian_Sparse_Synth(SNR,density,ModelType);
titleL='Synthetic Dataset with Absorbance Spectral Response of Hb, HbO2, Fat and Water';
[L,K]=size(Z);
        

Pu=P0(:,rNs);
tic;
[P2,A2,D2,S2,Zh2,V2,J2]=NESSEAE(Z,N,paramvec,Pu);
tnesseae=toc;
disp(['E_Z =' num2str(norm(Zh2-Z,'fro')/norm(Z,'fro'))])
disp(['E_a =' num2str(errorabundances(A0,A2))]);
disp(['E_p =' num2str(errorendmembers(P0,P2))]);
disp(['E_{SAM} =' num2str(errorSAM(P0,P2))]);
disp(['(Mean,STD) NIL =' num2str([mean(D2), std(D2)])]);
disp(['Computational time =' num2str(tnesseae)]);

wavelength=linspace(450,700,size(P0,1));
subplot(3,2,5)
semilogy(wavelength,P2,'linewidth',2); grid on;
legend('Fixed','Unknown-1','Unknown-2','Unknown-3','FontSize',12);
xlabel('wavelength (nm)','FontSize',12); axis([min(wavelength) max(wavelength) min(P1(:)) max(P1(:))]);
ylabel('Normalized absorbance','FontSize',12)
title('(e)','FontSize',14,'FontWeight','normal')

for i=1:4
    subplot(3,8,20+i);
    imagesc(reshape(A2(i,:),nRow,nCol),[0,1]); axis off;
    xlabel('pixels'); ylabel('pixels'); 
    if i>1 
        title(['Unknown-' num2str(i-1)],'FontWeight','normal');
    else
        title('Fixed','FontWeight','normal')
    end
end

annotation('textbox', [0.7, 0.3, 0.1, 0.1], 'string', '(f)', 'FontSize', 14, 'FontWeight', 'normal', 'EdgeColor', 'none');


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%

subplot(3,2,1)
semilogy(wavelength,P0,'linewidth',2); grid on;
legend('1','2','3','4','FontSize',12);
xlabel('wavelength (nm)','FontSize',12); axis([min(wavelength) max(wavelength) min(P0(:)) max(P0(:))]);
ylabel('Normalized absorbance','FontSize',12)
title('(a)','FontSize',14,'FontWeight','normal')

for i=1:4
    subplot(3,8,4+i);
    imagesc(reshape(A0(i,:),nRow,nCol),[0,1]); axis off;
    xlabel('pixels'); ylabel('pixels'); 
    title(['End-member-' num2str(i)],'FontWeight','normal');
end

subplot(3,8,6); 
annotation('textbox', [0.7, 0.9 0.1, 0.1], 'string', '(b)', 'FontSize', 14, 'FontWeight', 'normal', 'EdgeColor', 'none');
