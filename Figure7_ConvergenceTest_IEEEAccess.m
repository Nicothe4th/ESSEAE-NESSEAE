%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Analysis of Convergence in ESSEAE & NESSEAE with Synthetic Datasets
%
% Daniel Ulises Campos-Delgado & Nicolas Mendoza Chavarria
% October/2024
% UASLP-ULPGC
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;
clc;
addpath('EBEAE');
SNR=30; 
density=0.01;

N=4;                % Number of End-members
Nsamples=64;
nCol=Nsamples;
nRow=Nsamples;
Rep=30;
%%
initcond=6;             % Initial condition of end-members matrix: 6 (VCA) and 8 (SISAL).
rho=0.1;               % Similarity weight in end-members estimation
lambda=0.1;             % Entropy weight for abundance estimation
epsilon=1e-3;
maxiter=20;
parallel=1;
downsampling=0.0;       % Downsampling in end-members estimation
display_iter=0;            % Display partial performance in BEAE
lm=0.1;
paramvec=[initcond,rho,lambda, lm, epsilon,maxiter,downsampling,parallel,display_iter];
%%
%
jl1=zeros(Rep,maxiter);
jl2=zeros(Rep,maxiter);
jl3=zeros(Rep,maxiter);
jm1=zeros(Rep,maxiter);
jm2=zeros(Rep,maxiter);
jm3=zeros(Rep,maxiter);

%% LINEAR PART
ModelType=0;        % 0 --> Linear Mixing Model and 5 --> Multilinear Model
EndMembersSynth=1;  % 0--> USGS Spectral Library Version & 1 --> Craneotomy 
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp('Synthetic Datasets');
    disp('Absorbance HSI');
    disp('Semi-Supervised Unmixing Estimation');
    disp(['SNR =' num2str(SNR) ' dB']);
    disp(['density =' num2str(density) ]);
    disp(['Number of end-members=' num2str(N)]);
    for j=1:Rep 
        [Z,P0,A0,V0,D0]=MatternGaussian_Sparse_Synth(SNR,density,ModelType);
        titleL='Synthetic Dataset with Absorbance Spectral Response of Hb, HbO2, Fat and Water';
        [L,K]=size(Z);
        P0=P0./sum(P0);
        disp(['Iteration=' num2str(j)])
        for n=1:3
                rNs =  sort(randperm(4, n));
                Pu=P0(:,rNs);
                [P1,A1,S1,Zh1,V1,J1]=HybridEBEAESN(Z,N,paramvec,Pu);
                if n==1
                    jl1(j,:)=J1;
                end
                if n==2
                    jl2(j,:)=J1;
                end
                if n==3
                    jl3(j,:)=J1;
                end
        end
    end
%% MULTI-LINEAR PART
ModelType=5;        % 0 --> Linear Mixing Model and 5 --> Multilinear Model
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp('Synthetic Datasets');
    disp('Absorbance HSI');
    disp('Semi-Supervised Unmixing Estimation');
    disp(['SNR =' num2str(SNR) ' dB']);
    disp(['density =' num2str(density) ]);
    disp(['Number of end-members=' num2str(N)]);
    for j=1:Rep 
        [Z,P0,A0,V0,D0]=MatternGaussian_Sparse_Synth(SNR,density,ModelType);
        titleL='Synthetic Dataset with Absorbance Spectral Response of Hb, HbO2, Fat and Water';
        [L,K]=size(Z);
        P0=P0./sum(P0);
        disp(['Iteration=' num2str(j)])
        for n=1:3
            rNs =  sort(randperm(4, n));
                Pu=P0(:,rNs);
                [P1,A1,D1,S1,Zh1,V1,J1]=HybridNEBEAESN(Z,N,paramvec,Pu);
                if n==1
                    jm1(j,:)=J1;
                end
                if n==2
                    jm2(j,:)=J1;
                end
                if n==3
                    jm3(j,:)=J1;
                end
        end
    end
%%
meanl1=mean(jl1,1);
stdl1=std(jl1,1);

meanl2=mean(jl2,1);
stdl2=std(jl2,1);

meanl3=mean(jl3,1);
stdl3=std(jl3,1);

meanm1=mean(jm1,1);
meanm2=mean(jm2,1);
meanm3=mean(jm3,1);

stdm1=std(jm1,1);
stdm2=std(jm2,1);
stdm3=std(jm3,1);
%%
x = 1:maxiter;  % Assuming x-axis values from 1 to 20

figure(1);
clf;

% First subplot for meanl and stdl
subplot(2,1,1);  % 2 rows, 1 column, this is the first plot
hold on;

% Plot meanl1 with shading for stdl1
fill([x fliplr(x)], [meanl1+stdl1 fliplr(meanl1-stdl1)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');  % Shaded area
h1 = plot(x, meanl1, 'r', 'LineWidth', 2);  % Line plot for meanl1

% Plot meanl2 with shading for stdl2
fill([x fliplr(x)], [meanl2+stdl2 fliplr(meanl2-stdl2)], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');  % Shaded area
h2 = plot(x, meanl2, 'g', 'LineWidth', 2);  % Line plot for meanl2

% Plot meanl3 with shading for stdl3
fill([x fliplr(x)], [meanl3+stdl3 fliplr(meanl3-stdl3)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');  % Shaded area
h3 = plot(x, meanl3, 'b', 'LineWidth', 2);  % Line plot for meanl3

title('(a) ESSEAE','FontSize',11, 'FontWeight','normal');
xlabel('iteration');
ylabel('$\frac{|J^{l}-J^{l+1}|}{J^l}\; (53)$','Interpreter', 'latex');
% Add legend only for the mean lines (h1, h2, h3)
legend([h1, h2, h3], {'U=1', 'U=2', 'U=3'});
grid on;
hold off;

% Second subplot for meanm and stdm
subplot(2,1,2);  % This is the second plot
hold on;

% Plot meanm1 with shading for stdm1
fill([x fliplr(x)], [meanm1+stdm1 fliplr(meanm1-stdm1)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');  % Shaded area
h4 = plot(x, meanm1, 'r', 'LineWidth', 2);  % Line plot for meanm1

% Plot meanm2 with shading for stdm2
fill([x fliplr(x)], [meanm2+stdm2 fliplr(meanm2-stdm2)], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');  % Shaded area
h5 = plot(x, meanm2, 'g', 'LineWidth', 2);  % Line plot for meanm2

% Plot meanm3 with shading for stdm3
fill([x fliplr(x)], [meanm3+stdm3 fliplr(meanm3-stdm3)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');  % Shaded area
h6 = plot(x, meanm3, 'b', 'LineWidth', 2);  % Line plot for meanm3

title('(b) NESSEAE','FontSize',11, 'FontWeight','normal');xlabel('iteration');
ylabel('$\frac{|J^{l}-J^{l+1}|}{J^l}\; (53)$','Interpreter', 'latex');
% Add legend only for the mean lines (h4, h5, h6)
legend([h4, h5, h6], {'U=1', 'U=2', 'U=3'});
grid on;
hold off;
