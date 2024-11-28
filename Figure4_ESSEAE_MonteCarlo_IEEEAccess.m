%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Experimental Evaluation of ESSEAE with Synthetic Dataset
% 
%
% October/2024
% JNMC-DUCD-UASLP
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 
addpath('EBEAE');
sSNR=[40 35 30];   
pDensity=[0.005 0.0075 0.01];

N=4;                % Number of End-members
Nsamples=64;
nCol=Nsamples;
nRow=Nsamples;
ModelType=0;        % 0 --> Linear Mixing Model and 5 --> Multilinear Model
Rep=50;
%%
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
%%
ResultsYh =zeros(length(sSNR),Rep,5);
ResultsAh =zeros(length(sSNR),Rep,5);
ResultsPh =zeros(length(sSNR),Rep,5);
ResultsTh =zeros(length(sSNR),Rep,5);
ResultsPh2=zeros(length(sSNR),Rep,5);

%%
for index=1:length(sSNR)

    SNR=sSNR(index);
    density=pDensity(index);
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
        disp(['Iteration=' num2str(j)])
        
        for n=0:4
            if n==0
                tic;
                [P2,A2,S2,Zh2,V2,J2]=EBEAESN(Z,N,paramvec,P0,1);
                tsup=toc;
                ResultsYh(index,j,1)=norm(Zh2-Z,'fro')/norm(Z,'fro');
                ResultsAh(index,j,1)=errorabundances(A0,A2);
                ResultsPh(index,j,1)=errorendmembers(P0,P2);
                ResultsPh2(index,j,1)=errorSAM(P0,P2);
                ResultsTh(index,j,1)=tsup;
            elseif n>0 && n<4
                rNs =  sort(randperm(4, 4-n));
                Pu=P0(:,rNs);
                tic;
                [P1,A1,S1,Zh1,V1,J1]=ESSEAE(Z,N,paramvec,Pu);
                tesseae=toc;
                ResultsYh(index,j,n+1)=norm(Zh1-Z,'fro')/norm(Z,'fro');
                ResultsAh(index,j,n+1)=errorabundances(A0,A1);
                ResultsPh(index,j,n+1)=errorendmembers(P0,P1);
                ResultsPh2(index,j,n+1)=errorSAM(P0,P1);
                ResultsTh(index,j,n+1)=tesseae;
            elseif n==4
                [P3,A3,S3,Zh3,V3,J3]=EBEAESN(Z,N,paramvec);
                tblind=toc;
                ResultsYh(index,j,5)=norm(Zh3-Z,'fro')/norm(Z,'fro');
                ResultsAh(index,j,5)=errorabundances(A0,A3);
                ResultsPh(index,j,5)=errorendmembers(P0,P3);
                ResultsPh2(index,j,5)=errorSAM(P0,P3);
                ResultsTh(index,j,5)=tblind;
            end
        end
        
    end
end
%%
figure
clf
subplot(321)
cla()  % Clear the current axes

a = squeeze(ResultsYh(1,:,:));  % Data for the first noise scenario
b = squeeze(ResultsYh(2,:,:));  % Data for the second noise scenario
c = squeeze(ResultsYh(3,:,:));  % Data for the third noise scenario

hold on
nConditions = size(a, 2);  % Number of conditions (should be 5, based on ResultsYh)

% Offset positions to avoid overlap: Create enough space between groups
positions_a = (1:nConditions) - 0.2;  % Slightly left
positions_b = (1:nConditions);        % Centered
positions_c = (1:nConditions) + 0.2;  % Slightly right

% Plot the boxplots with specified positions to avoid overlap
boxplot(a, 'Colors', 'r', 'positions', positions_a, 'widths', 0.1,'Notch','on','symbol','+');
boxplot(b, 'Colors', 'b', 'positions', positions_b, 'widths', 0.1,'Notch','on','symbol','+');
boxplot(c, 'Colors', 'g', 'positions', positions_c, 'widths', 0.1,'Notch','on','symbol','+');


% Customize x-tick labels to display 0 to 4
xticks(1:nConditions);         % Adjust the ticks to correspond to the middle set (b)
xticklabels(0:nConditions-1);  % X labels as 0 to 4

% Add title and labels
title('(a) \it E_Z','FontSize', 11, 'FontWeight','normal');

% Adjust the y-axis to fit the data range
allData = [a(:); b(:); c(:)];  % Combine all data into one array
ylim([min(allData) max(allData)]);  % Set the y-axis limits based on data range

% Create dummy plot handles to be used for the legend
hLegend = zeros(3,1);
hLegend(1) = plot(NaN,NaN,'r');
hLegend(2) = plot(NaN,NaN,'b');
hLegend(3) = plot(NaN,NaN,'g');
grid on;
% Add the legend
legend(hLegend, {'40/0.005' , '35/0.0075', '30/ 0.01'}, 'Location', 'northwest');
hold off

subplot(322)
cla()  % Clear the current axes

a = squeeze(ResultsAh(1,:,:));  % Data for the first noise scenario
b = squeeze(ResultsAh(2,:,:));  % Data for the second noise scenario
c = squeeze(ResultsAh(3,:,:));  % Data for the third noise scenario

hold on
nConditions = size(a, 2);  % Number of conditions (should be 5, based on ResultsYh)

% Offset positions to avoid overlap: Create enough space between groups
positions_a = (1:nConditions) - 0.2;  % Slightly left
positions_b = (1:nConditions);        % Centered
positions_c = (1:nConditions) + 0.2;  % Slightly right

% Plot the boxplots with specified positions to avoid overlap
boxplot(a, 'Colors', 'r', 'positions', positions_a, 'widths', 0.1,'Notch','on','symbol','+');
boxplot(b, 'Colors', 'b', 'positions', positions_b, 'widths', 0.1,'Notch','on','symbol','+');
boxplot(c, 'Colors', 'g', 'positions', positions_c, 'widths', 0.1,'Notch','on','symbol','+');


% Customize x-tick labels to display 0 to 4
xticks(1:nConditions);         % Adjust the ticks to correspond to the middle set (b)
xticklabels(0:nConditions-1);  % X labels as 0 to 4

% Add title and labels
title('(b) \it E_a','FontSize', 11, 'FontWeight','normal');
%xlabel('Unkown end-members');
%ylabel('');

allData = [a(:); b(:); c(:)];  % Combine all data into one array
ylim([min(allData) max(allData)]);  % Set the y-axis limits based on data range

% Create dummy plot handles to be used for the legend
hLegend = zeros(3,1);
hLegend(1) = plot(NaN,NaN,'r');
hLegend(2) = plot(NaN,NaN,'b');
hLegend(3) = plot(NaN,NaN,'g');
grid on;
% Add the legend
legend(hLegend, {'40/0.005' , '35/0.0075', '30/ 0.01'}, 'Location', 'northwest');
hold off

subplot(323)
cla()  % Clear the current axes

a = squeeze(ResultsPh(1,:,:));  % Data for the first noise scenario
b = squeeze(ResultsPh(2,:,:));  % Data for the second noise scenario
c = squeeze(ResultsPh(3,:,:));  % Data for the third noise scenario

hold on
nConditions = size(a, 2);  % Number of conditions (should be 5, based on ResultsYh)

% Offset positions to avoid overlap: Create enough space between groups
positions_a = (1:nConditions) - 0.2;  % Slightly left
positions_b = (1:nConditions);        % Centered
positions_c = (1:nConditions) + 0.2;  % Slightly right

% Plot the boxplots with specified positions to avoid overlap
boxplot(a, 'Colors', 'r', 'positions', positions_a, 'widths', 0.1,'Notch','on','symbol','+');
boxplot(b, 'Colors', 'b', 'positions', positions_b, 'widths', 0.1,'Notch','on','symbol','+');
boxplot(c, 'Colors', 'g', 'positions', positions_c, 'widths', 0.1,'Notch','on','symbol','+');


% Customize x-tick labels to display 0 to 4
xticks(1:nConditions);         % Adjust the ticks to correspond to the middle set (b)
xticklabels(0:nConditions-1);  % X labels as 0 to 4

% Add title and labels
title('(c) \it E_p','FontSize',11, 'FontWeight','normal');
%xlabel('Unkown end-members');
%ylabel('');

allData = [a(:); b(:); c(:)];  % Combine all data into one array
ylim([min(allData) max(allData)]);  % Set the y-axis limits based on data range

% Create dummy plot handles to be used for the legend
hLegend = zeros(3,1);
hLegend(1) = plot(NaN,NaN,'r');
hLegend(2) = plot(NaN,NaN,'b');
hLegend(3) = plot(NaN,NaN,'g');
grid on;
% Add the legend
legend(hLegend, {'40/0.005' , '35/0.0075', '30/ 0.01'}, 'Location', 'northwest');
hold off

subplot(324)
cla()  % Clear the current axes

a = squeeze(ResultsPh2(1,:,:));  % Data for the first noise scenario
b = squeeze(ResultsPh2(2,:,:));  % Data for the second noise scenario
c = squeeze(ResultsPh2(3,:,:));  % Data for the third noise scenario

hold on
nConditions = size(a, 2);  % Number of conditions (should be 5, based on ResultsYh)

% Offset positions to avoid overlap: Create enough space between groups
positions_a = (1:nConditions) - 0.2;  % Slightly left
positions_b = (1:nConditions);        % Centered
positions_c = (1:nConditions) + 0.2;  % Slightly right

% Plot the boxplots with specified positions to avoid overlap
boxplot(a, 'Colors', 'r', 'positions', positions_a, 'widths', 0.1,'Notch','on','symbol','+');
boxplot(b, 'Colors', 'b', 'positions', positions_b, 'widths', 0.1,'Notch','on','symbol','+');
boxplot(c, 'Colors', 'g', 'positions', positions_c, 'widths', 0.1,'Notch','on','symbol','+');


% Customize x-tick labels to display 0 to 4
xticks(1:nConditions);         % Adjust the ticks to correspond to the middle set (b)
xticklabels(0:nConditions-1);  % X labels as 0 to 4

% Add title and labels
title('(d) \it E_{SAM}','FontSize',11, 'FontWeight','normal');

allData = [a(:); b(:); c(:)];  % Combine all data into one array
ylim([min(allData) max(allData)]);  % Set the y-axis limits based on data range

% Create dummy plot handles to be used for the legend
hLegend = zeros(3,1);
hLegend(1) = plot(NaN,NaN,'r');
hLegend(2) = plot(NaN,NaN,'b');
hLegend(3) = plot(NaN,NaN,'g');
grid on;

% Add the legend
legend(hLegend, {'40/0.005' , '35/0.0075', '30/ 0.01'}, 'Location', 'northwest');
hold off

subplot(325)
cla()  % Clear the current axes

a = squeeze(ResultsTh(1,:,:));  % Data for the first noise scenario
b = squeeze(ResultsTh(2,:,:));  % Data for the second noise scenario
c = squeeze(ResultsTh(3,:,:));  % Data for the third noise scenario

hold on
nConditions = size(a, 2);  % Number of conditions (should be 5, based on ResultsYh)

% Offset positions to avoid overlap: Create enough space between groups
positions_a = (1:nConditions) - 0.2;  % Slightly left
positions_b = (1:nConditions);        % Centered
positions_c = (1:nConditions) + 0.2;  % Slightly right

% Plot the boxplots with specified positions to avoid overlap
boxplot(a, 'Colors', 'r', 'positions', positions_a, 'widths', 0.1,'Notch','on','symbol','+');
boxplot(b, 'Colors', 'b', 'positions', positions_b, 'widths', 0.1,'Notch','on','symbol','+');
boxplot(c, 'Colors', 'g', 'positions', positions_c, 'widths', 0.1,'Notch','on','symbol','+');

% Customize x-tick labels to display 0 to 4
xticks(1:nConditions);         % Adjust the ticks to correspond to the middle set (b)
xticklabels(0:nConditions-1);  % X labels as 0 to 4

allData = [a(:); b(:); c(:)];  % Combine all data into one array
ylim([min(allData) max(allData)]);  % Set the y-axis limits based on data range

% Add title and labels
title('(e) Computational time (s)','FontSize',11, 'FontWeight','normal');
%ylabel('');

% Create dummy plot handles to be used for the legend
hLegend = zeros(3,1);
hLegend(1) = plot(NaN,NaN,'r');
hLegend(2) = plot(NaN,NaN,'b');
hLegend(3) = plot(NaN,NaN,'g');
grid on;
% Add the legend
legend(hLegend, {'40/0.005' , '35/0.0075', '30/ 0.01'}, 'Location', 'northwest');
hold off


annotation('textbox', [0.4, 0, 0.1, 0.1], 'string', 'Number of unknown end-members', 'FontSize', 11, 'FontWeight', 'normal', 'EdgeColor', 'none');
annotation('textbox', [0.1, 0, 0.1, 0.1], 'string', 'Number of unknown end-members', 'FontSize', 11, 'FontWeight', 'normal', 'EdgeColor', 'none');