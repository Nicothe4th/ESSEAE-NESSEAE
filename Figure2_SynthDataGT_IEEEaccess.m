%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Display the Ground-Truth Synthetic Dataset
% 
%
% October/2024
% JNMC-DUCD-UASLP
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=4;                % Number of End-members
Nsamples=64;       % Size of the Squared Image Nsamples x Nsamples 
SNR=40;             % Additive Gaussian noise
density=0.005;       % Density of sparse noise component
ModelType=5;        % 0 --> Linear Mixing Model and 5 --> Multilinear Model
EndMembersSynth=1;  % 0--> USGS Spectral Library Version & 1 --> Craneotomy 


[Z,P0,A0,V0,D0]=MatternGaussian_Sparse_Synth(SNR,density,ModelType);


nRow=Nsamples;
nCol=Nsamples;
x=linspace(450,700,281);
%%
figure(1)
subplot(2,1,1);
plot(x,P0,'LineWidth',2);
set(gca, 'YScale', 'log');
xlabel('Wavelength (nm)','FontSize', 10); axis tight; 
title('(a) End-members','FontSize', 12); grid on;
ylabel('Normalized Intensity','FontSize', 12);
lgd = legend('#1','#2','#3','#4'); % Create the legend
set(lgd, 'FontSize', 8);         % Set the font size

% Second row title
annotation('textbox', [0.15, 0.45, 0.7, 0.05], 'string', '(b) Abundance maps', ...
    'FontSize', 12, 'FontWeight', 'Bold', 'HorizontalAlignment', 'center', 'EdgeColor', 'none');

% Second row: Abundance maps with color bar
subplot(2,4,5);
imshow(reshape(A0(1,:), nRow, nCol), [0,1]);
title('End-member 1', 'FontSize', 10);
xlabel('pixels');
ylabel('pixels');
axis on;

subplot(2,4,6);
imshow(reshape(A0(2,:), nRow, nCol), [0,1]);
title('End-member 2', 'FontSize', 10);
xlabel('pixels');
ylabel('pixels');
axis on;
subplot(2,4,7);
imshow(reshape(A0(3,:), nRow, nCol), [0,1]);
title('End-member 3', 'FontSize', 10);
xlabel('pixels');
ylabel('pixels');
axis on;
subplot(2,4,8);
imshow(reshape(A0(4,:), nRow, nCol), [0,1]);
title('End-member 4', 'FontSize', 10);
xlabel('pixels');
ylabel('pixels');
axis on;
% Add color bar next to the last image
colormap parula;
cb = colorbar;
cb.Position = [0.92, 0.15, 0.02, 0.26];  % Adjust position to align with the second row
cb.FontSize = 10;
