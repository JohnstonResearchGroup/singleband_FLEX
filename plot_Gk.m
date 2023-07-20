clear;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = 6;
n = 1.1;
T = 0.1;
data_path = ['../Data/Ueq' num2str(U) '_neq' num2str(n) '/'];
datafile1 = ['selfE_Nk=64_64_U=' num2str(U) '_T=' num2str(T) '.mat'];
datafile2 = ['selfE_Nk=64_64_T=' num2str(T) '_U=' num2str(U) '_freqFFT.mat'];
output_path = '../Analysis/';
Gfilename = ['Gwm=1_Ueq' num2str(U) '_neq' num2str(n) '_Teq' num2str(T) '.eps'];
load([data_path, datafile1]);
load([data_path, datafile2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the k-grid and bare band dispersion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nktot = 4*Nk(1)*Nk(2);
for ii = 1:2
    K{ii} = (0:1:2*Nk(ii)-1)*pi/Nk(ii);
end
[KX,KY] = meshgrid(K{2},K{1}); 

t = 1.0;                %hopping t as the energy units
tp = -0.2*t;            %hopping t'
% inline function defining the bare band.
xik = @(kx,ky) -2*t*(cos(kx)+cos(ky)) - 4*tp*cos(kx).*cos(ky); 
ek = xik(KX,KY);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the green's function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gn = zeros(size(Z));
for mm = 1:2*numwi
  den = -Z(:,:,mm).^2 - (ek(:,:) - mu + X(:,:,mm)).^2 - P(:,:,mm).^2;
  Gn(:,:,mm) = (1i*Z(:,:,mm) + ek(:,:) - mu + X(:,:,mm))./den(:,:);
  Fn(:,:,mm) = P(:,:,mm)./den(:,:);
end
gap = zeros(2*Nk); 
gap(:,:) = WN(numwi+1)*P(:,:,numwi+1)./Z(:,:,numwi+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Renderer', 'painters', 'Position', [10 10 600 900])
set(gcf,'color','white')

K = (0:1:Nk(1))/Nk(1);
subplot(2,1,1); hold on; imagesc(K,K,real(Gn(1:Nk(1),1:Nk(2),numwi+1)));
axis([0,1,0,1])
xlabel('$k_x~[\pi/a]$','FontSize',25,'Interpreter','latex');
ylabel('$k_y~[\pi/a]$','FontSize',25,'Interpreter','latex');
title('$\mathrm{Re}G(k,\pi/\beta)$','FontSize',25,'Interpreter','latex');
set(gca,'FontSize',25','XTick',[0:0.25:1],'YTick',[0:0.25:1])
contour(KX/pi,KY/pi,ek-mu,[0 0],'--w','linewidth',2)
colorbar

subplot(2,1,2); hold on; imagesc(K,K,imag(Gn(1:Nk(1),1:Nk(2),numwi+1)));
axis([0,1,0,1])
xlabel('$k_x~[\pi/a]$','FontSize',25,'Interpreter','latex');
ylabel('$k_y~[\pi/a]$','FontSize',25,'Interpreter','latex');
title('$\mathrm{Im}G(k,\pi/\beta)$','FontSize',25,'Interpreter','latex');
set(gca,'FontSize',25','XTick',[0:0.25:1],'YTick',[0:0.25:1])
contour(KX/pi,KY/pi,ek-mu,[0 0],'--w','linewidth',2)
colorbar

saveas(gcf,[output_path, Gfilename],'epsc');

