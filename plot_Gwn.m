clear;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 1.1;
T = 0.1;
output_path = '../Analysis/';
%output_file = ['G_k=X_n=' num2str(n) '_T=' num2str(T) '.eps'];
%output_file = ['G_k=M_n=' num2str(n) '_T=' num2str(T) '.eps'];
output_file = ['G_k=piby2_piby2_n=' num2str(n) '_T=' num2str(T) '.eps'];

t = 1.0;                %hopping t as the energy units
tp = -0.2*t;            %hopping t'

% inline function defining the bare band.
xik = @(kx,ky) -2*t*(cos(kx)+cos(ky)) - 4*tp*cos(kx).*cos(ky); 

figure('Renderer', 'painters', 'Position', [10 10 500 500])
set(gcf,'color','white')

Ulist = [4:1:6];
linestr = ['-ks';'-bv';'-ro'];
facecolorstr = ['k','b','r'];

miny = 1000;
maxy = 0;
for i = 1:numel(Ulist)
  U = Ulist(i);
  data_path = ['../Data/Ueq' num2str(U) '_neq' num2str(n) '/'];
  datafile1 = ['selfE_Nk=64_64_U=' num2str(U) '_T=' num2str(T) '.mat'];
  datafile2 = ['selfE_Nk=64_64_T=' num2str(T) '_U=' num2str(U) '_freqFFT.mat'];
  load([data_path, datafile1]);
  load([data_path, datafile2]);

  nktot = 4*Nk(1)*Nk(2);
  for ii = 1:2
    K{ii} = (0:1:2*Nk(ii)-1)*pi/Nk(ii);
  end
  [KX,KY] = meshgrid(K{2},K{1}); 
  ek = xik(KX,KY);
  
  Gn = zeros(size(Z));
  SE = zeros(size(WN));
  
  %idx = 1; idy = Nk(1)+1; kstr = 'X'; %X-point
  %idx = Nk(1)+1; idy = Nk(2)+1; kstr = 'M'; %X-point  
  idx = Nk(1)/2+1; idy = Nk(2)/2+1; kstr = '(pi/2,\pi/2)'; %X-point  

  for mm = 1:2*numwi
    SE(mm) = WN(mm) - Z(idx,idy,mm) + X(idx,idy,mm);
  end
  miny = min([miny,SE]);
  maxy = max([maxy,SE]);
  LegendStr = ['U = ' num2str(Ulist(i)) 't'];
  hold on; 
  plot(WN(:),real(SE), linestr(i,:), 'DisplayName',LegendStr, ...
       'MarkerFaceColor', facecolorstr(i));
end

xlabel('$\mathrm{i}\omega_n~[t]$','FontSize',25,'Interpreter','latex');
ylabel(['$\Sigma({\bf k} = ' kstr ',\mathrm{i}\omega_n)~[t]$'],'FontSize',25,'Interpreter','latex');
set(gca,'FontSize',25','XTick',[-20:5:20],'YTick',[-10:1:10])
axis([-20,20,miny*1.1,maxy*1.1])
box on;
legend('location','northeast');
legend boxoff;
saveas(gcf,[output_path,output_file],'epsc');