clc
clear
close all
% profile on

runInSerial = 1;    %switch between for-loop and parfor-loop
saveSelfE = 1;      %save self-energy
loadSelfE = 1;      %load self-energy at different T or U as the input
flySelfE = 0;       %load self-energy calculated on the fly
if (loadSelfE == 1) && (saveSelfE ~= 1)
    flySelfE = 1;
end

checkFreqSymm = 0;  %check frequency even/odd symmetry
%Nota Bene: Eliashberg equation for normal self-energy S = 1i*wn*(1-Z) + X is used.
%Although there is no explicit frequency even/odd symmetry, the even/odd
%symmetry for real(S)/imag(S) is imposed every iteration. Most error
%(<1e-14) related to the symmetry might come from the fft transformation,
%which depends on the version of fftw (C library) called by Matlab, hence
%the Matlab version as well. Another error related to the fft is the
%approximation of Fourier integral with fft in the fourier_bwd when
%transforming from \tau-space to \omega_n-space. This can be reduced by
%using a higher Newton-Cotes formula.

useSymmetry = 1;    %C_2(inversion) symmetry is used
NCorder = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the momentum grid and all of the constants that we need.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cpus = feature('numCores');
% pool = parpool(cpus)
Nk = [32 32];     %[nky nkx]

nktot = 4*Nk(1)*Nk(2);
for ii = 1:2
    K{ii} = (0:1:2*Nk(ii)-1)*pi/Nk(ii);
end
[KX,KY] = meshgrid(K{2},K{1});      %row: y, column: x
%--------------------------------------------------------------------------
kb = 1.0;
%kb = 8.617e-5;         %boltzmann
%kb = 8.6173303e-5;     %Source: 2014 CODATA
%--------------------------------------------------------------------------
t = 1.0;              %hopping t as the energy units
tp = -t/6;            %hopping t'
tpp = t/5;
xik = @(kx,ky) -2*t*(cos(kx)+cos(ky)) - 4*tp*cos(kx).*cos(ky) ...
               -2*tpp*(cos(2*kx)+cos(2*ky));  %dispersion
Ncut = 3; %5
nuscale = Ncut*(8*t)/2;
%Matsubara freq. cut-off w_c = Ncut*(band width 8t)
%w_c = (2*n_c-1)*pi*T --> n_c~(w_c/2)/(pi*T)~nuscale/wn

U = 5.0*t;
T = 0.01*t/kb;
beta = 1/(kb*T);%
% Tlist = [0.02:-0.001:0.001];
Tlist = 0.05;
Ulist = U;
%Ulist = [0.5:0.5:4];
%--------------------------------------------------------------------------
sig = 0.05*t;           %~v_F*\Delta k ~ 4*t*pi/nk
%mu0 = -0.60;            %128*128 k-grid, T=0.015t, U=0;
%mu0 = -0.65;            %3ls2*32   k-grid, T=0.1t,   U=4t;
mu0 = -2;
eps = 1e-5;
mixing = 1;
wgt = 0.3;
mix_method = 2;
%mix_method: 0, mix self-energy; 1, mix Green's function; 
%            2, mix self-energy -- Anderson acceleration

maxiter = 500;
%--------------------------------------------------------------------------
ek = xik(KX,KY);
% mu = mu0;
% fk = fermi(ek-mu,beta);
% filling0 = 2*sum(fk(:))/nktot;  %2 for two spins
filling0 = 1.15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup the file names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileDir = './';
filenamestr = ['_Nk=' num2str(Nk(1)) 'x' num2str(Nk(2)) '_U=' num2str(U) '.dat'];
fileout = ['out' filenamestr];
filegaps = ['gaps' filenamestr];
filechi = ['chi' filenamestr];
filelambda = ['lambda' filenamestr];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start outputting to the screen and file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen([fileDir,fileout],'w');
fclose(fid);
diary([fileDir, fileout]);
fprintf('\n');
fprintf('Gap calculation\n');
fprintf('  grid [nky nkx] = [%4d,%4d], convergence criterion = %g [t], smearing = %g [t]\n',...
    Nk(1),Nk(2),eps,sig)
fprintf('  U = %g [t]\n',U)
fprintf('  mixing = %d, weight = %g\n',mixing,wgt)
fprintf('  maxiter = %d\n',maxiter)
fprintf('  useSymmetry = %1d\n',useSymmetry)
fprintf('  NCorder = %1d\n',NCorder)
fprintf('  [loadSelfE, saveSelfE, flySelfE] = [%1d, %1d, %1d]\n',loadSelfE,saveSelfE,flySelfE)
fprintf('  runInSerial = %1d\n',runInSerial)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%tic;
fidgaps = fopen([fileDir filegaps],'a');
fidchi = fopen([fileDir filechi],'a'); 
fidlambda = fopen([fileDir filelambda],'a'); 
%Tlist = [0.4 0.3 0.2 0.1];

gap_input = 0.1*(cos(KX)-cos(KY));        %initial value for the gap
%gap_input = 0.001;
for nt = 1:numel(Tlist)
    tic;
    for iiu = 1:numel(Ulist)
        Told = T;
        Uold = U;
        T = Tlist(nt);
        U = Ulist(iiu);
        beta = 1/(kb*T);
        wn = pi/beta;
        numwi = round(nuscale/wn);

        fprintf('\n')
        fprintf('\n')
        fprintf('**************************************************************\n')
        fprintf('**************************************************************\n')
        fprintf('  T = %g [t], U = %g [t]\n',T,U)
        fprintf('  w(n=0) = %12.8f [t], number of wn>0 numwi =%6d\n',wn,numwi)
        if nt>1
            WNold = WN;
        end
        %define the freqency grid
        WN = pi*(2*(-numwi:numwi-1)+1)/beta;
        WNU = pi*2*(-numwi:numwi-1)/beta;
        if saveSelfE == 1
            fileSelfE{nt,iiu} = ['selfE_Nk=' num2str(Nk(1)) '_' num2str(Nk(2)) '_U=' num2str(U) '_T=' num2str(T) '.mat'];
            fprintf('  self-energy %s%s will be saved.\n',fileDir,fileSelfE{nt,iiu})
        end
        if (saveSelfE == 1) && exist([fileDir fileSelfE{nt,iiu}],'file')
            fprintf('  self-energy %s%s exists. Go to next T and U.\n',fileDir,fileSelfE{nt,iiu})
            continue
        end
        if (loadSelfE == 1) && (nt>1 || iiu>1)
            if numel(Tlist)>1 && numel(Ulist)>1 && saveSelfE~=1
                error('  set saveSelfE = 1!')
            end
            if numel(Tlist) == 1 && numel(Ulist) > 1
                if flySelfE ~= 1
                    load([fileDir fileSelfE{nt,iiu-1}],'Z','X','P','WN','mu','-mat')
                end
            elseif  numel(Tlist) > 1 && numel(Ulist) == 1
                if flySelfE == 1
                    [Z X P mu] = interpSelfE(WN,1,Z,X,P,WNold,mu);
                else
                    [Z X P mu] = interpSelfE(WN,0,[fileDir fileSelfE{nt-1,1}]);
                end
            else
                if (Told~=T)
                    [Z X P mu] = interpSelfE(WN,0,[fileDir fileSelfE{nt-1,iiu}]);
                elseif (Uold~=U)
                    load([fileDir fileSelfE{nt,iiu-1}],'Z','X','P','WN','mu','-mat')
                end
            end
            wn = pi/beta;
            gap = max(max( abs(real(wn*P(:,:,numwi+1)./Z(:,:,numwi+1))) ));
            if gap < 0.003 %0.001
                for nn = 1:length(WN)
                    P(:,:,nn) = gap_input;
                end
            end
        else            
            %setup the arrays to store the self-energies on the imaginary
            %axis and initialize the self-energies
            Z(1:2*Nk(1),1:2*Nk(2),1:2*numwi) = 0;
            X(1:2*Nk(1),1:2*Nk(2),1:2*numwi) = 0;
            P(1:2*Nk(1),1:2*Nk(2),1:2*numwi) = 0;
            S(1:2*Nk(1),1:2*Nk(2),1:2*numwi) = 0;
            for nn = 1:2*numwi
                Z(:,:,nn) = WN(nn);
                P(:,:,nn) = gap_input;
                S(:,:,nn) = 0;
            end
            %mu = get_mu(Z,X,P,WN,ek,Nk,beta,filling0);
            mu = mu0;
        end
        filling = get_filling(Z,X,P,ek,WN,Nk,beta,mu);
        fprintf('  filling0 = %12.8f\n',filling0)
        fprintf('  filling  = %12.8f, mu = %12.8f [t]\n',filling,mu)

        dk = gauss(ek-mu,sig);
        Nf0 = sum(dk(:))/nktot;
        fprintf('  N(E_f,T =%12.8f [t]) = %12.8f (1/[t]/per spin/Volume)\n',sig,Nf0)
        dfk = -beta*exp(beta*(ek-mu))./((exp(beta*(ek-mu))+1).^2);
        Nf = -sum(dfk(:))/nktot;
        fprintf('  N(E_f,T =%12.8f [t]) = %12.8f (1/[t]/per spin/Volume)\n',T,Nf)

        %%{
        
        calculate_im_axis;
        wn = pi/beta;

        %if (U == 0), minvsf = 1; end
        if (U == 0), maxchi = 0; end
        
        %if (minvsf >= 0)        
        if (maxchi < 1)
            gapwn = real(wn*P(:,:,numwi+1)./Z(:,:,numwi+1));
            gapmax = max(max(gapwn));
            gapmin = min(min(gapwn));
            fprintf(fidgaps,'%12.8f %12.8f %12.8f %12.8f %12.8f\n', ...
                T,min(min(ek(:,:)-mu)),mu,gapmax,gapmin);
        end
        
        %if (U==0), clear('minvsf'); end
        if (U==0), clear('maxchi'); end
        
        %}
        %fprintf(fidgaps,'\n');
        if saveSelfE == 1
            save([fileDir fileSelfE{nt,iiu}],'Z','X','P','WN','mu','-mat')
        end
        
        vtm = toc;
        
%         %Solve the BSE
%         for mm = 1:2*numwi
%           den = -Z(:,:,mm).^2 - (ek(:,:) - mu + X(:,:,mm)).^2 - P(:,:,mm).^2;
%           Gn(:,:,mm) = (1i*Z(:,:,mm) + ek(:,:) - mu + X(:,:,mm))./den(:,:);
%         end
%         Vnm = Vsf + Vcf;
%         Solve_BSE;

        
        % save(['T=' num2str(T) '_U=' num2str(U) '.mat'],'nk','numwi','T','mu','WN','Z','X','P',...
        %     'WNU','chis0','chic0','Vsf','Vcf','-mat')
        % q dependence
        chis00 = chis0(:,:,numwi+1);
        chic00 = chic0(:,:,numwi+1);
        Vsf0 = Vsf(:,:,numwi+1);
        Vcf0 = Vcf(:,:,numwi+1);
        Z0 = Z(:,:,numwi+1);
        X0 = X(:,:,numwi+1);
        P0 = P(:,:,numwi+1);
        
        % w dependence q = (0 0)
        chis0w1 = reshape(chis0(1,1,:),1,[]);
        chic0w1 = reshape(chic0(1,1,:),1,[]);
        Vsfw1 = reshape(Vsf(1,1,:),1,[]);
        Vcfw1 = reshape(Vcf(1,1,:),1,[]);
        Zw1 = reshape(Z(1,1,:),1,[]);
        Xw1 = reshape(X(1,1,:),1,[]);
        Pw1 = reshape(P(1,1,:),1,[]);
        % w dependence q = (pi pi)
        chis0w2 = reshape(chis0(Nk(1)+1,Nk(2)+1,:),1,[]);
        chic0w2 = reshape(chic0(Nk(1)+1,Nk(2)+1,:),1,[]);
        Vsfw2 = reshape(Vsf(Nk(1)+1,Nk(2)+1,:),1,[]);
        Vcfw2 = reshape(Vcf(Nk(1)+1,Nk(2)+1,:),1,[]);
        Zw2 = reshape(Z(Nk(1)+1,Nk(2)+1,:),1,[]);
        Xw2 = reshape(X(Nk(1)+1,Nk(2)+1,:),1,[]);
        Pw2 = reshape(P(Nk(1)+1,Nk(2)+1,:),1,[]);
        % w dependence q = (pi 0)
        chis0w3 = reshape(chis0(Nk(1)+1,1,:),1,[]);
        chic0w3 = reshape(chic0(Nk(1)+1,1,:),1,[]);
        Vsfw3 = reshape(Vsf(Nk(1)+1,1,:),1,[]);
        Vcfw3 = reshape(Vcf(Nk(1)+1,1,:),1,[]);
        Zw3 = reshape(Z(Nk(1)+1,1,:),1,[]);
        Xw3 = reshape(X(Nk(1)+1,1,:),1,[]);
        Pw3 = reshape(P(Nk(1)+1,1,:),1,[]);
        save([fileDir 'selfE_Nk=' num2str(Nk(1)) '_' num2str(Nk(2))  '_T=' num2str(T) '_U=' num2str(U) '_freqFFT.mat'],...
            'Nk','numwi','T','mu','WN','WNU',...
            'Z0','X0','P0','chis00','chic00','Vsf0','Vcf0',...
            'Zw*','Xw*','Pw*','chis0w*','chic0w*','Vsfw*','Vcfw*','-mat')
        % profile viewer
        % p = profile('info');
        % profsave(p,'profile_results')
        fprintf(fidchi,'%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n',...
              T, U, real(chis0(Nk(1)+1,Nk(2)+1,numwi+1)), real(chic0(Nk(1)+1,Nk(2)+1,numwi+1)), ...
                    real(Vsf(Nk(1)+1,Nk(2)+1,numwi+1)), real(Vcf(Nk(1)+1,Nk(2)+1,numwi+1)));
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('Done gap calculation. Total Time = %.2f s\n',sum(vtm));
        %save([fileDir fileout(1:end-3) 'mat'],'-mat')
        
   end
end
fclose(fidgaps);
fclose(fidchi);
diary off
return

