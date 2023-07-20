fprintf('\n')
fprintf('**************************************************************\n')
fprintf('Begin imaginary axis calculation.\n')

if (2*numwi) ~= numel(WN)
    error('Wrong Matsubara frequency grid!')
end
fmWt = getNCwt(beta,numwi,NCorder,1);   %fermion freq weight
bsWt = getNCwt(beta,numwi,NCorder,0);   %boson freq weight

pre = 1/(nktot*beta);
filling = get_filling(Z,X,P,ek,WN,Nk,beta,mu);
fprintf('  filling  = %12.8f, mu = %12.8f [t]\n',filling,mu)

%array for Fourier transforms
clear('Gn','Fn','chis0','chic0','Vsf','Vcf','Vnm','Van')
Gn(1:2*Nk(1),1:2*Nk(2),1:2*numwi) = 0;
Fn(1:2*Nk(1),1:2*Nk(2),1:2*numwi) = 0;
chis0(1:2*Nk(1),1:2*Nk(2),1:2*numwi) = 0;
chic0(1:2*Nk(1),1:2*Nk(2),1:2*numwi) = 0;
Vsf(1:2*Nk(1),1:2*Nk(2),1:2*numwi) = 0;
Vcf(1:2*Nk(1),1:2*Nk(2),1:2*numwi) = 0;
Vnm(1:2*Nk(1),1:2*Nk(2),1:2*numwi) = 0;
Van(1:2*Nk(1),1:2*Nk(2),1:2*numwi) = 0;
%array for end point
Gbeta(1:2*Nk(1),1:2*Nk(2)) = 0;


loop = 1;
iter = 0;
recomp_chi = 0; %2;
count = recomp_chi;
fprintf(['  Iter.  _______dZ_______  _______dX_______  _______dP_______',...
    '  _______mu_______  _______<n>______',...
    '  _______gap______\n'])
doneiter = 1;
delmu = 0.5;
maxchi_cutoff = 0.9995;
%Anderson acceleration
if (mixing == 1) && (mix_method == 2)
    nptmatZ = 2*Nk(1)*2*Nk(2)*2*numwi;
    xx(1:nptmatZ*2,1) = 0;  %store Zold, Xold, Pold in one column
    gg(1:2*Nk(1)*2*Nk(2)*2*numwi*2,1) = 0;  %store Z, X, P in one column
    %mMax = 5;       %mMax>0, 5
    %mMax = 2;
    mMax = 3;
    % maximum number of stored residuals (non-negative integer); should not
    % exceed number of rows of xx vectors (so we have a tall/thin matrix in
    % QR decompostion)
    droptol = 1.e10;
    % tolerance for dropping stored residual vectors to improve
    % conditioning: If droptol > 0, drop residuals if the condition number
    % exceeds droptol; if droptol <= 0, do not drop residuals.    
    %dampbt = @(x)0.2*(x<20) + 0.2*(x>=20 && x<25) + 0.5*(x>=25); %0.2
    %dampbt = @(x)0.5*(x<30) + 0.5*(x>=30);
    dampbt = 1;
    % damping factor: If dampbt > 0 (and dampbt ~= 1), then the step is damped
    % by dampbt; otherwise, the step is not damped. NOTE: dampbt can be a
    % function handle; form dampbt(iter), where iter is the iteration number
    % and 0 < dampbt(iter) <= 1.
    AAstart = 1;    %AAstart>=1, 5
    % acceleration delay factor: If AAstart > 0, start acceleration when
    % iter = AAstart.
    res_hist = [];
    % residual history matrix (iteration numbers and residual norms).
    DG = [];
    % Storage of g-value differences.
    R = []; Q = [];
    mAA = 0;
    % Initialize the number of stored residuals.
end

while (loop == 1)
    iter = iter + 1;

    %updata old self-energies (with mixing)
    if (doneiter == 1)
        if (mixing == 1) && (iter>1)
            if (mix_method == 0)
                Zold = wgt*Z + (1-wgt)*Zold;
                Xold = wgt*X + (1-wgt)*Xold;
                Pold = wgt*P + (1-wgt)*Pold;
            elseif mix_method == 1
                %mix Green's function Gn and Fn
                %(a) update Gn/Fn with pre-mixed self-energies
                for mm = 1:2*numwi
                    den = -Z(:,:,mm).^2 - (ek(:,:) - mu + X(:,:,mm)).^2 - P(:,:,mm).^2;
                    Gn(:,:,mm) = (1i*Z(:,:,mm) + ek(:,:) - mu + X(:,:,mm))./den(:,:);
                    Fn(:,:,mm) = P(:,:,mm)./den(:,:);
                end
                %(b) mix the Green's function
                Gnold = wgt*Gn + (1-wgt)*Gnold;
                Fnold = wgt*Fn + (1-wgt)*Fnold;
                %(c) find the self-energies from the mixed Green's function [not necessary]
                Zold = imag(Gnold);
                Xold = real(Gnold);
                Pold = real(Fnold);
                denom = -1./(Zold.^2 + Xold.^2 + Fnold.^2);
                Zold = Zold.*denom;
                Xold = Xold.*denom - repmat(ek(:,:)-mu, [1, 1, 2*numwi]);
                Pold = Pold.*denom;                
            else
                if mix_method ~= 2
                    error('mix_method = %d is not supported!',mix_method)
                end
            end
        else
            Zold = Z;
            Xold = X;
            Pold = P;
        end
    end

    %update Green's function G_{1,1} = Gn, G_{1,2} = Fn; updata old normal
    %self-energy: S (old anomalous self-energy: P)
    Fn(:) = 0;      %set Fn real, deallocate imaginary part at once (old version Matlab problem)
    %parfor
    if runInSerial == 1
        for mm = 1:2*numwi
            den = -Zold(:,:,mm).^2 - (ek(:,:) - mu + Xold(:,:,mm)).^2 - Pold(:,:,mm).^2;
            Gn(:,:,mm) = (1i*Zold(:,:,mm) + ek(:,:) - mu + Xold(:,:,mm))./den(:,:);
            Fn(:,:,mm) = Pold(:,:,mm)./den(:,:);
            %Sold(:,:,mm) = 1i*(repmat(WN(mm),Nk(1),Nk(2)) - Zold(:,:,mm)) + Xold(:,:,mm);
        end
        if (mix_method == 1) && (mixing == 1) && (iter>1)
            tmp = abs(Gn - Gnold);
            err1 = max(tmp(:));
            tmp = abs(Fn - Fnold);
            err2 = max(tmp(:));
            if max([err1,err2]) > 1e-12
                fprintf('    %g, %g\n',err1,err2)
            end
        end
    else
        parfor mm = 1:2*numwi
            den = -Zold(:,:,mm).^2 - (ek(:,:) - mu + Xold(:,:,mm)).^2 - Pold(:,:,mm).^2;
            Gn(:,:,mm) = (1i*Zold(:,:,mm) + ek(:,:) - mu + Xold(:,:,mm))./den(:,:);
            Fn(:,:,mm) = Pold(:,:,mm)./den(:,:);
            %Sold(:,:,mm) = 1i*(repmat(WN(mm),Nk(1),Nk(2)) - Zold(:,:,mm)) + Xold(:,:,mm);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %save old Green's function Gn(k,i\omega_n) and Fn(k,i\omega_n) for mixing
    if (mix_method == 1 || mix_method == 2) && (mixing == 1) && (iter==1)
        Gnold = Gn;
        Fnold = Fn;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %update Hartree-Fock self-energies
    %It should not be folded into FLEX interaction because U is
    %instantaneous so the Fourier transform to \tau-space will be
    %inaccurate. (FLEX interaction decays faster than 1/w_n while U is
    %constant for any w_n.)
    %S_HF = U*sum(Gn(:))*pre;
    S_HF = 0;   %temporarily, we set this to 0 to compare with Monthoux and Scalapino
    P_HF = U*sum(Fn(:))*pre;
   
    %compute the self-energy ----------------------------------------------
    %forward Fourier transform Green's function (k,i\omega_n) to (r,\tau)
    Gn = fourier_fwd(Gn,ek-mu,beta,Nk,numwi,1,1,useSymmetry,runInSerial);
    Fn = fourier_fwd(Fn,ek-mu,beta,Nk,numwi,0,1,useSymmetry,runInSerial);      
    Gbeta = -Gn(:,:,1);
    Gbeta(1,1) = Gbeta(1,1) - 1;
    chis0(:,:,1) = Gn(:,:,1).*Gbeta - Fn(:,:,1).*Fn(:,:,1);
    chic0(:,:,1) = Gn(:,:,1).*Gbeta + Fn(:,:,1).*Fn(:,:,1);
    %parfor
    if runInSerial == 1
        for mm = 2:2*numwi
            chis0(:,:,mm) = Gn(:,:,mm).*Gn(:,:,2*numwi+2-mm) - Fn(:,:,mm).*Fn(:,:,mm);
            chic0(:,:,mm) = Gn(:,:,mm).*Gn(:,:,2*numwi+2-mm) + Fn(:,:,mm).*Fn(:,:,mm);
        end
    else
        parfor mm = 2:2*numwi
            chis0(:,:,mm) = Gn(:,:,mm).*Gn(:,:,2*numwi+2-mm) - Fn(:,:,mm).*Fn(:,:,mm);
            chic0(:,:,mm) = Gn(:,:,mm).*Gn(:,:,2*numwi+2-mm) + Fn(:,:,mm).*Fn(:,:,mm);
        end
    end
    %backward Fourier transform susceptibility (r,\tau) to (k,i\omega_n)
    chis0 = fourier_bwd(chis0,0,NCorder,bsWt,beta,Nk,numwi,0,useSymmetry,runInSerial);
    chic0 = fourier_bwd(chic0,0,NCorder,bsWt,beta,Nk,numwi,0,useSymmetry,runInSerial);
    if (U ~= 0)
        %{
        minvsf = min(min(real(Vsf(:,:,numwi+1))));
        if (minvsf < 0)
            %loop = 0;
            %fprintf('  Warning: min Vsf = %g\n',minvsf)
            %Vsf(Vsf<0) = 0;
            surf(KX,KY,real(Vsf(:,:,numwi+1)))
            shading flat
            error('  Warning: min Vsf = %g\n',minvsf)
            %break
        end
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        maxchi = U*max(max(real(chis0(:,:,numwi+1))));
        if (maxchi >= maxchi_cutoff)
            fprintf('  Warning: max U*chis0 = %g, decrease mu = %g --> %g\n',maxchi,mu,mu-delmu)
            mu = mu - delmu;
            iter = iter - 1;
            doneiter = 0;
            continue
        else
            doneiter = 1;
        end
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%{
        maxchi = U*max(max(real(chis0(:,:,numwi+1))));
        if (maxchi > maxchi_cutoff)
            fprintf('  Warning: max U*chis0 = %g,  decrease U*chis0 to U*chi_max = %g\n',maxchi,maxchi_cutoff)
            chis0(real(chis0) > maxchi_cutoff/U) = maxchi_cutoff/U; 
            
            %%fprintf('%d\n',isempty(find(real(chis0) > maxchi_cutoff)))
            %%chis0(real(chis0) > maxchi_cutoff) = maxchi_cutoff;   %bug code: 1/U should be included
            
            %fprintf('  Warning: max U*chis0 = %g\n',maxchi)
        end
        %}        
    end    
    %compute the full spin/charge-fluctuation propagator (RPA interactions)
    Vsf = 1.5*U*U*chis0./(1 - U*chis0) - 0.5*U*U*chis0;
    Vcf = 0.5*U*U*chic0./(1 + U*chic0) - 0.5*U*U*chic0;
    Vnm = Vsf + Vcf;
    Van = Vsf - Vcf;
    %forward Fourier transform effective interaction (k,i\omega_n) to (r,\tau)
    Vnm = fourier_fwd(Vnm,ek-mu,beta,Nk,numwi,0,0,useSymmetry,runInSerial);
    Van = fourier_fwd(Van,ek-mu,beta,Nk,numwi,0,0,useSymmetry,runInSerial);
    S = Vnm.*Gn;
    P = Van.*Fn;
    %backward Fourier transform self-energies (r,\tau) to (k,i\omega_n)
    S = fourier_bwd(S,-Vnm(1,1,1),NCorder,fmWt,beta,Nk,numwi,1,useSymmetry,runInSerial);
    P = fourier_bwd(P,0,NCorder,fmWt,beta,Nk,numwi,1,useSymmetry,runInSerial);
    S = S + S_HF;
    P = P + P_HF;
    %calculate Z,X; symmetrize +/- frequency; even-k symmetry used (singlet pairing)
    %parfor
    if runInSerial == 1
        for mm = (numwi+1):(2*numwi)
            Z(:,:,mm) = WN(mm);
            %Z(:,:,mm) = Z(:,:,mm) - imag(S(:,:,mm)-S(:,:,2*numwi+1-mm))/2;
            %X(:,:,mm) = real(S(:,:,mm)+S(:,:,2*numwi+1-mm))/2;
            Z(:,:,mm) = Z(:,:,mm) - imag(S(:,:,mm));
            X(:,:,mm) = real(S(:,:,mm));
        end
    else
        parfor mm = (numwi+1):(2*numwi)
            %Z(:,:,mm) = WN(mm);        %this is not allowed in parfor
            Z(:,:,mm) = repmat(WN(mm),2*Nk(1),2*Nk(2));
            %Z(:,:,mm) = Z(:,:,mm) - imag(S(:,:,mm)-S(:,:,2*numwi+1-mm))/2;
            %X(:,:,mm) = real(S(:,:,mm)+S(:,:,2*numwi+1-mm))/2;
            Z(:,:,mm) = Z(:,:,mm) - imag(S(:,:,mm));
            X(:,:,mm) = real(S(:,:,mm));
        end
    end
    P = real(P);
    %----------------------------------------------------------------------
    %impose d_{x^2-y^2}-wave symmetry
    %{
    for mm = (numwi+1):(2*numwi)
        Ptmp = triu(P(1:(end/2+1),1:(end/2+1),mm));
        Ptmp = Ptmp - transpose(Ptmp);        
        P(:,:,mm) = Ptmp([1:end, end-1:-1:2],[1:end, end-1:-1:2]);
    end
    %}
    %----------------------------------------------------------------------    
    for mm = (numwi+1):(2*numwi)
        Z(:,:,2*numwi+1-mm) =-Z(:,:,mm);
        X(:,:,2*numwi+1-mm) = X(:,:,mm);
        P(:,:,2*numwi+1-mm) = P(:,:,mm);
    end
    if checkFreqSymm == 1
        tmp = imag(S);
        errZ = min(min(min(abs(tmp(:,:,1:end) + tmp(:,:,end:-1:1)))));
        tmp = real(S);
        errX = min(min(min(abs(tmp(:,:,1:end) - tmp(:,:,end:-1:1)))));
        tmp = imag(P);
        errP = min(min(min(abs(P))));
        if errZ>1e-16 || errX>1e-16 || errP>1e-16
            error('  Error of frequency symmetry is too large! [errZ errX errP] =[%g %g %g]',...
                errZ,errX,errP)
        end
    end
    %----------------------------------------------------------------------

    % Need to insert code to compute the new filling ect
    mu = get_mu(Z,X,P,WN,ek,Nk,beta,filling0);
    filling = get_filling(Z,X,P,ek,WN,Nk,beta,mu);

    wn = pi/beta;
    gap = max(max( abs(real(wn*P(:,:,numwi+1)./Z(:,:,numwi+1))) ));

    % Decide if we have to exit
    diffz = 0;
    diffp = 0;
    diffx = 0;
    for nn = 1:numwi
        diffz = max([diffz, max(max(abs(Z(:,:,nn)-Zold(:,:,nn))))]);
        diffx = max([diffx, max(max(abs(X(:,:,nn)-Xold(:,:,nn))))]);
        diffp = max([diffp, max(max(abs(P(:,:,nn)-Pold(:,:,nn))))]);
    end
    if (diffz < eps && diffx < eps && diffp < eps), loop = 0; end
    fprintf('  %5d  %16.12f  %16.12f  %16.12f',iter,diffz,diffx,diffp)
    fprintf('  %16.12f  %16.12f',mu,filling)
    fprintf('  %16.12f\n',gap)
    %Decide if we need to exit the loop.
    if (iter >= maxiter)
        loop = 0;
        fprintf('  Maximal iterations =%5d reached. Check convergence!\n',maxiter)
        %continue
    end
    %Anderson acceleration
    if (loop ~= 0) && (mixing == 1) && (mix_method == 2)
        % Compute the current residual norm.
        %xx = real([Zold(:); Xold(:); Pold(:)]);
        %gval = real([Z(:); X(:); P(:)]);

        for mm = 1:2*numwi
            den = -Z(:,:,mm).^2 - (ek(:,:) - mu + X(:,:,mm)).^2 - P(:,:,mm).^2;
            Gn(:,:,mm) = (1i*Z(:,:,mm) + ek(:,:) - mu + X(:,:,mm))./den(:,:);
            Fn(:,:,mm) = P(:,:,mm)./den(:,:);
        end
        xx = [Gnold(:); Fnold(:)];
        gval = [Gn(:); Fn(:)];       
                
        fval = gval - xx;
        res_norm = max(abs(fval));   % = norm(A,inf)
        %fprintf('     -- %d %g \n', iter, res_norm);
        res_hist = [res_hist;[iter,res_norm]];
        
        if iter < AAstart %|| mMax == 0 %note we set mMax>0
            % Without acceleration, update x <- g(x) to obtain the next
            % approximate solution.
            
            %Zold = wgt*Z + (1-wgt)*Zold;
            %Xold = wgt*X + (1-wgt)*Xold;
            %Pold = wgt*P + (1-wgt)*Pold;
            
            Gnold = wgt*Gn + (1-wgt)*Gnold;
            Fnold = wgt*Fn + (1-wgt)*Fnold;
            %find the self-energies from the mixed Green's function [not necessary]
            Zold = imag(Gnold);
            Xold = real(Gnold);
            Pold = real(Fnold);
            denom = -1./(Zold.^2 + Xold.^2 + Fnold.^2);
            Zold = Zold.*denom;
            Xold = Xold.*denom - repmat(ek(:,:)-mu, [1, 1, 2*numwi]);
            Pold = Pold.*denom;
            
            %xx = gval;
        else
            % Apply Anderson acceleration.
            
            % Update the df vector and the DG array.
            if iter > AAstart
                df = fval-f_old;
                if mAA < mMax
                    DG = [DG gval-g_old];
                else
                    DG = [DG(:,2:mAA) gval-g_old];
                end
                mAA = mAA + 1;
            end
            f_old = fval;
            g_old = gval;

            if mAA == 0     %iter == AAstart
                % If mAA == 0, update x <- g(x) to obtain the next approximate solution.
                fprintf('  %d, start Anderson acceleration and store %d steps of self-energies\n', iter, mMax);
                
                %Zold = Z;
                %Xold = X;
                %Pold = P;
                
                Gnold = Gn;
                Fnold = Fn;
                %find the self-energies from the mixed Green's function [not necessary]
                Zold = imag(Gnold);
                Xold = real(Gnold);
                Pold = real(Fnold);
                denom = -1./(Zold.^2 + Xold.^2 + Fnold.^2);
                Zold = Zold.*denom;
                Xold = Xold.*denom - repmat(ek(:,:)-mu, [1, 1, 2*numwi]);
                Pold = Pold.*denom;
                
                %xx = gval;
            else
                % If mAA > 0, solve the least-squares problem and update the solution.
                if mAA == 1
                    % If mAA == 1, form the initial QR decomposition.
                    R(1,1) = norm(df);
                    Q = R(1,1)\df;
                else
                    % If mAA > 1, update the QR decomposition.
                    if mAA > mMax
                        % If the column dimension of Q is mMax, delete the first column and
                        % update the decomposition.
                        [Q,R] = qrdelete(Q,R,1);
                        mAA = mAA - 1;
                        if size(R,1) ~= size(R,2)
                            Q = Q(:,1:mAA-1); R = R(1:mAA-1,:);
                        end
                    end
                    % Now update the QR decomposition to incorporate the new column.
                    for j = 1:mAA - 1
                        R(j,mAA) = Q(:,j)'*df;
                        df = df - R(j,mAA)*Q(:,j);
                    end
                    R(mAA,mAA) = norm(df);
                    Q = [Q, R(mAA,mAA)\df];
                end
                if droptol > 0
                    % Drop residuals to improve conditioning if necessary.
                    condDF = cond(R);
                    while condDF > droptol && mAA > 1
                        fprintf(' cond(D) = %e, reducing mAA to %d \n', condDF, mAA-1);
                        [Q,R] = qrdelete(Q,R,1);
                        DG = DG(:,2:mAA);
                        mAA = mAA - 1;
                        % The following treats the qrdelete quirk described above.
                        if size(R,1) ~= size(R,2)
                            Q = Q(:,1:mAA); R = R(1:mAA,:);
                        end
                        condDF = cond(R);
                    end
                end
                % Solve the least-squares problem.
                gamma = R\(Q'*fval);
                % Update the approximate solution.
                xx = gval - DG*gamma;
                % Apply damping if dampbt is a function handle or if dampbt > 0
                % (and dampbt ~= 1).
                if isa(dampbt,'function_handle')
                    xx = xx - (1-dampbt(iter))*(fval - Q*R*gamma);
                else
                    if dampbt > 0 && dampbt ~= 1
                        xx = xx - (1-dampbt)*(fval - Q*R*gamma);
                    end
                end
                
                %Zold(:) = xx(1:nptmatZ);
                %Xold(:) = xx(1+nptmatZ:nptmatZ*2);
                %Pold(:) = xx(1+nptmatZ*2:nptmatZ*3);
                
                Gnold(:) = xx(1:nptmatZ);
                Fnold(:) = xx(1+nptmatZ:nptmatZ*2);
                %find the self-energies from the mixed Green's function [not necessary]
                Zold = imag(Gnold);
                Xold = real(Gnold);
                Pold = real(Fnold);
                denom = -1./(Zold.^2 + Xold.^2 + Fnold.^2);
                Zold = Zold.*denom;
                Xold = Xold.*denom - repmat(ek(:,:)-mu, [1, 1, 2*numwi]);
                Pold = Pold.*denom;                
            end
        end %if iter < AAstart
    end %if (mixing == 1) && (mix_method == 2)
end %self-consistency loop
