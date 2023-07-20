disp('Solving the BSE.')

nk = 2*Nk(1);
N = nk*nk*2*numwi;
mat = zeros(N,N);
pre = 1/(nk*nk*beta);

%Loop over the k elements of the matrix
for nkx = 0:nk-1
  for nky = 0:nk-1
    for nwk = 1:length(WN)
      %Index for the k' variable
      idxk = (nwk-1)*nk*nk + (nky)*nk + nkx + 1;
      
      %Loop over the k elements of the matrix
      for npx = 0:nk-1
        for npy = 0:nk-1
          for nwp = 1:length(WN)
            idxp = (nwp-1)*nk*nk + (npy)*nk + npx + 1;
            
            %Determine the momentum transfer
            nqx = mod(nkx - npx,nk);
            nqy = mod(nky - npy,nk);
            nwq = find(WNU == WN(nwk) - WN(nwp));
            %Index for the k variable
            
            if(~isempty(nwq))
              G2 = Gn(npx+1,npy+1,nwp)*conj(Gn(npx+1,npy+1,nwp));
              mat(idxk,idxp) = mat(idxk,idxp) - pre*G2*Vnm(nqx+1,nqy+1,nwq);
            end
            
          end
        end
      end
    end
  end
end

%Now diagonalize to get the leading eigenvector
[eigenvec,eigenval] = eigs(mat,3,'la');

fprintf(fidlambda,'%12.8f %12.8f %12.8f %12.8f %12.8f \n',...
              T, U, eigenval(1,1), eigenval(2,2), eigenval(3,3));

phi1 = reshape(eigenvec(:,1),[nk,nk,2*numwi]);
phi2 = reshape(eigenvec(:,2),[nk,nk,2*numwi]);
phi3 = reshape(eigenvec(:,3),[nk,nk,2*numwi]);

disp(['lambda = ' num2str(eigenval(1,1))]);

K = (0:1:nk-1)*2*pi/nk;
close all; figure(2); 
subplot(2,2,1); hold on;
imagesc(K/pi,K/pi,real( phi1(:,:,numwi+1)))
title('$\phi_1({\bf k},\pi/\beta)$','Interpreter','latex');
axis([0,2,0,2])

subplot(2,2,2); hold on;
imagesc(K/pi,K/pi,real( phi2(:,:,numwi+1)))
title('$\phi_2({\bf k},\pi/\beta)$','Interpreter','latex');
axis([0,2,0,2])

subplot(2,2,3); hold on;
imagesc(K/pi,K/pi,real( phi3(:,:,numwi+1)))
title('$\phi_3({\bf k},\pi/\beta)$','Interpreter','latex');
axis([0,2,0,2])

subplot(2,2,4); hold on;
imagesc(K/pi,K/pi,(cos(KX)-cos(KY))/2)
title('$\cos(k_xa)-\cos(k_ya)$','Interpreter','latex');
axis([0,2,0,2])

fidvec = fopen(['eigenvector_T=' num2str(T) filenamestr],'w'); 
for nkx = 1:nk
  for nky = 1:nk
    fprintf(fidvec,'%12.8f %12.8f %12.8f \n',...
              KX(nkx,nky), KY(nkx,nky), phi1(nkx,nky)); 
  end
end
fclose(fidvec);