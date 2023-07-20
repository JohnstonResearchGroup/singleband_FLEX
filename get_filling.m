function out = get_filling(Z,X,P,ek,WN,Nk,beta,mu)

nktot = 4*Nk(1)*Nk(2);
numwi = length(WN)/2;

fill = 0;       %for two spins
fill = fill + 2*sum(sum(fermi(ek-mu,beta)))/nktot;

%slower -------------------------------------------------------------------
% green = zeros(2*Nk(1),2*Nk(2));
% den = zeros(2*Nk(1),2*Nk(2));
% for nn = (numwi+1):(numwi*2)
%     wn = WN(nn);
%     den(:,:) = Z(:,:,nn).^2 + (ek(:,:)-mu+X(:,:,nn)).^2 + P(:,:,nn).^2;
%     green(:,:) = (ek(:,:)-mu+X(:,:,nn))./den;
%     fill = fill - 4*sum(sum(green))/(nktot*beta);       %4=2(from +/-wn) * 2(from trace<-->spin)
% 
%     den(:,:) =  wn^2 + (ek(:,:)-mu).^2;
%     green(:,:) = (ek(:,:)-mu)./den;
%     fill = fill + 4*sum(sum(green))/(nktot*beta);
% end

%faster -------------------------------------------------------------------
green(1:nktot,1) = 0;
den(1:nktot,1) = 0;
Z = reshape(Z,nktot,[]);
X = reshape(X,nktot,[]);
P = reshape(P,nktot,[]);
ek = reshape(ek-mu,nktot,[]);
for nn = (numwi+1):(numwi*2)
    wn = WN(nn);
    den = Z(:,nn).^2 + (ek+X(:,nn)).^2 + P(:,nn).^2;
    green = (ek+X(:,nn))./den;
    fill = fill - 4*sum(green)/(nktot*beta);       %4=2(from +/-wn) * 2(from trace<-->spin)

    den =  wn^2 + ek.^2;
    green = ek./den;
    fill = fill + 4*sum(green)/(nktot*beta);
end

out = fill;
