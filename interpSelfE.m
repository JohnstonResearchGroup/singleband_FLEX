function [Zut Xut Put mu] = interpSelfE(WNin,varargin)
%INTERPSELFE Interpolate the self-energy on different Matsubara frequency
%grid WN = pi*(2*(-numwi:numwi-1)+1)/beta.

flySelfE = varargin{1};
if flySelfE == 1
    Z = varargin{2};
    X = varargin{3};
    P = varargin{4};
    WN = varargin{5};
    mu = varargin{6};
else
    dataFileName = varargin{2};
    load(dataFileName,'Z','X','P','WN','mu','-mat')
end
[nk1 nk2 nw] = size(Z);
nwin = numel(WNin);
if mod(nw,2)~=0 || mod(nwin,2)~=0
    error('Number of Matsubara frequencies must be even!')
end
if nw~=numel(WN)
    error('Wrong grid for self-energy data!')
end

Z = shiftdim(Z,2);          %now Z becomes [nw nk1 nk2]
X = shiftdim(X,2);
P = shiftdim(P,2);
Zut(1:nwin,1:nk1,1:nk2) = 0;
Xut(1:nwin,1:nk1,1:nk2) = 0;
Put(1:nwin,1:nk1,1:nk2) = 0;

%1D interpolation with spline-extrapolated values
Zut = interp1(WN,Z,WNin,'spline','extrap');
Xut = interp1(WN,X,WNin,'spline','extrap');
Put = interp1(WN,P,WNin,'spline','extrap');

%Replace spline-extrapolated values by 'nearest' values
idxL = find(WNin<WN(1));
idxR = find(WNin>WN(end));
if ~isempty(idxL)
    idxL = idxL(end);
    Zut(1:idxL,:,:) = Zut(repmat(idxL+1,1,idxL),:,:);
    Xut(1:idxL,:,:) = Xut(repmat(idxL+1,1,idxL),:,:);
    Put(1:idxL,:,:) = Put(repmat(idxL+1,1,idxL),:,:);
end
if ~isempty(idxR)
    idxR = idxR(1);
    Zut(idxR:nwin,:,:) = Zut(repmat(idxR-1,1,nwin-idxR+1),:,:);
    Xut(idxR:nwin,:,:) = Xut(repmat(idxR-1,1,nwin-idxR+1),:,:);
    Put(idxR:nwin,:,:) = Put(repmat(idxR-1,1,nwin-idxR+1),:,:);
end

Zut = shiftdim(Zut,1);      %now Zut becomes [nk1 nk2 nwin]
Xut = shiftdim(Xut,1);
Put = shiftdim(Put,1);