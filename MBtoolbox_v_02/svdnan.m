function [U,D,V,Xrec] = svdnan(X,nF,options)
% function [U,D,V,Xrec] = svdnan(X,nF,options)
% 040706 FvdB
% Computes a Singular Value Decomposition for data tables with missing values via imputation.
%
% in:
% X (objects x variables) data table containing missing values as NaN
% nF (1 x 1) number of singular values (factors) used for reconstructing "X"
%    (CAUTION: this is an important decision for the quality of reconstruction, should be equal to the (assumed) rank!)
% options (1 x 3) centering, tolerance for convergence and maximum number of iterations (default 0=no, 1e-8 and 2000)
%    centering means center during computation, not compute mean centered solution
%
% out:
% U (objects x nF) left singular vectors
% D (nF x nF) singular values
% V (variables x nF) right singular vectors
% Xrec (objects x variables) data table "X" reconstructed via imputation
%   if options(1) = 1 (centering) then U, D and V are based on centered "X"
%
% based on:
% B.Walczak and D.L.Massart 'Dealing with missing data Part I' Chemo. Lab. 58(2001)15-27

if (nargin < 1)
    help svdnan
    return
end
[nX,mX] = size(X);
if nargin == 1
    options = [0 1e-8 2000];
    nF = min([nX mX]);
elseif nargin == 2
    options = [0 1e-8 2000];
else
    if length(options) == 1
        options = [options 1e-8 2000];
    elseif length(options) == 2
        options = [options 2000];
    end
    if isempty(nF)
        nF = min([nX mX]);
    end    
end

MV = find(isnan(X));
iter = 1;
if isempty(MV)
    if options(1)
        [X,mx] = meanc(X);
    end
    [U,D,V] = svd(X,0);
    U = U(:,1:nF);
    D = D(1:nF,1:nF);
    V = V(:,1:nF);
    Xrec = U*D*V';
    if options(1)
        Xrec = meanc(Xrec,mx,1);
    end
else
    Xnan = isnan(X);
    [MVi,MVj] = find(Xnan);
    X(MV) = 0;
    colmean = sum(X)./(ones(1,mX)*nX-sum(Xnan));
    rowmean = sum(X')./(ones(1,nX)*mX-sum(Xnan'));
    for a=1:length(MVi)
        X(MVi(a),MVj(a)) = (rowmean(MVi(a)) + colmean(MVj(a)))/2;
    end
    ssqnan = sum(sum(X(MV).^2));
    ssq_old = 10*ssqnan;
    while (abs((ssq_old - ssqnan)/ssq_old) > options(2)) & (iter < options(3))
        ssq_old = ssqnan;
        if options(1)
            [X,mx] = meanc(X);
        end
        [U,D,V] = svd(X,0);
        U = U(:,1:nF);
        D = D(1:nF,1:nF);
        V = V(:,1:nF);
        Xrec = U*D*V';
        if options(1)
            X = meanc(X,mx,1);
            Xrec = meanc(Xrec,mx,1);
        end
        X(MV) = Xrec(MV);
        ssqnan = sum(sum(X(MV).^2));
        iter = iter + 1;
    end
end
if iter == options(3)
   s = ['WARNING: maximum number of iterations (' num2str(options(3)) ') reached before convergence'];
   disp(s)
end

function [Zm,mz] = meanc(Z,mz,rev)
% 020827 FvdB
[n,m] = size(Z);
if nargin == 1
   if ~sum(sum(isnan(Z)))
      mz = mean(Z);
      Zm = Z-mz(ones(n,1),:);
   else
      Zmv = sparse(isnan(Z));
      zmv = n./(n-sum(Zmv));
      Z(Zmv) = 0;
      mz = mean(Z).*zmv;
      Zm = Z-mz(ones(n,1),:);
      Zm(Zmv) = NaN;
   end
elseif nargin == 2
   Zm = Z-mz(ones(n,1),:);
elseif nargin == 3
   Zm = Z+mz(ones(n,1),:);
end