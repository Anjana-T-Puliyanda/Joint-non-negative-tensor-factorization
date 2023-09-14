function [Xn,Yn,diagnos] = procrusND(X,Y,options)
% function [Xn,Yn,diagnos] = procrusND(X,Y,options)
% 040111 FvdB
% Procrustes (full) matching for two 'images' X and Y.
% Order: 1) center 2) rotate 3) resize
% Uses zero padding if number of columns (variables) in X and Y are different. 
% In case of padding uniform scaling (see options) is turned of.
%
% in:
% X (n x m) coordinates of landmarks x dimensionality (variables) in X-matrix (the target matrix)
% Y (n x m) coordinates of landmarks x dimensionalety (variables) in Y-matrix 
% options (1 x 2) 1: perform uniform scaling (default 1 = yes)
%                 2: norm matrices X and Y to 2-norm = 1 (default 0 = no)
%
% out:
% Xn (n x 2) new X-matrix (column-centered)
% Yn (n x 2) new Y-matrix (column-centered, rotated and 'match-scaled')
% diagnos (struct) Procrustes step-diagnostics
%
% based on: G.H. Golub and C.F. van Loan 'Matrix Computations' 3rd(1996)
%           J.C. Gower 'Generalized Procrustes Analysis' Psychometrika 40/1(1975)33-51

if nargin < 2
    help procrusND
    return
elseif nargin == 2
    options = [1 0];
end

[nX,mX] = size(X);
[nY,mY] = size(Y);

if nX ~= nY
    s = ['ERROR: number of rows in X ' num2str(nX) ' and Y ' num2str(nY) ' must be the same'];
    error(s)
end
if (mX < mY)
    X = [X zeros(nX,mY-mX)];
    options = 0;
elseif (mX > mY)
    Y = [Y zeros(nX,mX-mY)];
    options = 0;
end
if any(isnan([X Y]))
    error('ERROR: function does not work with missing values (NaN)');
end

diagnos.options = options;
diagnos.cX = mean(X);
diagnos.cY = mean(Y);
Xn = (X-ones(nX,1)*diagnos.cX);
Yn = (Y-ones(nY,1)*diagnos.cY);
if options(2)
    diagnos.sXpre = norm(X);
    diagnos.sYpre = norm(Y);
    Xn = Xn/diagnos.sXpre;
    Yn = Yn/diagnos.sYpre;
end

[U,D,V] = svd(Yn'*Xn);
diagnos.Q = U*V';
if options(1)
    diagnos.r = trace(Yn*diagnos.Q*Xn')/trace(Yn*Yn');
    Yn = Yn*diagnos.Q*diagnos.r;
else
    Yn = Yn*diagnos.Q;
end

for a=1:mX
    diagnos.pre_congruence(a) = (X(:,a)'*Y(:,a))/(sqrt(X(:,a)'*X(:,a))*sqrt(Y(:,a)'*Y(:,a))); 
    diagnos.post_congruence(a) = (Xn(:,a)'*Yn(:,a))/(sqrt(Xn(:,a)'*Xn(:,a))*sqrt(Yn(:,a)'*Yn(:,a))); 
end

if (det(diagnos.Q) > 0)
    if options(1)
        diagnos.Procrus_operation = 'rotation and uniform scaling';
    else
        diagnos.Procrus_operation = 'rotation and no scaling';
    end
else
    if options(1)
        diagnos.Procrus_operation = 'rotation, reflection and uniform scaling';
    else
        diagnos.Procrus_operation = 'rotation, reflection and no scaling';
    end
end