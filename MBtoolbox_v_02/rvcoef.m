function RV = rvcoef(X1,X2)
% function RV = rvcoef(X1,X2)
% 040630 FvdB
% Computes RV-coefficient, the correlation coefficient between two data
% tables
%
% in:
% X1 (objects x variables1) data table (augmented if X2 is structure)
% X2 (objects x variables2 or x-index cell) data tables or index into augmented data table "X1"
%    e.g. is X1 = (10 x 25) and X2 = {1:5 6:10 11:25} "RV" will be 3x3 for the three matrices
%
% out:
% RV (1 x 1 or n x n) RV-coefficient or RV-coef.-matrix of size n = number of indexes in X2
%    e.g. RV = 0 is no overlap in subspaces of X1 and X2, RV = 1 same subspaces
%
% based on:
% E.Qannari, I.Wakeling and J.MacFie 'A hierarchy of models for analysing sensory data' 
% Food Quality and Preference 6(1995)309-314
% uses:
% meanc.m

if nargin < 2
    help rvcoef
    return
end

if any(isnan(X1(:)))
    error('ERROR: "rvcoef" can not handle missing values in "X1", try "matrixcorr" function instead');
end

if iscell(X2)
    [nX1,mX1] = size(X1);
    n = length(X2);
    for a=1:n
        index = X2{n};
        if any(index < 0) | any(index > mX1)
            error('ERROR: index in "X2" (struct) is outside of the column range for "X1"')
        end
    end
    RV = eye(n);
    for a=1:n-1
        for aa=2:n
            index1 = X2{a};
            index2 = X2{aa};
            W1 = meanc(X1(:,index1));
            W2 = meanc(X1(:,index2));
            W1 = W1*W1';
            W2 = W2*W2';
            RV(a,aa) = trace(W1*W2)/(sqrt(trace(W1*W1))*sqrt(trace(W2*W2)));
            RV(aa,a) = RV(a,aa);
        end
    end
else
    if any(isnan(X2(:)))
        error('ERROR: "rvcoef" can not handle missing values in "X2", try "matrixcorr" function instead');
    end
    [nX1,mX1] = size(X1);
    [nX2,mX2] = size(X2);
    if nX1 ~= nX2
        error('ERROR: number of objects in "X1" and "X2" must be equal')
    end
    X1 = meanc(X1);
    X2 = meanc(X2);
    W1 = X1*X1';
    W2 = X2*X2';
    RV = trace(W1*W2)/(sqrt(trace(W1*W1))*sqrt(trace(W2*W2)));
end