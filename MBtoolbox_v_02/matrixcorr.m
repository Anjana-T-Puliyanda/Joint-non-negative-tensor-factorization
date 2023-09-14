function [R,nF] = matrixcorr(X,Y,corr,nF)
% function [R,nF] = matrixcorr(X,Y,corr,nF)
% 040630 FvdB
% Computes different matrix correlation coefficients
%
% in:
% X (objects x m-variables) data table one
% Y (objects x p-variables or cell array) data table two or cell array with "n" indices in "X"
% corr (1 x 1) type of correlation to
%    1 : spectrum and orientation dependent: spectrum, row and column space
%    2 : orientation independent: spectrum and column space
%    3 : spectrum independent: row and column space
%    4 : spectrum independent: column space
%    5 : RV-coefficient: association matrix reconstructed from max nF (to eliminate missing values) (default)
% nF (1 x 1) number of singular values used in (truncated rank) matrix correlation computation 
%    (default min(objects,m1,m2))
%
% out:
% R (1 x 1 or n x n) correlation coefficient for two input matrices or correlation matrix for "n" entries in "X"
% nF (1 x 1) number of singular values used in computation matrix correlation
%
% uses:
% svdnan.m, meanc.m
%
% based on:
% J.O.Ramsay, J.ten Berge and G.P.H. Styan 'Matrix Correlation' Psychometrika 49/3(1984)403-423
% E.M.Qannari, I.Wakeling and H.J.H MacFie 'A Hierarchy of Models for Analysing Sensory Data' Food Quality and Prefrence 6(1995)309-314

if nargin < 2
    help matrixcorr
    return
elseif nargin == 2
    corr = 5;
    nF = 0;
elseif nargin == 3
    nF = 0;
elseif nargin == 4 & isempty(corr)
    corr = 5;
end

[nX,mX] = size(X);
if iscell(Y)
    nXin = length(Y);
    for a=1:nXin
        lXin(a) = length(Y{a});
        if any(Y{a} < 0) | any(Y{a} > mX)
            error('ERROR: index in "X" (in cell "Y") is out of range for "X"')
        end
    end
    if nF == 0
        nF = min([nX lXin]);
    elseif nF > min([nX lXin])
        error('ERROR: "nF" must be smaller than/equal to the minimum of the number of objects and variables in the matrices')
    end
    R = eye(nXin);
    for a=1:nXin-1
        if (corr == 5)
            [U1,D1,V1] = svdnan(X(:,Y{a}),min(size(X(:,Y{a}))),[1 1e-8 2000]);
        else
            [U1,D1,V1] = svdnan(X(:,Y{a}),nF,[1 1e-8 2000]);
        end
        for aa=a+1:nXin
            if (corr == 5)
                [U2,D2,V2] = svdnan(X(:,Y{aa}),min(size(X(:,Y{aa}))),[1 1e-8 2000]);
            else
                [U2,D2,V2] = svdnan(X(:,Y{aa}),nF,[1 1e-8 2000]);
            end
            switch corr
            case 1
                Rt = V1*D1*U1'*U2*D2*V2';
                if any(size(Rt)==1)
                    Rt = Rt(1,1);
                else
                    Rt = trace(Rt);
                end
                R(a,aa) = Rt/(trace(D1.^2)*trace(D2.^2))^0.5;
                R(a,aa) = abs(R(a,aa));
            case 2
                for aaa=1:nF
                    if (U1(:,aaa)'*U2(:,aaa)) < 0
                        U1(:,aaa) = -U1(:,aaa);
                    end
                end
                Rt = D1*U1'*U2*D2;
                if any(size(Rt)==1)
                    Rt = Rt(1,1);
                else
                    Rt = trace(Rt);
                end
                R(a,aa) = Rt/(trace(D1.^2)*trace(D2.^2))^0.5;
            case 3
                Rt = V1*U1'*U2*V2';
                if any(size(Rt)==1)
                    Rt = Rt(1,1);
                else
                    Rt = trace(Rt);
                end
                R(a,aa) = Rt/nF;
                R(a,aa) = abs(R(a,aa));
            case 4
                for aaa=1:nF
                    if (U1(:,aaa)'*U2(:,aaa)) < 0
                        U1(:,aaa) = -U1(:,aaa);
                    end
                end
                Rt = U1'*U2;
                if any(size(Rt)==1)
                    Rt = Rt(1,1);
                else
                    Rt = trace(Rt);
                end
                R(a,aa) = Rt/nF;
            case 5
                W1 = U1*D1*V1'*V1*D1'*U1';
                W2 = U2*D2*V2'*V2*D2'*U2';
                R(a,aa) = trace(W1*W2)/(sqrt(trace(W1*W1))*sqrt(trace(W2*W2)));
            otherwise
                error('ERROR: "corr" must be an integer in the range 1 to 5')
            end
            R(aa,a) = R(a,aa);
        end
    end
else
    [nY,pY] = size(Y);
    if nX ~= nY
        error('ERROR: Number of objects in the two data tables must be the same')
    end
    if nF == 0
        nF = min([nX mX pY]);
    elseif nF > min([nX mX pY])
        error('ERROR: "nF" must be smaller than/equal to the minimum of the number of objects and variables in the two matrices')
    end
    if (corr == 5)
        [U1,D1,V1] = svdnan(X,min(size(X)),[1 1e-8 2000]);
        [U2,D2,V2] = svdnan(Y,min(size(Y)),[1 1e-8 2000]);
    else
        [U1,D1,V1] = svdnan(X,nF,[1 1e-8 2000]);
        [U2,D2,V2] = svdnan(Y,nF,[1 1e-8 2000]);
    end
    switch corr
    case 1
        R = V1*D1*U1'*U2*D2*V2';
        if any(size(R)==1)
            R = R(1,1);
        else
            R = trace(R);
        end
        R = R/(trace(D1.^2)*trace(D2.^2))^0.5;
        R = abs(R);
    case 2
        for a=1:nF
            if (U1(:,a)'*U2(:,a)) < 0
                U1(:,a) = -U1(:,a);
            end
        end
        R = D1*U1'*U2*D2;
        if any(size(R)==1)
            R = R(1,1);
        else
            R = trace(R);
        end
        R = R/(trace(D1.^2)*trace(D2.^2))^0.5;
    case 3
        R = V1*U1'*U2*V2';
        if any(size(R)==1)
            R = R(1,1);
        else
            R = trace(R);
        end
        R = R/nF;
        R = abs(R);
    case 4
        for a=1:nF
            if (U1(:,a)'*U2(:,a)) < 0
                U1(:,a) = -U1(:,a);
            end
        end
        R = U1'*U2;
        if any(size(R)==1)
            R = R(1,1);
        else
            R = trace(R);
        end
        R = R/nF;
    case 5
        W1 = U1*D1*D1'*U1';
        W2 = U2*D2*D2'*U2';
        R = trace(W1*W2)/(sqrt(trace(W1*W1))*sqrt(trace(W2*W2)));
    otherwise
        error('ERROR: "corr" must be an integer in the range 1 to 5')
    end
end
if corr == 5
    nF = -1;
end