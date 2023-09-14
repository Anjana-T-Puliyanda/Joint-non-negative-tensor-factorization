function [T,P,ssq,Ro,Rv,Lo,Lv] = mypca(X,nF,options)
% function [T,P,ssq,Ro,Rv,Lo,Lv] = mypca(X,nF,options)
% 040706 FvdB
% Principal Component Analysis bilinear factor model.
%
% in :
% X (objects x variables) data block
% nF (1 x 1) number of factors/principal components
% options (1 x 3) tolerance for convergence, maximum number of iterations 
%    and use SVD/SVDNaN-imputation i.s.o. NIPALS to estimate parameters (default 1e-8, 2000 and 1=yes/SVDNaN)
%
% out : 
% T (objects x nF) scores
% P (variables x nF) loadings
% ssq (nF x 1) cumulative sum of squares
% Ro (objects x nF) object residuals
% Rv (variables x nF) variable residuals
% Lo (objects x nF) object leverages
% Lv (variables x nF) variable leverages
%
% uses:
% svdnan.m

if nargin < 2
    help mypca
    return
elseif nargin == 2
    options = [1e-8 2000 1];
else
    if (length(options) == 1)
        options = [options 2000 1];
    elseif (length(options) == 2)
        options = [options 1];
    end
end

[nX,mX] = size(X);
if nF > nX
    s = ['ERROR: number of objects (' int2str(nX) ') is to small to compute ' int2str(nF) ' factors'];
    error(s)
end
if nF > mX
    s = ['ERROR: number of variables (' int2str(mX) ') is to small to compute ' int2str(nF) ' factors'];
    error(s)
end

T = zeros(nX,nF);
P = zeros(mX,nF);
ssq = zeros(nF,1);
MV = sum(sum(isnan(X)));
if options(3)
    [T,D,P] = svdnan(X,nF,[0 options(1:2)]);
    T = T*D;
    if MV
        Xmv = sparse(isnan(X));
        X(Xmv) = 0;
    end
    ssqX = sum(sum(X.^2));
    for a=1:nF
        X = X - T(:,a)*P(:,a)';
        if MV
            X(Xmv) = 0;
        end
        ssq(a) = (ssqX - sum(sum(X.^2)))/ssqX;
        if (nargout > 3)
            Ro(:,a) = sqrt(sum(X.^2,2));
            Rv(:,a) = sqrt(sum(X.^2,1))';
            Lo(:,a) = diag(T(:,1:a)*pinv(T(:,1:a)'*T(:,1:a))*T(:,1:a)');
            Lv(:,a) = diag(P(:,1:a)*P(:,1:a)');
        end
    end
else
    if MV
        Xmv = sparse(isnan(X));
        X(Xmv) = 0;
    end
    ssqX = sum(sum(X.^2));
    for a=1:nF
        iter = 0;
        [aa,aaa] = max(sum(X.^2,1));
        T(:,a) = X(:,aaa);
        t_old = T(:,a)*10;
        if MV
            while (sum((t_old - T(:,a)).^2)/sum(t_old.^2) > options(1)) & (iter <= options(2))
                iter = iter + 1;
                t_old = T(:,a);
                P(:,a) = X'*T(:,a);
                for aa=1:mX
                    c = (T(~Xmv(:,aa),a)'*T(~Xmv(:,aa),a));
                    if (abs(c) > eps)
                        P(aa,a) = P(aa,a)/c;
                    end
                end
                P(:,a) = P(:,a)/norm(P(:,a));
                T(:,a) = X*P(:,a);
                for aa=1:nX
                    c = (P(~Xmv(aa,:),a)'*P(~Xmv(aa,:),a));
                    if (abs(c) > eps)
                        T(aa,a) = T(aa,a)/c;
                    end
                end
            end
        else
            while (sum((t_old - T(:,a)).^2)/sum(t_old.^2) > options(1)) & (iter < options(2))
                iter = iter + 1;
                t_old = T(:,a);
                P(:,a) = X'*T(:,a)/(T(:,a)'*T(:,a));
                P(:,a) = P(:,a)/norm(P(:,a));
                T(:,a) = X*P(:,a)/(P(:,a)'*P(:,a));
            end
        end
        if iter == options(2)
            s = ['WARNING: maximum number of iterations (' num2str(options(2)) ') reached before convergence'];
            disp(s)
        end
        X = X - T(:,a)*P(:,a)';
        if MV
            X(Xmv) = 0;
        end
        ssq(a) = (ssqX - sum(sum(X.^2)))/ssqX;
        if (nargout > 3)
            Ro(:,a) = sqrt(sum(X.^2,2));
            Rv(:,a) = sqrt(sum(X.^2,1))';
            Lo(:,a) = diag(T(:,1:a)*pinv(T(:,1:a)'*T(:,1:a))*T(:,1:a)');
            Lv(:,a) = diag(P(:,1:a)*P(:,1:a)');
        end
    end
end