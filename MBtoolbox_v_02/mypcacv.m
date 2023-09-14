function [rmspcv,diagnos] = mypcacv(X,nF,options)
% function [rmspcv,diagnos] = mypcacv(X,nF,options)
% 040706 FvdB
% Performs Principal Component Analysis cross validation to determine the 'true' rank of an X-block.
% Approximately 2.5% of missing values ('diagonal' segments) are used for reconstruction to separate structure from noise.
%
% in:
% X (object x variables) X-data table
% nF (1 x 1) number of factors/PC's used in validation
% options (1 x 6) tolerance for convergence, maximum number of iterations, c-segments,
%    use SVDNaN-imputation i.s.o. NIPALS to estimate missing values, output all cross validation scores and loadings
%    do exastive leave-one-element-out cv (default [1e-8 2000 40 1=yes/SVDNaN 0=no 0=no])
%
% out:
% rmspcv (2 x nF) mean and standard deviation over the 'segments'
% diagnos (structure) cross validation diagnostics: 
%    nseg (1 x 1) number of segments
%    rmspcv_seg (nseg x nF) mean reconstruction error per cv-segment
%    fracnan (nseg x 1) fraction missing values per cv-segment
%    ssq (nF x 1) nominal explained variance per factor
%    s_ssq (nF x 1) standard deviation in ssq
%    T (objects x nF) nominal scores per factor
%    sT (objects x 1) standard deviation in scores
%    P (variables x nF) nominal loadings per factor
%    sP (variables x nF) standard deviation in loadings
%    ssqcv, Tcv, Pcv (... x ... x nseg) output for all cv-segments (if option(5) = 1)
%
% uses:
% mypca.m, svdnan.m

if nargin < 1
    help mypcacv
    return
elseif nargin == 2
    options = [1e-8 2000 40 1 0]; 
end

[nX,mX] = size(X);
options(3) = min([nX mX options(3)]);

if nargout == 2
    if ~options(6)
        diagnos.nseg = options(3);
        diagnos.fracnan = zeros(options(3),1);
    end
    [diagnos.T,diagnos.P,diagnos.ssq] = mypca(X,nF,options([1 2 4]));
end

rmspcv = zeros(2,nF);

if options(6)
    Xnmv = find(~isnan(X));
    diagnos.nseg = length(Xnmv);
    rmspcv_seg = zeros(diagnos.nseg,nF);
    for a=1:diagnos.nseg
        if (a == 1)
            tic
        end
        x = X(Xnmv(a));
        X(Xnmv(a)) = NaN;
        [T,P,ssq] = mypca(X,nF,options([1 2 4]));
        X(Xnmv(a)) = x;
        for aa=1:nF
            Xr = T(:,1:aa)*P(:,1:aa)';
            rmspcv_seg(a,aa) = sqrt((X(Xnmv(a))-Xr(Xnmv(a))).^2);
        end
        if nargout == 2
            ssqcv(:,a+1) = ssq;
            Tcv(:,:,a+1) = T;
            Pcv(:,:,a+1) = P;
        end
        if (a == 1)
            t = toc;
            t = (t*length(Xnmv)/60);
            s = ['Estimated time for all ' num2str(diagnos.nseg) ' trails: ' num2str(t,1) ' minutes'];
            disp(s);
        end
    end
else
    rmspcv_seg = zeros(options(3),nF);
    for a=0:options(3)-1
        if (a == 0)
            tic
        end
        index = -nX+1+a:options(3):mX;
        A = find(spdiags(ones(nX,length(index)),index,nX,mX));
        Xmv = X;
        Xmv(A(:)) = NaN;
        [T,P,ssq] = mypca(Xmv,nF,options([1 2 4]));
        AA = A(find(~isnan(X(A))));
        for aa=1:nF
            Xr = T(:,1:aa)*P(:,1:aa)';
            rmspcv_seg(a+1,aa) = sqrt(mean((X(AA)-Xr(AA)).^2));
        end
        if nargout == 2
            diagnos.fracnan(a+1) = sum(isnan(Xmv(:)))/(nX*mX);
            ssqcv(:,a+1) = ssq;
            Tcv(:,:,a+1) = T;
            Pcv(:,:,a+1) = P;
        end
        if (a == 0)
            t = toc;
            t = (t*options(3)/60);
            if (options(3) > 1) & (t > 0.5)
                s = ['Estimated time for all ' num2str(options(3)) ' trails: ' num2str(t,1) ' minutes'];
                disp(s);
            end
        end
    end
end

if nargout == 2
    diagnos.rmspcv_seg = rmspcv_seg;
    diagnos.s_ssq = std(ssqcv,0,2);
    diagnos.sT = std(Tcv,0,3);
    diagnos.sP = std(Pcv,0,3);
    if options(5)
        diagnos.ssqcv = ssqcv;
        diagnos.Tcv = Tcv;
        diagnos.Pcv = Pcv;
    end
end

rmspcv(1,:) = mean(rmspcv_seg);
rmspcv(2,:) = std(rmspcv_seg);