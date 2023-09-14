function MB = mbpca(X,nF,Xin,Xpp,options)
% function MB = mbpca(X,nF,Xin,Xpp,options)
% 040927 FvdB
% Computes (consensus) multi-block PCA model.
%
% in:
% X (objects x all variables) single, augmented data-block
% nF (1 x 1) number of factors/principal components extracted
% Xin (list) structure index for variables per X-block:
%    e.g. Xin={1:50 51:100} is 2 blocks: X1(1:50,:), X2(51:100,:)
% Xpp (number of blocks x 2) per row (X-block) 
%    1: block scaling (-1 = interactive, 0 = none, 1 = mean center, 2 = autoscale, 3 = range 0 to 1 scale)
%    2: scaled X-block norm (-1 = interactive, 0 = do not change, number = new norm)
% options (1 x 5)
%    1: cross validation (0 = no, default; 1 = yes)
%    2: convergence tolerance (default 1e-8)
%    3: maximum number of iterations (defauls 2000)
%    4: number of segments in cross-validation (default 40)
%    5: override PCA with SVDNaN imputation solution for cross validation and missing values (default 1 = yes)
%
% out:
% MB (struct) Multi-Block record with results: model (type of model), nF (number of factors), Xin (X block index), 
%    Xpp (X preprocessing per block), options, Tt (super scores objects x nF), ssq (explained variance nF x blocks), 
%    Pb (block loadings variables x nF), Tb (block scores objects x nF x blocks), Wt (block weights blocks x nF),
%    Lto (leverage objects super level objects x nF), rmspcv (cross validation prediction errors + std), cvdiag (cv diagnostics)
%
% uses:
% blocknorm.m, meanc.m, autosc.m, rangesc.m, mypcacv.m, mypca.m, svdnan.m

if (nargin < 3)
    help mbpca
    return
end
MB.model = 'mbpca';
MB.nF = nF;
MB.Xin = Xin;
[nX,mX] = size(X);
nbX = size(MB.Xin,2);

if nargin == 3
    MB.Xpp = -ones(nbX,2);
    MB.options = [0 1e-8 2000 40 1];
elseif nargin == 4
    MB.Xpp = Xpp;
    MB.options = [0 1e-8 2000 40 1];
    if size(MB.Xpp,1) ~= nbX
        error('ERROR: number of X-blocks in "Xin" and "Xpp" is not the same')
    elseif size(MB.Xpp,2) ~= 2
        error('ERROR: number of columns in "Xpp" must be two')
    end
else
    MB.Xpp = Xpp;
    MB.options = [0 1e-8 2000 40 1];
    MB.options(1:length(options)) = options;
end
clear nF Xin Xpp options

MV = 0;
ssqX = zeros(1,nbX+1);
Xcoli = [];
for a=1:nbX
    disp(' ');
    coli = MB.Xin{a};
    Xcoli = [Xcoli coli];
    s = ['Block #' num2str(a) ', ' num2str(length(coli)) ' variables'];
    disp(s);
    if any(coli > mX)
        error('ERROR: block index is outside of X-block range')
    end
    Xmv = sparse(isnan(X(:,coli)));
    pXmv = sum(sum(Xmv))/(nX*length(coli))*100;
    s = [num2str(pXmv) '% missing values'];
    disp(s);
    if pXmv ~= 0
        MV = 1;
    end
    if MB.Xpp(a,1) == -1
        inp = input('block scaling (0=none, 1=mean center, 2=autoscale, 3=range(0-1)scale)? : ');
        MB.Xpp(a,1) = inp;
    else
        inp = MB.Xpp(a,1);
    end
    switch inp
        case 0
            disp('no block scaling')
        case 1
            X(:,coli) = meanc(X(:,coli));
            disp('block mean centering')
        case 2
            X(:,coli) = autosc(X(:,coli));
            disp('block autoscaling')
        case 3
            X(:,coli) = rangesc(X(:,coli));
            disp('block range scaling')
        otherwise
            error('ERROR: block scaling must be 0(none), 1(mean center), 2(autoscale) or 3(range(0-1)scale)');
    end
    if MB.Xpp(a,2) == -1
        inp = input('new block norm (0 = do not change)? : ');
        MB.Xpp(a,2) = inp;
    end
    if MB.Xpp(a,2) > 0
        X(:,coli) = blocknorm(X(:,coli),MB.Xpp(a,2));
    end
    if pXmv ~= 0
        Xtemp = X(:,coli);
        ssqX(a+1) = sum(sum(Xtemp(~Xmv).^2));
    else
        ssqX(a+1) = sum(sum(X(:,coli).^2));
    end
    ssqX(1) = ssqX(1) + ssqX(a+1);
    s = ['block-norm (sum-of-squares) ' num2str(ssqX(a+1))];
    disp(s)
end

if MV
    Xmv = sparse(isnan(X));
end
if (MB.options(5) == -1) & MV
    disp(' ')
    inp = input('X-data contains missing values. Do you want to use SVDNaN in stead of NIPALS-PCA (0=no, 1=yes)? : ');
    if inp == 1
        MB.options(5) = 1;
    else
        MB.options(5) = 0;
    end
else
    MB.options(5) = 1;
end

Xstore = X;
if MB.options(1) == 0
    disp('No cross validation')
else
    if MB.options(5) == 1
        disp('Cross validation with SVDNaN imputations')
    else
        disp('Cross validation')
    end
    [MB.rmspcv,MB.cvdiag] = mypcacv(X(:,Xcoli),MB.nF,[MB.options(2:5) 1 0]);
end

[MB.Tt,P,MB.ssq] = mypca(X(:,Xcoli),MB.nF,MB.options([2 3 5]));
for a=1:MB.nF
    for aa=1:nbX
        coli = MB.Xin{aa};
        if MV
            X(Xmv) = 0;
            MB.Pb(coli,a) = X(:,coli)'*MB.Tt(:,a);
            for aaa=coli
                c = (MB.Tt(~Xmv(:,aaa),a)'*MB.Tt(~Xmv(:,aaa),a));
                if (abs(c) > eps)
                    MB.Pb(aaa,a) = MB.Pb(aaa,a)/c;
                end
            end
            tempnorm = norm(MB.Pb(coli,a));
            MB.Pb(coli,a) = MB.Pb(coli,a)/tempnorm;
            MB.Tb(:,a,aa) = X(:,coli)*MB.Pb(coli,a);
            for aaa=1:nX
                c = (MB.Pb(~Xmv(aaa,coli),a)'*MB.Pb(~Xmv(aaa,coli),a));
                if (abs(c) > eps)
                    MB.Tb(aaa,index) = MB.Tb(aaa,a,aa)/c;
                end
            end
            MB.Pb(coli,a) = MB.Pb(coli,a)*tempnorm;
            X(:,coli) = X(:,coli) - MB.Tt(:,a)*MB.Pb(coli,a)';
            X(Xmv) = 0;
        else
            MB.Pb(coli,a) = X(:,coli)'*MB.Tt(:,a)/(MB.Tt(:,a)'*MB.Tt(:,a));
            tempnorm = norm(MB.Pb(coli,a));
            MB.Pb(coli,a) = MB.Pb(coli,a)/tempnorm;
            MB.Tb(:,a,aa) = X(:,coli)*MB.Pb(coli,a)/(MB.Pb(coli,a)'*MB.Pb(coli,a));
            MB.Pb(coli,a) = MB.Pb(coli,a)*tempnorm;
            X(:,coli) = X(:,coli) - MB.Tt(:,a)*MB.Pb(coli,a)';
        end
        MB.ssq(a,aa+1) = (ssqX(aa+1) - sum(sum(X(:,coli).^2)))/ssqX(aa+1);
    end
    index = (a-1)*nbX+1:(a-1)*nbX+nbX;
    MB.Wt(:,a) = MB.Tb(:,index)'*MB.Tt(:,a)/(MB.Tt(:,a)'*MB.Tt(:,a));
    MB.Wt(:,a) = MB.Wt(:,a)/norm(MB.Wt(:,a)); %!!! super scores are not automatically of length one
    if MV
        X(Xmv) = 0;
    end
    MB.Lto(:,a) = diag(MB.Tt(:,1:a)*pinv(MB.Tt(:,1:a)'*MB.Tt(:,1:a))*MB.Tt(:,1:a)');
end

if MB.options(1)
    for acv=1:MB.cvdiag.nseg
        X = Xstore;
        for a=1:MB.nF
            for aa=1:nbX
                coli = MB.Xin{aa};
                if MV
                    X(Xmv) = 0;
                    Pb(coli,a) = X(:,coli)'*MB.cvdiag.Tcv(:,a,acv);
                    for aaa=coli
                        c = (MB.cvdiag.Tcv(~Xmv(:,aaa),a,acv)'*MB.cvdiag.Tcv(~Xmv(:,aaa),a,acv));
                        if (abs(c) > eps)
                            Pb(aaa,a) = Pb(aaa,a)/c;
                        end
                    end
                    tempnorm = norm(Pb(coli,a));
                    Pb(coli,a) = Pb(coli,a)/tempnorm;
                    index = (a-1)*nbX+aa;
                    Tb(:,index) = X(:,coli)*Pb(coli,a);
                    for aaa=1:nX
                        c = (Pb(~Xmv(aaa,coli),a)'*Pb(~Xmv(aaa,coli),a));
                        if (abs(c) > eps)
                            Tb(aaa,index) = Tb(aaa,index)/c;
                        end
                    end
                    Pb(coli,a) = Pb(coli,a)*tempnorm;
                    X(:,coli) = X(:,coli) - MB.cvdiag.Tcv(:,a,acv)*Pb(coli,a)';
                    X(Xmv) = 0;
                else
                    Pb(coli,a) = X(:,coli)'*MB.cvdiag.Tcv(:,a,acv)/(MB.cvdiag.Tcv(:,a,acv)'*MB.cvdiag.Tcv(:,a,acv));
                    tempnorm = norm(Pb(coli,a));
                    Pb(coli,a) = Pb(coli,a)/tempnorm;
                    index = (a-1)*nbX+aa;
                    Tb(:,index) = X(:,coli)*Pb(coli,a)/(Pb(coli,a)'*Pb(coli,a));
                    Pb(coli,a) = Pb(coli,a)*tempnorm;
                    X(:,coli) = X(:,coli) - MB.cvdiag.Tcv(:,a,acv)*Pb(coli,a)';
                end
            end
            index = (a-1)*nbX+1:(a-1)*nbX+nbX;
            MB.cvdiag.Wt_cv(:,a,acv) = Tb(:,index)'*MB.cvdiag.Tcv(:,a,acv)/(MB.cvdiag.Tcv(:,a,acv)'*MB.cvdiag.Tcv(:,a,acv));
            MB.cvdiag.Wt_cv(:,a,acv) = MB.cvdiag.Wt_cv(:,a,acv)/norm(MB.cvdiag.Wt_cv(:,a,acv)); %!!! super scores are not automatically of length one
            if MV
                X(Xmv) = 0;
            end
        end
    end
    MB.cvdiag.sWt = std(MB.cvdiag.Wt_cv,0,3);
end