function MB = mbpls(X,Y,nF,Xin,Xpp,Ypp,options)
% function MB = mbpls(X,Y,nF,Xin,Xpp,Ypp,options)
% 040927 FvdB
% Computes multi-block PLS model
%
% in:
% X (objects x all variables) single, augmented data-block
% Y (objects x Y-variables) Y-data block
% nF (1 x 1) number of factors/latent variables extracted
% Xin (list) structure index for variables per X-block:
%    e.g. Xin={1:50 51:100} is 2 blocks: X1(1:50,:), X2(51:100,:)
% Xpp (number of blocks x 2) per row (X-block) 
%    1: block scaling (-1 = interactive, 0 = none, 1 = mean center, 2 = autoscale, 3 = range 0 to 1 scale)
%    2: scaled X-block norm (-1 = interactive, 0 = do not change, number = new norm)
% Ypp (1 x 1) Y-block scaling
% options (1 x 5)
%    1: cross validation (0 = no, default; 1 = yes)
%    2: convergence tolerance (default 1e-8)
%    3: maximum number of iterations (defauls 2000)
%    4: not used
%    5: not used
%
% out:
% MB (struct) Multi-Block record with results: model (type of model), nF (number of factors), Xin (X block index), 
%    Xpp (X preprocessing per block), options, Tt (super scores objects x nF), ssq (explained variance nF x blocks), 
%    Pb (block loadings variables x nF), Tb (block scores objects x nF x blocks), Wt (block weights blocks x nF),
%    Lto (leverage objects super level objects x nF), rmspcv (cross validation prediction errors), diagnos (cv diagnostics)
%
% uses:
% blocknorm.m, meanc.m, autosc.m, rangesc.m, myplscv.m, mypls.m

if (nargin < 4)
    help mbpls
    return
end
MB.model = 'mbpls';
MB.nF = nF;
MB.Xin = Xin;
[nX,mX] = size(X);
[nY,pY] = size(Y);
MB.Yin = {1:pY};
if nX ~= nY
    s = ['ERROR: number of objects in X (' num2str(nX) ') and Y-block (' num2str(nY) ') is not the same'];
    error(s);
end
nbX = size(MB.Xin,2);
nbY = 1;
if nargin == 4
    MB.Xpp = -ones(nbX,2);
    MB.Ypp = -ones(nbY,1);
    MB.options = [0 1e-8 2000 40 0];
elseif nargin == 5
    MB.Xpp = Xpp;
    MB.Ypp = -ones(nbY,1);
    MB.options = [0 1e-8 2000 40 0];
    if size(Xpp,1) ~= nbX
        error('ERROR: number of X-blocks in "Xin" and "Xpp" is not the same')
    elseif size(Xpp,2) ~= 2
        error('ERROR: number of columns in "Xpp" must be two')
    end
elseif nargin == 6
    MB.Xpp = Xpp;
    MB.Ypp = Ypp;
    MB.options = [0 1e-8 2000 40 0];
    if size(Ypp,1) ~= nbY
        error('ERROR: number of Y-blocks in "Yin" and "Ypp" is not the same')
    elseif size(Ypp,2) ~= 1
        error('ERROR: number of columns in "Ypp" must be one')
    end
else
    MB.Xpp = Xpp;
    MB.Ypp = Ypp;
    MB.options = [0 1e-8 2000 40 0];
    MB.options(1:length(options)) = options;
end
clear nF Xin Xpp Ypp options

MV = 0;
ssqX = zeros(1,nbX+1);
ssqY = zeros(1,nbY+1);
Xcoli = [];
Ycoli = [];
for a=1:nbX
    disp(' ');
    coli = MB.Xin{a};
    Xcoli = [Xcoli coli];
    s = ['X-Block #' num2str(a)];
    disp(s);
    s = [num2str(length(coli)) ' variables'];
    disp(s)
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
        inp = input('X-block scaling (0=none, 1=mean center, 2=autoscale, 3=range(0-1)scale)? : ');
        MB.Xpp(a,1) = inp;
    else
        inp = MB.Xpp(a,1);
    end
    switch inp
        case 0
            disp('no X-block scaling')
        case 1
            X(:,coli) = meanc(X(:,coli));
            disp('X-block mean centering')
        case 2
            X(:,coli) = autosc(X(:,coli));
            disp('X-block autoscaling')
        case 3
            X(:,coli) = rangesc(X(:,coli));
            disp('X-block range scaling')
        otherwise
            error('ERROR: X-block scaling must be 0(none), 1(mean center), 2(autoscale) or 3(range(0-1)scale)');
    end
    if MB.Xpp(a,2) == -1
        inp = input('new X-block norm (0 = do not change)? : ');
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
    s = ['X-block-norm (sum-of-squares) ' num2str(ssqX(a+1))];
    disp(s)
end

for a=1:nbY
    disp(' ');
    coli = MB.Yin{a};
    Ycoli = [Ycoli coli];
    s = ['Y-Block #' num2str(a)];
    disp(s);
    s = [num2str(length(coli)) ' variables'];
    disp(s)
    if any(coli > pY)
        error('ERROR: Y-block index is outside of Y-block range')
    end
    Ymv = sparse(isnan(Y(:,coli)));
    pYmv = sum(sum(Ymv))/(nX*length(coli))*100;
    s = [num2str(pYmv) '% missing values'];
    disp(s);
    if pYmv ~= 0
        MV = 1;
    end
    if MB.Ypp(a,1) == -1
        inp = input('Y-block scaling (0=none, 1=mean center, 2=autoscale, 3=range(0-1)scale)? : ');
        MB.Ypp(a,1) = inp;
    else
        inp = MB.Ypp(a,1);
    end
    switch inp
        case 0
            disp('no Y-block scaling')
        case 1
            [Y(:,coli),my] = meanc(Y(:,coli));
            disp('Y-block mean centering')
        case 2
            [Y(:,coli),my,sy] = autosc(Y(:,coli));
            disp('Y-block autoscaling')
        case 3
            [Y(:,coli),ry] = rangesc(Y(:,coli));
            disp('Y-block range scaling')
        otherwise
            error('ERROR: Y-block scaling must be 0(none), 1(mean center), 2(autoscale) or 3(range(0-1)scale)');
    end
    if pYmv ~= 0
        Ytemp = Y(:,coli);
        ssqY(a+1) = sum(sum(Ytemp(~Ymv).^2));
    else
        ssqY(a+1) = sum(sum(Y(:,coli).^2));
    end
    ssqY(1) = ssqY(1) + ssqY(a+1);
end

if MB.options(1) == 0
    disp('No cross validation')
else
    disp('Cross validation')
    [MB.rmspcv,MB.Ypcv] = myplscv(X(:,Xcoli),Y(:,Ycoli),MB.nF,[0 0],MB.options(2:4));
    for a=1:MB.nF
        coli = (a-1)*pY+1:(a-1)*pY+pY;
        if MB.Ypp == 1
            MB.Ypcv(:,coli) = meanc(MB.Ypcv(:,coli),my,1);
        elseif MB.Ypp == 2
            MB.Ypcv(:,coli) = autosc(MB.Ypcv(:,coli),my,sy,1);
        elseif MB.Ypp == 3
            MB.Ypcv(:,coli) = rangesc(MB.Ypcv(:,coli),ry,1);
        end
    end
end

[MB.Tt,MB.Pb,MB.W,MB.U,MB.Q,MB.B,MB.ssq] = mypls(X(:,Xcoli),Y(:,Ycoli),MB.nF,MB.options(2:3));

for a=1:MB.nF
    if MV
        X(Xmv) = 0;
        Y(Ymv) = 0;
    end
    for aa=1:nbX
        coli = MB.Xin{aa};
        if MV
            MB.Pb(coli,a) = X(:,coli)'*MB.Tt(:,a);
            for aaa=coli
                MB.Pb(aaa,a) = MB.Pb(aaa,a)/(MB.Tt(~Xmv(:,aaa),a)'*MB.Tt(~Xmv(:,aaa),a));
            end
            MB.Wb(coli,a) = X(:,coli)'*MB.U(:,a);
            for aaa=coli
               MB.Wb(aaa,a) = MB.Wb(aaa,a)/(MB.U(~Ymv(:,aaa),a)'*MB.U(~Ymv(:,aaa),a));
            end
        else
            MB.Pb(coli,a) = X(:,coli)'*MB.Tt(:,a)/(MB.Tt(:,a)'*MB.Tt(:,a));
            MB.Wb(coli,a) = X(:,coli)'*MB.U(:,a)/(MB.U(:,a)'*MB.U(:,a));
        end
        MB.Tb(:,(a-1)*nbX+aa) = X(:,coli)*MB.Wb(coli,a);
        X(:,coli) = X(:,coli) - MB.Tt(:,a)*MB.Pb(coli,a)';
        Y = Y - MB.Tt(:,a)*MB.Q(:,a)';
        MB.ssq(a,aa+2) = (ssqX(aa+1) - sum(sum(X(:,coli).^2)))/ssqX(aa+1);
    end
    MB.Wt(:,a) = MB.Tb(:,(a-1)*nbX+1:(a-1)*nbX+nbX)'*MB.U(:,a)/(MB.U(:,a)'*MB.U(:,a));
    MB.Wt(:,a) = MB.Wt(:,a)/norm(MB.Wt(:,a)); %!!!
    
    for aa=1:nbX
        rowi = MB.Xin{aa};
        coli = (a-1)*nbX+aa;
        MB.Rbo(:,coli) = sqrt(sum(X(:,rowi).^2,2));
        index = aa:nbX:(a-1)*nbX+aa;
        MB.Lbo(:,coli) = diag(MB.Tb(:,index)*pinv(MB.Tb(:,index)'*MB.Tb(:,index))*MB.Tb(:,index)');
    end
    MB.Ry(:,coli) = sqrt(sum(Y.^2,2));
    MB.Rbv(:,a) = sqrt(sum(X.^2,1))';
    MB.Lbv(:,a) = diag(MB.Pb(:,1:a)*MB.Pb(:,1:a)');
    MB.Lto(:,a) = diag(MB.Tt(:,1:a)*pinv(MB.Tt(:,1:a)'*MB.Tt(:,1:a))*MB.Tt(:,1:a)');
end