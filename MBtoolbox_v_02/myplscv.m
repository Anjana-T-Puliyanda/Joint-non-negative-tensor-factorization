function [rmspcv,Ypcv,diagnos] = myplscv(X,Y,nF,prepro,options,cvseg)
% function [rmspcv,Ypcv,diagnos] = myplscv(X,Y,nF,prepro,options,cvseg)
% 041104 FvdB
% Partial least Squares regression bi-linear model cross-validation
%
% in : 
% X (objects x m-variables) data-block
% Y (objects x p-variables) data-block
% nF (1 x 1) number of factors (latent variables)
% prepro (1 x 1) or (1 x 2) preprocessing of (X+Y)-block or X and Y-block:
%    0 = none, 1 = mean center, 2 = auto scale, 3 = 0-1 range scale
% options (1 x 2) 1: tolerance for convergence (default 1e-8)
%                 2: maximum number of iterations (default 2000)
% cvseg (objects x 1) defines segments for cross-validation, starting with 1, using consequtive numbers
%    if 'cvseg' is not specified: leave-one-out-cv is used, otherwise e.g. cvseg = [1 1 2 3 3 2]' 
%    predicts samples 1 and 2 in the first cv-loop, 3 and 6 in the second and 4 and 5 in the third 
%    (use e.g. to leave out replicates in one cross validation cycle).
%
% out : 
% rmspcv (1 x nF) root mean squared error of prediction for all Y-block variables
% Ypcv (objects x p.nF) cross validation Y-block predictions
% diagnos (structure) cross validation diagnostics
%    nseg (1 x 1) number of segments
%    s_rmspcv (1 x nF or p x nF) standard deviation in rmspcv (deviation from nominal)
%    ssq (nF x 1) nominal explained variance per factor
%    r_ssq (nF x 4) X and Y range in ssq over cvsegments (nF x [Xmin Ymin Xmax Ymax])
%    B (m-variables x p.nF) nominal b-regression coefficients for all objects (p.nF entries)
%    s_B (m-variables x p.nF) Jack-knife based b-regression coefficient standard deviations over cvsegments
%    T (n-object x nF) nominal scores 
%    s_T (n-objects x nF) Jack-knife based score standard deviations over cvsegments
%    P (m-variables x nF) nominal loadings
%    s_P (m-variables x nF) Jack-knife based loading standard deviations over cvsegments
%
% uses: meanc.m, autosc.m, rangesc.m, mypls.m

if (nargin < 3)
    help myplscv
    return
end
if length(prepro) == 1
    prepro = [prepro prepro];
end
if nargin == 4
    options = [1e-8 2000];
    loo_cv = 1;
end
if nargin == 5
    loo_cv = 1;
end
if nargin == 6
    loo_cv = 0;
end
if isempty(options)
    options = [1e-8 2000];
end

[nX,mX] = size(X);
[nY,pY] = size(Y);
if nX ~= nY
    s = ['number of objects in X (' int2str(nX) ') and Y (' int2str(nY) ') must be the same'];
    error(s)
end
if ~loo_cv & (length(cvseg) ~= nX)
    s = ['number of objects (' int2str(nX) ') and entries in "cvsegments" (' int2str(length(cvseg)) ') must be the same'];
    error(s)
end

if prepro(1) == 1
    Xf = meanc(X);
elseif prepro(1) == 2
    Xf = autosc(X);
elseif prepro(1) == 3
    Xf = rangesc(X);
else
    Xf = X;
end
if prepro(2) == 1
    Yf = meanc(Y);
elseif prepro(2) == 2
    Yf = autosc(Y);
elseif prepro(2) == 3
    Yf = rangesc(Y);
else
    Yf = Y;
end
[Tf,Pf,Wf,Uf,Qf,Bf,ssqf] = mypls(Xf,Yf,nF,options);

if loo_cv
    nseg = nX;
    Xc = zeros(nX-1,mX);
    Yc = zeros(nX-1,pY);
    Ypcv = zeros(nX,pY*nF);
    Btemp = zeros(mX,pY*nF,nseg);
    r_ssq = zeros(nF,2,nseg);
    Ttemp = zeros(nX,nF,nseg);
    Ptemp = zeros(mX,nF,nseg);
    Qtemp = zeros(pY,nF,nseg);
    for a=1:nX
        indexc = [1:a-1 a+1:nX];
        Xc = X(indexc,:);
        Yc = Y(indexc,:);
        xp = X(a,:);
        yp = Y(a,:);
        if prepro(1) == 1
            [Xc,mx] = meanc(Xc);
            xp = meanc(xp,mx);
        elseif prepro(1) == 2
            [Xc,mx,stdx] = autosc(Xc);
            xp = autosc(xp,mx,stdx);
        elseif prepro(1) == 3
            [Xc,rx] = rangesc(Xc);
            xp = rangesc(xp,rx);
        end
        if prepro(2) == 1
            [Yc,my] = meanc(Yc);
        elseif prepro(2) == 2
            [Yc,my,stdy] = autosc(Yc);
        elseif prepro(2) == 3
            [Yc,ry] = rangesc(Yc);
        end
        [Ttemp(indexc,:,a),Ptemp(:,:,a),W,U,Qtemp(:,:,a),Btemp(:,:,a),r_ssq(:,:,a)] = mypls(Xc,Yc,nF,options);
        xpmv = isnan(xp);
        xp(xpmv) = 0;
        for aa=1:nF
            index = (aa-1)*pY+1:aa*pY;
            if aa > 1
                Ypcv(a,index) = Ypcv(a,index-pY);
            end
            tp = xp*W(:,aa)/(W(~xpmv,aa)'*W(~xpmv,aa));
            Ypcv(a,index) = Ypcv(a,index) + tp*Qtemp(:,aa,a)';
            xp = xp - tp*Ptemp(:,aa,a)';
            Ttemp(a,aa,a) = tp;
        end
        for aa=1:nF
            index = (aa-1)*pY+1:aa*pY;
            if prepro(2) == 1
                Ypcv(a,index) = meanc(Ypcv(a,index),my,1);
            elseif prepro(2) == 2
                Ypcv(a,index) = autosc(Ypcv(a,index),my,stdy,1);
            elseif prepro(2) == 3
                Ypcv(a,index) = rangesc(Ypcv(a,index),ry,1);
            end
        end
    end
else
    Ypcv = zeros(nX,pY*nF);
    nseg = max(cvseg);
    Btemp = zeros(mX,pY*nF,nseg);
    r_ssq = zeros(nF,2,nseg);
    Ttemp = zeros(nX,nF,nseg);
    Ptemp = zeros(mX,nF,nseg);
    Qtemp = zeros(pY,nF,nseg);
    for a=1:nseg
        indexc = find(cvseg ~= a);
        indexp = find(cvseg == a);
        if isempty(indexp)
            error('ERROR: empty cv-segmen; index in "cvseg" must be consecutive sequence with whole numbers, starting with 1');
        end
        Xc = X(indexc,:);
        Yc = Y(indexc,:);
        Xp = X(indexp,:);
        Yp = Y(indexp,:);
        if prepro(1) == 1
            [Xc,mx] = meanc(Xc);
            Xp = meanc(Xp,mx);
        elseif prepro(1) == 2
            [Xc,mx,stdx] = autosc(Xc);
            Xp = autosc(Xp,mx,stdx);
        elseif prepro(1) == 3
            [Xc,rx] = rangesc(Xc);
            Xp = rangesc(Xp,rx);
        end
        if prepro(2) == 1
            [Yc,my] = meanc(Yc);
        elseif prepro(2) == 2
            [Yc,my,stdy] = autosc(Yc);
        elseif prepro(2) == 3
            [Yc,ry] = rangesc(Yc);
        end
        [Ttemp(indexc,:,a),Ptemp(:,:,a),W,U,Qtemp(:,:,a),Btemp(:,:,a),r_ssq(:,:,a)] = mypls(Xc,Yc,nF,options);
        for aa=1:nF
            index = (aa-1)*pY+1:aa*pY;
            if aa > 1
                Ypcv(indexp,index) = Ypcv(indexp,index-pY);
            end
            for aaa=1:length(indexp)
                xp = Xp(aaa,:);
                xpmv = isnan(xp);
                xp(xpmv) = 0;
                tp = xp*W(:,aa)/(W(~xpmv,aa)'*W(~xpmv,aa));
                Ypcv(indexp(aaa),index) = Ypcv(indexp(aaa),index) + tp*Qtemp(:,aa,a)';
                Xp(aaa,:) = Xp(aaa,:) - tp*Ptemp(:,aa,a)';
                Ttemp(indexp(aaa),aa) = tp;
            end
        end
        for aa=1:nF
            index = (aa-1)*pY+1:aa*pY;
            if prepro(2) == 1
                Ypcv(indexp,index) = meanc(Ypcv(indexp,index),my,1);
            elseif prepro(2) == 2
                Ypcv(indexp,index) = autosc(Ypcv(indexp,index),my,stdy,1);
            elseif prepro(2) == 3
                Ypcv(indexp,index) = rangesc(Ypcv(indexp,index),ry,1);
            end
        end
    end
end

Ymv = isnan(Y);
for a=1:nF
    index = (a-1)*pY+1:a*pY;
    Yp = Ypcv(:,index);
    rmspcv(1,a) = sqrt(mean((Yp(~Ymv)-Y(~Ymv)).^2));
end

if nargout == 3
    diagnos.nseg = nseg;
    g = (nseg - 1)/nseg;
    
    diagnos.s_rmspcv = zeros(size(rmspcv));
    if loo_cv
        for a=1:nX
            for aa=1:nF
                index = (aa-1)*pY+1:aa*pY;
                Yp = Ypcv(a,index);
                Yp = Yp(~Ymv(a,:));
                Yf = Y(~Ymv(a,:));
                diagnos.s_rmspcv(1,aa) = diagnos.s_rmspcv(1,aa) + (sqrt(mean((Yp-Yf).^2)) - rmspcv(1,aa))^2 * g;
            end
        end
    else
        for a=1:nseg
            indexp = find(cvseg == a);
            for aa=1:nF
                index = (aa-1)*pY+1:aa*pY;
                Yp = Ypcv(indexp,index);
                Yp = Yp(~Ymv(indexp,:));
                Yf = Y(indexp,:);
                Yf = Yf(~Ymv(indexp,:));
                diagnos.s_rmspcv(1,aa) = diagnos.s_rmspcv(1,aa) + (sqrt(mean((Yp(:)-Yf(:)).^2)) - rmspcv(1,aa))^2 * g;
            end
        end
    end
    diagnos.s_rmspcv = (diagnos.s_rmspcv).^0.5;
    
    diagnos.B = Bf;
    diagnos.s_B = zeros(size(Bf));
    for a=1:nseg
        diagnos.s_B = diagnos.s_B + (Bf-Btemp(:,:,a)).^2 * g;
    end
    diagnos.s_B = (diagnos.s_B).^0.5;
    
    diagnos.ssq = ssqf;
    diagnos.r_ssq(:,[1 2]) = min(r_ssq,[],3);
    diagnos.r_ssq(:,[3 4]) = max(r_ssq,[],3);
    
    diagnos.T = Tf;
    diagnos.s_T = zeros(size(Tf));
    diagnos.P = Pf;
    diagnos.s_P = zeros(size(Pf));
    diagnos.Q = Qf;
    diagnos.s_Q = zeros(size(Qf));
    for a=1:nseg
        C = pinv(Ttemp(:,:,a)'*Ttemp(:,:,a))*Ttemp(:,:,a)'*Tf;
        diagnos.s_T = diagnos.s_T + (Tf-Ttemp(:,:,a)*C).^2 * g;
        Cinv = pinv(C);
        diagnos.s_P = diagnos.s_P + (Pf-Ptemp(:,:,a)*Cinv).^2 * g;
        diagnos.s_Q = diagnos.s_Q + (Qf-Qtemp(:,:,a)*Cinv).^2 * g;
    end
    diagnos.s_T = (diagnos.s_T).^0.5;
    diagnos.s_P = (diagnos.s_P).^0.5;
    diagnos.s_Q = (diagnos.s_Q).^0.5;
end