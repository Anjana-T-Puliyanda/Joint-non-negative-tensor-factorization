function [Xn,Yn,diagnos] = procrus2D(X,Y,options)
% function [Xn,Yn,diagnos] = procrus2D(X,Y,options)
% 030805 FvdB
% Procrustes (full) matching for two 'images' X and Y.
% Order: 1) center 2) resize 3) rotate 4) resize
%
% in:
% X (n x 2) coordinates of landmarks in X-image (the target image)
% Y (n x 2) coordinates of landmarks in Y-image
% options (1 x 2) default [0 0]
%       1: try to change column signs (after procrusted rotation) to maximize simularity (default = 0 is no)
%       2: k=(1 x 1) first rotate object "k" towards coordinate (x,y,z)=(1,1,1)
%          (default k=0 for no rotation; "k" must be in the range of "n")
%
% out:
% Xn (n x 2) new X-image (centered and 'center-scaled')
% Yn (n x 2) new Y-image (centered, 'center-scaled', rotated and 'match-scaled')
% diagnos (struct) Procrustes step-diagnostics
%
% based on: Robinson et al. 'Planar Procrustes analysis of tooth shape' Archives of Oral Biology 46(2001)191-199 (ERRORS!!!)

if nargin < 2
    help procrus2D
    return
elseif nargin == 2
    options = [0 0];
end

[nX,mX] = size(X);
[nY,mY] = size(Y);

if nX ~= nY
    s = ['ERROR: number of landmarks in X ' num2str(nX) ' and Y ' num2str(nY) ' must be the same'];
    error(s)
end
if (mX ~= 2) | (mY ~= 2)
    error('ERROR: number of coordinates per X- and Y-landmark must be two')
end
if any(isnan([X Y]))
    error('ERROR: function does not work with missing values (NaN)');
end
if options(2)
    if (options(2) < 0) | (options(2) > nX)
        error('ERROR: options(2) must be an object in "X"')
    end
end

diagnos.cX = mean(X);
diagnos.cY = mean(Y);
diagnos.sXpre = sum(sum((X-ones(nX,1)*diagnos.cX).^2))^.5;
diagnos.sYpre = sum(sum((Y-ones(nY,1)*diagnos.cY).^2))^.5;
Xn = (X-ones(nX,1)*diagnos.cX)/diagnos.sXpre;
Yn = (Y-ones(nY,1)*diagnos.cY)/diagnos.sYpre;
diagnos.rX = 0;

if options(2)
    k = options(2);
    rX = atan(sum(Xn(k,1)-Xn(k,2))/sum(Xn(k,1)+Xn(k,2)));
    norms(1) = norm([(Xn(k,1)*cos(rX)-Xn(k,2)*sin(rX)) (Xn(k,1)*sin(rX)+Xn(k,2)*cos(rX))] - [1 1]);
    norms(2) = norm([(Xn(k,1)*cos(rX+pi)-Xn(k,2)*sin(rX+pi)) (Xn(k,1)*sin(rX+pi)+Xn(k,2)*cos(rX+pi))] - [1 1]);
    if norms(2) < norms(1)
        diagnos.rX = rX + pi;
    else
        diagnos.rX = rX;
    end
    Xn = [(Xn(:,1)*cos(diagnos.rX)-Xn(:,2)*sin(diagnos.rX)) (Xn(:,1)*sin(diagnos.rX)+Xn(:,2)*cos(diagnos.rX))];
end

diagnos.signcol1 = 0;
diagnos.signcol2 = 0;
if options(1)
    % No sign flipping
    Yn_t(:,1) = Yn(:,1); Yn_t(:,2) = Yn(:,2);
    rY{1} = atan(sum(Yn_t(:,1).*Xn(:,2)-Yn_t(:,2).*Xn(:,1))/sum(Yn_t(:,1).*Xn(:,1)+Yn_t(:,2).*Xn(:,2)));
    rsY{1} = (sum(Yn_t(:,1).*Xn(:,1)+Yn_t(:,2).*Xn(:,2))^2 + sum(Yn_t(:,1).*Xn(:,2)+Yn_t(:,2).*Xn(:,1))^2)^0.5;
    Yn_t = [rsY{1}*(Yn_t(:,1)*cos(rY{1})-Yn_t(:,2)*sin(rY{1})) rsY{1}*(Yn_t(:,1)*sin(rY{1})+Yn_t(:,2)*cos(rY{1}))];
    sYpost{1} = sum(sum(Yn_t.^2))^.5;
    Yn_t = Yn_t/sYpost{1};
    norms(1) = norm(Xn-Yn);
    
    % Flip first coordinate
    Yn_t(:,1) = -Yn(:,1); Yn_t(:,2) = Yn(:,2);
    rY{2} = atan(sum(Yn_t(:,1).*Xn(:,2)-Yn_t(:,2).*Xn(:,1))/sum(Yn_t(:,1).*Xn(:,1)+Yn_t(:,2).*Xn(:,2)));
    rsY{2} = (sum(Yn_t(:,1).*Xn(:,1)+Yn_t(:,2).*Xn(:,2))^2 + sum(Yn_t(:,1).*Xn(:,2)+Yn_t(:,2).*Xn(:,1))^2)^0.5;
    Yn_t = [rsY{2}*(Yn_t(:,1)*cos(rY{2})-Yn_t(:,2)*sin(rY{2})) rsY{2}*(Yn_t(:,1)*sin(rY{2})+Yn_t(:,2)*cos(rY{2}))];
    sYpost{2} = sum(sum(Yn_t.^2))^.5;
    Yn_t = Yn_t/sYpost{2};
    norms(2) = norm(Xn-Yn_t);
    
    % Flip second coordinate
    Yn_t(:,1) = Yn(:,1); Yn_t(:,2) = -Yn(:,2);
    rY{3} = atan(sum(Yn_t(:,1).*Xn(:,2)-Yn_t(:,2).*Xn(:,1))/sum(Yn_t(:,1).*Xn(:,1)+Yn_t(:,2).*Xn(:,2)));
    rsY{3} = (sum(Yn_t(:,1).*Xn(:,1)+Yn_t(:,2).*Xn(:,2))^2 + sum(Yn_t(:,1).*Xn(:,2)+Yn_t(:,2).*Xn(:,1))^2)^0.5;
    Yn_t = [rsY{3}*(Yn_t(:,1)*cos(rY{3})-Yn_t(:,2)*sin(rY{3})) rsY{3}*(Yn_t(:,1)*sin(rY{3})+Yn_t(:,2)*cos(rY{3}))];
    sYpost{3} = sum(sum(Yn_t.^2))^.5;
    Yn_t = Yn_t/sYpost{3};
    norms(3) = norm(Xn-Yn_t);
    
    % Flip both coordinates
    Yn_t(:,1) = -Yn(:,1); Yn_t(:,2) = -Yn(:,2);
    rY{4} = atan(sum(Yn_t(:,1).*Xn(:,2)-Yn_t(:,2).*Xn(:,1))/sum(Yn_t(:,1).*Xn(:,1)+Yn_t(:,2).*Xn(:,2)));
    rsY{4} = (sum(Yn_t(:,1).*Xn(:,1)+Yn_t(:,2).*Xn(:,2))^2 + sum(Yn_t(:,1).*Xn(:,2)+Yn_t(:,2).*Xn(:,1))^2)^0.5;
    Yn_t = [rsY{4}*(Yn_t(:,1)*cos(rY{4})-Yn_t(:,2)*sin(rY{4})) rsY{4}*(Yn_t(:,1)*sin(rY{4})+Yn_t(:,2)*cos(rY{4}))];
    sYpost{4} = sum(sum(Yn_t.^2))^.5;
    Yn_t = Yn_t/sYpost{4};
    norms(4) = norm(Xn-Yn_t);
    
    [dummy,minnorm] = min(norms);
    switch minnorm
        case 2
            diagnos.signcol1 = 1;
            Yn(:,1) = -Yn(:,1);
        case 3
            diagnos.signcol2 = 1;
            Yn(:,2) = -Yn(:,2);
        case 4
            diagnos.signcol1 = 1;
            Yn(:,1) = -Yn(:,1);
            diagnos.signcol2 = 1;
            Yn(:,2) = -Yn(:,2);
    end
    diagnos.rY = rY{minnorm};
    diagnos.rsY = rsY{minnorm};
    diagnos.sYpost = sYpost{minnorm};
    Yn = [diagnos.rsY*(Yn(:,1)*cos(diagnos.rY)-Yn(:,2)*sin(diagnos.rY)) diagnos.rsY*(Yn(:,1)*sin(diagnos.rY)+Yn(:,2)*cos(diagnos.rY))];
    Yn = Yn/diagnos.sYpost;
else
    diagnos.rY = atan(sum(Yn(:,1).*Xn(:,2)-Yn(:,2).*Xn(:,1))/sum(Yn(:,1).*Xn(:,1)+Yn(:,2).*Xn(:,2)));
    diagnos.rsY = (sum(Yn(:,1).*Xn(:,1)+Yn(:,2).*Xn(:,2))^2 + sum(Yn(:,1).*Xn(:,2)+Yn(:,2).*Xn(:,1))^2)^0.5;
    Yn = [diagnos.rsY*(Yn(:,1)*cos(diagnos.rY)-Yn(:,2)*sin(diagnos.rY)) diagnos.rsY*(Yn(:,1)*sin(diagnos.rY)+Yn(:,2)*cos(diagnos.rY))];
    diagnos.sYpost = sum(sum(Yn.^2))^.5;
    Yn = Yn/diagnos.sYpost;
end
