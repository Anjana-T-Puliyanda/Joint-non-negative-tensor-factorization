function [X,Xin2,Xin4,XVarLab,Y,Yin2,YVarLab,ObjLab,T,P] = MBdata(N)
% function [X,Xin2,Xin4,XVarLab,Y,Yin2,YVarLab,ObjLab,T,P] = MBdata(N)
% 040630 FvdB
% Generate "toy" data for the multi-block toolbox
%
% in:
% N (1 x 1) number of objects/samples in the data set (default = 10)
%
% out:
% X (N x 100) predictor data table
% Xin2 (cell array) index into "X" assuming two blocks (2 x 50 variables)
% Xin4 (cell array) index into "X" assuming four blocks (4 x 25 variables)
% XVarLab (1 x 100) labels for all "X" variables
% Y (N x 3) response data table
% Yin2 (cell array) index into "Y" assuming two blocks (1 + 2 variables)
% YVarLab (1 x 3) labels for all "Y" variables
% Objlab (N x 1) labels for all objects/samples
% T (N x 4) original score values for "X"
% P (40 x 4) original loading vectors for "X"

if (nargin == 0)
    N = 10;
end

T = rand(N,4)*.9+.1;
P(1,:) = 0.6*exp(-0.5*(([1:50]-15)/5.5).^2)./(15*(5.5*pi)^0.5);
P(2,:) = 0.6*exp(-0.5*(([1:50]-35)/7.5).^2)./(35*(7.5*pi)^0.5);
P(3,:) = 1.1*exp(-0.5*(([1:50]-20)/4).^2)./(20*(4*pi)^0.5);
P(4,:) = 2.1*exp(-0.5*(([1:50]-40)/2.5).^2)./(40*(2.5*pi)^0.5);
X = [T(:,1:2)*P(1:2,:) T(:,3:4)*P(3:4,:)];
X = X + randn(size(X))*1e-4;
Xin2{1} = 1:50; Xin2{2} = 51:100;
Xin4{1} = 1:25; Xin4{2} = 26:50; Xin4{3} = 51:75; Xin4{4} = 76:100;
XVarLab = 0:.1:9.9;
Y(:,1) = 0.7*T(:,1) + 0.5*T(:,4) + randn(size(T(:,1)))*3e-2;
Y(:,2) = 0.4*T(:,1).^0.5 + 0.6*T(:,2) + randn(size(T(:,1)))*4e-2;
Y(:,3) = T(:,1) + randn(size(T(:,1)))*5e-2;
Yin2{1} = 1; Yin2{2} = 2:3;
YVarLab = [0.2 0.8 1.0];

a = 97; b = 1;
for c=1:N
    ObjLab(c,:) = [char(a) num2str(b)];
    a = a + 1;
    if (a > 122)
        a = 97;
        b = b + 1;
    end
end