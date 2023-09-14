function [Za,mz,stdz] = autosc(Z,mz,stdz,rev)
% function [Za,mz,stdz] = autosc(Z)
% function [Za] = autosc(Z,mz,stdz)
% function [Za] = autosc(Z,mz,stdz,1)
% 040630 FvdB
% Auto-scales data table Z (new operation or based on know parameters),
% or reverses auto-scaling operation.
%
% in: 
% Z (objects x variables) data block
% mz (1 x variables) column means
% stdz (1 x variables) column standard deviations
% rev (1 x 1) trigger for reversing auto scale operation
%
% out:
% Za (objects x variables) auto scaled data block, or if
%    input 'rev' is activated, reversed auto-scaled data block
% mz (1 x variables) column means
% stdz (1 x variables) column standard deviations
%
% note:
% if a column in Z has standard deviation smaller than 'eps', 
% the scaled column in Za is set to zero.

if (nargin == 0) | (nargin == 2)
    help autosc
    return
end   

[nZ,mZ] = size(Z);
if nargin == 1
    if ~sum(sum(isnan(Z)))
        mz = mean(Z);
        stdz = std(Z);
        epsstdz = find(stdz < eps);
        lepsstdz = length(epsstdz);
        stdz(epsstdz) = ones(1,lepsstdz);
        Za = (Z-mz(ones(nZ,1),:))./stdz(ones(nZ,1),:);
        Za(:,epsstdz) = zeros(nZ,lepsstdz);
        stdz(epsstdz) = zeros(1,lepsstdz);
    else
        Zmv = sparse(isnan(Z));
        zmv = nZ./(nZ-sum(Zmv));
        Z(Zmv) = 0;
        mz = mean(Z).*zmv;
        for aa=1:mZ
            stdz(aa) = std(Z(~Zmv(:,aa),aa));
        end
        epsstdz = find(stdz < eps);
        lepsstdz = length(epsstdz);
        stdz(epsstdz) = ones(1,lepsstdz);
        Za = (Z-mz(ones(nZ,1),:))./stdz(ones(nZ,1),:);
        Za(:,epsstdz) = zeros(nZ,lepsstdz);
        stdz(epsstdz) = zeros(1,lepsstdz);
        Za(Zmv) = NaN;
    end
elseif nargin == 3
    epsstdz = find(stdz == 0);
    lepsstdz = length(epsstdz);
    stdz(epsstdz) = ones(1,lepsstdz);
    Za = (Z-mz(ones(nZ,1),:))./stdz(ones(nZ,1),:);
elseif nargin == 4
    epsstdz = find(stdz == 0);
    lepsstdz = length(epsstdz);
    stdz(epsstdz) = ones(1,lepsstdz);
    Za = (Z.*stdz(ones(nZ,1),:))+mz(ones(nZ,1),:);
end
if lepsstdz & (nargin < 4)
    s = ['WARNING: column(s) ' num2str(epsstdz) ' is(are) set to zeros due to lack of variance'];
    disp(s);
elseif lepsstdz & (nargin == 4)
    s = ['WARNING: column(s) ' num2str(epsstdz) ' is(are) set to mean-value(s) due to lack of variance'];
    disp(s);
end