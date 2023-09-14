function [Zr,rz] = rangesc(Z,rz,rev)
% function [Zr,rz] = RangeSC(Z)
% function [Zr] = RangeSC(Z,rz)
% function [Zr] = RangeSC(Z,rz,1)
% 030630 FvdB
% Range-scales data table Z (new operation or based on know parameters),
% or reverses range-scaling operation.
% 
% in:
% Z (objects x variables) data block
% rz (2 x variables) column minimum and maximum
% rev (1 x 1) trigger for reversing range-scaling operation 
%
% out:
% Zr (objects x variables) '0 to 1' range scaled data block, or if 
%    input 'rev' is activated, reverses range block scaling
% rz (2 x variables) column minimum and maximum
%
% note:
% if a column in Z has entry-range smaller than 'eps', the scaled column in Zr is set to zero.

if (nargin == 0)
   help rangesc
   return
end   

[n,m] = size(Z);
if nargin == 1
   az = min(Z);
   bz = max(Z);
   rz = [az; bz];
   bz = bz-az;
   epsbz = find(bz < eps);
   lepsbz = length(epsbz);
   bz(epsbz) = ones(1,lepsbz);
   Zr = (Z-az(ones(n,1),:))./bz(ones(n,1),:);
   Zr(:,epsbz) = zeros(n,lepsbz);
elseif nargin == 2
   az = rz(1,:);
   bz = rz(2,:)-az;
   epsbz = find(bz < eps);
   lepsbz = length(epsbz);
   bz(epsbz) = ones(1,lepsbz);
   Zr = (Z-az(ones(n,1),:))./bz(ones(n,1),:);
   Zr(:,epsbz) = zeros(n,lepsbz);
elseif nargin == 3
   az = rz(1,:);
   bz = rz(2,:)-az;
   epsbz = find(bz == 0);
   lepsbz = length(epsbz);
   bz(epsbz) = ones(1,lepsbz);
   Zr = (Z.*bz(ones(n,1),:))+az(ones(n,1),:);
end
if lepsbz & (nargin < 3)
   s = ['WARNING: column(s) ' num2str(epsbz) ' is(are) set to zero due to lack of variance'];
   disp(s);
elseif lepsbz & (nargin == 3)
   s = ['WARNING: column(s) ' num2str(epsbz) ' is(are) set to mean-value(s) due to lack of variance'];
   disp(s);
end
