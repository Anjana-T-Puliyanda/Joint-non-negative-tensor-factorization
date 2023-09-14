function [Zn,f] = blocknorm(Z,Znorm,Zin,options)
% function [Zn,f] = blocknorm(Z,Znorm,Zin,options)
% 040630 FvdB
% Rescales data table Z (possibly with index Zin) to make the sum-of-squares norm 
% of the table (or separate parts of the table) equal to target norm Znorm.
%
% in:
% Z (objects x variables) data table
% Znorm (1 x 1) or (1 x number of blocks) desired/target sum-of-squares norm for Z or blocks in Z
% Zin (cell array) structure index for variables per Z-block:
%    e.g. Zin = {1:50 51:100} is 2 blocks: Z1(1:50,:), Z2(51:100,:)
% options (1 x 2) convergence tolerance and maximum number of iterations (default = [1e-8 2000])
%
% out:
% Zn (objects x variables) normalized (sum-sum-squared of entries) data table
% f (1 x 1) or (1 x number of blocks) factor for desired normalization

if nargin < 2
    help blocknorm
    return
elseif nargin == 2
    if length(Znorm) ~= 1
        error('ERROR: for two input arguments "Znorm" must have length one');
    end
    Zin = {1:size(Z,2)};
elseif (nargin == 3) & length(Znorm) ~= length(Zin)
    error('ERROR: number of entries in "Znorm" and "Zin" must be the same');
end
if (nargin == 3)
    options(1) = 1e-8; % convergence tolerance
    options(2) = 2000; % maximum number of iterations
end

for a=1:length(Znorm)
    if Znorm(a) < 0
        s = ['ERROR: target sum-of-squares norm must be larger than zero (block #' num2str(a) ')'];
        error(s)
    end
    if abs(Znorm(a)) < eps
        s = ['WARNING: target sum-of-squares norm block #' num2str(a) ' set to "eps"'];
        disp(s)
        Znorm(a) = eps;
    end
end

for a=1:length(Znorm)
    Zt = Z(:,Zin{a});
    Zmv = sparse(isnan(Zt));
    Zt(Zmv) = 0;
    if Znorm(a) == 1
        f(a) = sqrt(sum(sum(Zt.^2)));
    else
        fmax = sum(sum(Zt.^2));
        fmaxn = sum(sum((Zt./fmax).^2));
        if fmaxn > Znorm(a)
            fmin = fmax;
            fminn = fmaxn;
            while fminn > Znorm(a)
                fmin = fmin*10;
                fminn = sum(sum((Zt./fmin).^2));
            end
        else
            fmin = fmax;
            fminn = fmaxn;
            while fmaxn < Znorm(a)
                fmax = fmax/10;
                fmaxn = sum(sum((Zt./fmax).^2));
            end
        end
        n = fmaxn;
        iter = 0;
        while (abs((Znorm(a)-n)/n) > options(1) & iter < options(2))
            iter = iter + 1;
            f(a) = (fmin + fmax)/2;
            n = sum(sum((Zt./f(a)).^2));
            if n > Znorm(a)
                fmax = f(a);
            else
                fmin = f(a);
            end
        end
        if iter == options(2)
            s = ['WARNING: maximum number of iterations (' num2str(options(2)) ') reached before convergence'];
            disp(s)
        end
    end
    Zt = Zt./f(a);
    Zt(Zmv) = NaN;
    Zn(:,Zin{a}) = Zt;
end