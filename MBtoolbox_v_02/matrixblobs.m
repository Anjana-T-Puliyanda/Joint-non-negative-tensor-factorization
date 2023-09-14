function matrixblobs(X,Labels,options)
% function matrixblobs(X,Labels,options)
% 040701 FvdB
% Plots entries in matrix "X" as circles, where diameter equal the relative size.
% If X is symmetric and option 1 is active (e.g. a correlation matrix) only half the entries are plotted.
%
% in:
% X (n x m) data table
% Labels (cell array n if n=m or n + m if n~=m) Labels places on the plot
% options (1 x 1) 1 : symmetric correlation plotting, entries assumed in 0 to 1 or -1 to 1
%                     1 = yes, 0 = no (default)
%                 2 : plot 'blob' index; 1 = yes, 0 = no (default)
%                 3 : trigger 'fancy fill' color; 1 = yes, 0 = no (default)

if (nargin == 0)
    help matrixblobs
    return
elseif (nargin == 1)
    Labels = [];
    options = [0 0 0];
elseif (nargin == 2)
    options = [0 0 0];
end

[nX,mX] = size(X);
if isempty(Labels)
    for a=1:nX
        XLab{a} = num2str(a);        
    end
    for a=1:mX
        YLab{a} = num2str(a);        
    end
elseif (nX==mX) & (length(Labels) == nX)
    XLab = Labels;
    YLab = Labels;
elseif (length(Labels) == nX+mX)
    XLab = Labels(1:nX);
    YLab = Labels(nX+1:nX+mX);
else
    error('ERROR: number of entries in "Labels" must be "n" or "n+m"');
end

t = (0:.01:1)'*2*pi;
x = sin(t);
y = cos(t);
minX = min(min(abs(X)));
maxX = max(max(abs(X)));

if options(3)
    color = t;
    mincolor = -t;
else
    color = 'b';
    mincolor = 'c';
end

if options(1)
    if (nX ~= mX) | (minX < -1) | (maxX > 1)
        error('ERROR: for option(1) the matrix "X" must be square and entries must be between -1 and 1')
    end
    axis([0.5 mX+0.5 0.5 nX+0.5])
    R = (abs(X)-0)./(1-0)*0.47 + 0.02;
    R = R.*sign(X);
    hold on
    for a=1:nX
        for aa=1:a
            if R(a,aa) > 0
                fill(x*R(a,aa)+aa,y*R(a,aa)+(nX-a+1),color)
            else
                fill(x*R(a,aa)+aa,y*R(a,aa)+(nX-a+1),mincolor)
            end
        end
    end
    if options(2) 
        fill(x*0.02+mX-1,y*0.02+nX,color)
        fill(x*0.47+mX,y*0.47+nX,color)
        text(mX-1,nX,num2str(0,'%0.2f'))
        text(mX,nX,num2str(1,'%0.2f'))
        if any(R(:)<0)
            fill(x*-0.47+mX-2,y*-0.47+nX,mincolor)
            text(mX-2,nX,num2str(-1,'%0.2f'))
        end
    end
    axis equal
else
    axis([0.5 mX+0.5 0.5 nX+0.5])
    R = (abs(X)-minX)./(maxX-minX)*0.47 + 0.02;
    R = R.*sign(X);
    hold on
    for a=1:nX
        for aa=1:mX
            if R(a,aa) > 0
                fill(x*R(a,aa)+aa,y*R(a,aa)+(nX-a+1),color)
            else
                fill(x*R(a,aa)+aa,y*R(a,aa)+(nX-a+1),mincolor)
            end
        end
    end
    if options(2) 
        axis([0.5 mX+0.5 0.5 nX+2.5])
        minX = min(min(X));
        maxX = max(max(X));
        [mina,minaa] = find(X == minX);
        mina = mina(1);
        minaa = minaa(1);
        [maxa,maxaa] = find(X == maxX);
        maxa = maxa(1);
        maxaa = maxaa(1);
        if R(mina,minaa) > 0
            fill(x*R(mina,minaa)+mX-2,y*R(mina,minaa)+(nX+2),color)
        else
            fill(x*R(mina,minaa)+mX-2,y*R(mina,minaa)+(nX+2),mincolor)
        end    
        if R(maxa,maxaa) > 0
            fill(x*R(maxa,maxaa)+mX,y*R(maxa,maxaa)+(nX+2),color)
        else
            fill(x*R(maxa,maxaa)+mX,y*R(maxa,maxaa)+(nX+2),mincolor)
        end
        text(mX-2,nX+2,num2str(minX,'%0.2f'))
        text(mX,nX+2,num2str(maxX,'%0.2f'))
        minabsX = min(min(abs(X)));
        if minX ~= minabsX
            [minabsa,minabsaa] = find(abs(X) == minabsX);
            minabsa = minabsa(1);
            minabsaa = minabsaa(1);
            if R(minabsa,minabsaa) > 0
                fill(x*R(minabsa,minabsaa)+mX-1,y*R(minabsa,minabsaa)+(nX+2),color)
            else
                fill(x*R(minabsa,minabsaa)+mX-1,y*R(minabsa,minabsaa)+(nX+2),mincolor)
            end
            text(mX-1,nX+2,num2str(X(minabsa,minabsaa),'%0.2f'))            
        end
    end
end
axis off
axis equal
for a=1:nX
    text(0.2,nX-a+1,XLab{a})
end
for a=1:mX
    text(a,0.2,YLab{a})
end    
hold off