% FvdB 041004

clear all
close all

load MB_NIR_temp
nB = length(Xin); % number of blocks

%X = Xd; % "Xd" are Sag-Gol(window=5, order=1, deriv=1) first derivative spectra

figure % plot the five spectral blocks
plot(Xin{1},X(:,Xin{1}),'b',Xin{2},X(:,Xin{2}),'g',Xin{3},X(:,Xin{3}),'r',Xin{4},X(:,Xin{4}),'c',Xin{5},X(:,Xin{5}),'m');
blocklabels{1} = '30C'; blocklabels{2} = '40C'; blocklabels{3} = '50C'; blocklabels{4} = '60C'; blocklabels{5} = '70C';
grid; 
title('5 NIR blocks');
for a=1:nB
    text(mean(Xin{a}),0,blocklabels{a});
end

R = matrixcorr(X,Xin)

figure % plot the RV-coefficients
matrixblobs(R,blocklabels,[0 1 0]);

nF = 6; % number of factors

% MB = mbpca(Xd,2,X,Xin); % compute a multi-block PCA
MB = mbpca(Xd,nF,Xin,ones(nB,2)); % no user interaction, mean centering and equal weights of 1

figure % plot percentage explained
plot(1:nF,MB.ssq(:,1)*100,'o-k',1:nF,MB.ssq(:,2)*100,'b',1:nF,MB.ssq(:,3)*100,'g',1:nF,MB.ssq(:,4)*100,'r',1:nF,MB.ssq(:,nF)*100,'c',1:nF,MB.ssq(:,6)*100,'m')
axis([1 nF 75 100]);
grid; xlabel('factors'); ylabel('% expl.');

figure % plot block loaings
subplot(3,1,1);
plot(Xin{1},MB.Pb(Xin{1},1),'b',Xin{2},MB.Pb(Xin{2},1),'g',Xin{3},MB.Pb(Xin{3},1),'r',Xin{4},MB.Pb(Xin{4},1),'c',Xin{5},MB.Pb(Xin{5},1),'m');
grid; title('loading 1');
subplot(3,1,2);
plot(Xin{1},MB.Pb(Xin{1},2),'b',Xin{2},MB.Pb(Xin{2},2),'g',Xin{3},MB.Pb(Xin{3},2),'r',Xin{4},MB.Pb(Xin{4},2),'c',Xin{5},MB.Pb(Xin{5},2),'m');
grid; title('loading 2');
subplot(3,1,3);
plot(Xin{1},MB.Pb(Xin{1},3),'b',Xin{2},MB.Pb(Xin{2},3),'g',Xin{3},MB.Pb(Xin{3},3),'r',Xin{4},MB.Pb(Xin{4},3),'c',Xin{5},MB.Pb(Xin{5},3),'m');
grid; title('loading 3');
legend(blocklabels);

figure % plot super weights/loadings
subplot(2,2,1);
plot(MB.Wt(1,1),MB.Wt(1,2),'ob',MB.Wt(2,1),MB.Wt(2,2),'og',MB.Wt(3,1),MB.Wt(3,2),'or',MB.Wt(4,1),MB.Wt(4,2),'oc',MB.Wt(5,1),MB.Wt(5,2),'om');
grid; title('block weights'); xlabel('w 1'); ylabel('w 2');
legend(blocklabels);
subplot(2,2,2);
plot(MB.Wt(1,1),MB.Wt(1,3),'ob',MB.Wt(2,1),MB.Wt(2,3),'og',MB.Wt(3,1),MB.Wt(3,3),'or',MB.Wt(4,1),MB.Wt(4,3),'oc',MB.Wt(5,1),MB.Wt(5,3),'om');
grid; title('block weights'); xlabel('w 1'); ylabel('w 3');
subplot(2,2,3);
plot(MB.Wt(1,2),MB.Wt(1,3),'ob',MB.Wt(2,2),MB.Wt(2,3),'og',MB.Wt(3,2),MB.Wt(3,3),'or',MB.Wt(4,2),MB.Wt(4,3),'oc',MB.Wt(5,2),MB.Wt(5,3),'om');
grid; title('block weights'); xlabel('w 2'); ylabel('w 3');

figure % plot super and block scores
plot3(MB.Tt(:,1),MB.Tt(:,2),MB.Tt(:,3),'ok');
text(MB.Tt(:,1),MB.Tt(:,2),MB.Tt(:,3),num2str(Y(:,2),2));
grid; title('scores labels with water'); xlabel('t 1'); ylabel('t 2'); zlabel('t 3');
hold on;
plot3(MB.Tb(:,1,1),MB.Tb(:,2,1),MB.Tb(:,3,1),'.b'); % block 1 / 30C
plot3(MB.Tb(:,1,2),MB.Tb(:,2,2),MB.Tb(:,3,2),'.g'); % block 1 / 40C
plot3(MB.Tb(:,1,3),MB.Tb(:,2,3),MB.Tb(:,3,3),'.r'); % block 1 / 50C
plot3(MB.Tb(:,1,4),MB.Tb(:,2,4),MB.Tb(:,3,4),'.c'); % block 1 / 60C
plot3(MB.Tb(:,1,5),MB.Tb(:,2,5),MB.Tb(:,3,5),'.m'); % block 1 / 70C
hold off