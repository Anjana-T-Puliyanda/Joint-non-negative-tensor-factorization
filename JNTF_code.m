set(groot,'defaultLineLineWidth',2)

clear all; close all;clc;

%% add path to these toolbox imports according to your system 
addpath('G:\tensor_toolbox-master');
addpath('G:\L-BFGS-B-C-master\Matlab');
addpath('G:\poblano_toolbox_1.1');


addpath('G:\MBtoolbox_v_02')
addpath('G:\nway331')

%% load the data
filename_ftir = 'bitumenfinaldata.xlsx';
D_ftir= xlsread(filename_ftir);  %data matrix
T_ftir=D_ftir(1,2:end); %temperature
time_ftir=D_ftir(2,2:end); %time
X_ftir=D_ftir(3:end,:);%intensity
lam_ftir=X_ftir(:,1);%wavenumber
X11_ftir=X_ftir(:,2:end); %intensity
D_ftir=2-log10(X11_ftir);% change to absorbance which is interpreted as conc

filename_hnmr='hnmr_finaldata.xlsx';
D_hnmr= xlsread(filename_hnmr);  %data matrix
T_hnmr=D_hnmr(1,2:end); %temperature
time_hnmr=D_hnmr(2,2:end); %time
X_hnmr=D_hnmr(3:end,:);%intensity
conc_hnmr=X_hnmr(:,1);%wavenumber
X11_hnmr=X_hnmr(:,2:end); %intensity
D_hnmr=X11_hnmr;


%baseline and background correction of the data
Dback_ftir=msbackadj(lam_ftir,D_ftir);
Data_ftir=mssgolay(lam_ftir,Dback_ftir,'DEGREE',2,'Span',5);
Noise_ftir=Dback_ftir-Data_ftir;

Dback_hnmr=msbackadj(conc_hnmr,D_hnmr);
Data_hnmr=mssgolay(conc_hnmr,Dback_hnmr,'DEGREE',2,'Span',5);
Noise_hnmr=Dback_hnmr-Data_hnmr;



%% create the data tensors
global Z_ftir; global Z_hnmr; global P_ftir; global P_hnmr;
Z_ftir=NaN(8,23,1764); % 8 temperatures, 23 res times, 1764 sspectral channels
Z_hnmr=NaN(8,23,683); % 8 temperatures, 23 res times, 683 spectral channels

% now start populating the tensor blocks with the individual spectra

%spectra at 20 degrees
Z_ftir(1,1,:)=Data_ftir(:,1)';
%spectra at 150
Z_ftir(2,7,:)=Data_ftir(:,2)';%2:9
Z_ftir(2,10,:)=Data_ftir(:,3)';%2:9
Z_ftir(2,12:2:22,:)=Data_ftir(:,4:9)';%2:9
%spectra at 200
Z_ftir(3,7,:)=Data_ftir(:,10)';%10:15
Z_ftir(3,10:2:16,:)=Data_ftir(:,11:14)';%10:15
Z_ftir(3,22,:)=Data_ftir(:,15)';%10:15
%spectra at 250
Z_ftir(4,14,:)=Data_ftir(:,16)';
%spectra at 300
Z_ftir(5,10:2:22,:)=Data_ftir(:,17:23)';%17:23
%spectra at 340
Z_ftir(6,2,:)=Data_ftir(:,24)';%24:28
Z_ftir(6,7,:)=Data_ftir(:,25)';%24:28
Z_ftir(6,10,:)=Data_ftir(:,26)';%24:28
Z_ftir(6,14,:)=Data_ftir(:,27)';%24:28
Z_ftir(6,22,:)=Data_ftir(:,28)';%24:28
%spectra at 360
Z_ftir(7,2:5,:)=Data_ftir(:,29:32)';%29:35
Z_ftir(7,7,:)=Data_ftir(:,33)';%29:35
Z_ftir(7,14,:)=Data_ftir(:,34)';%29:35
Z_ftir(7,23,:)=Data_ftir(:,35)';%29:35
%spectra at 400
Z_ftir(8,2:5,:)=Data_ftir(:,36:39)';
Z_ftir(8,7,:)=Data_ftir(:,40)';
Z_ftir(8,8,:)=Data_ftir(:,41)';
Z_ftir(8,10,:)=Data_ftir(:,42)';
Z_f=Z_ftir;
Z_ftir(isnan(Z_ftir))=0;

P_ftir=(Z_ftir~=0);
P_ftir=tensor(P_ftir);
    
Z_ftir=tensor(Z_ftir);
%spectra at 150
Z_hnmr(2,6,:)=Data_hnmr(:,1)';%1:8
Z_hnmr(2,9:2:21,:)=Data_hnmr(:,2:8)';%1:8
%spectra at 200
Z_hnmr(3,6,:)=Data_hnmr(:,9)';%9:16
Z_hnmr(3,9:2:21,:)=Data_hnmr(:,10:16)';%9:16
%spectra at 250
Z_hnmr(4,6,:)=Data_hnmr(:,17)';%17:24
Z_hnmr(4,9:2:21,:)=Data_hnmr(:,18:24)';%17:24
%spectra at 300
Z_hnmr(5,6,:)=Data_hnmr(:,25)';%25:32
Z_hnmr(5,9:2:21,:)=Data_hnmr(:,26:32)';%25:32
Z_h=Z_hnmr;
Z_hnmr(isnan(Z_hnmr))=0;

P_hnmr=(Z_hnmr~=0);
P_hnmr=tensor(P_hnmr);
Z_hnmr=tensor(Z_hnmr);



% Use corcondia to det rank of FTIR
LOF_F=[];CF_F=[];
for i=1:8
    [Factors_ftir,it_f,lof_f,cf_f]=parafac(Z_f,i);
    %M_f = nmodel(Factors_ftir);
    LOF_F=[LOF_F;lof_f];
    %cf_f = corcond(Z_f,Factors_ftir);
    CF_F=[CF_F;cf_f];
end

figure()
subplot(1,2,1)
plot(1:length(LOF_F),LOF_F,'-BX')
axis tight
xlabel('Number of components','fontweight','bold','FontSize',20)
ylabel('Lack of fit (LOF)','fontweight','bold','FontSize',20)
set(gca,'FontSize',20,'fontweight','bold')
subplot(1,2,2)
plot(1:length(CF_F),CF_F,'-BX')
axis tight
xlabel('Number of components','fontweight','bold','FontSize',20)
ylabel('Core consistency','fontweight','bold','FontSize',20)
set(gca,'FontSize',20,'fontweight','bold')
[A_f,B_f,C_f]=fac2let(Factors_ftir);
M_f = nmodel(Factors_ftir);

R_f=4; % rank of the FTIR matrix from core consistency

for i=1:(length(CF_F)-1)

    if abs((CF_F(i+1)-CF_F(i))/CF_F(i))>10e-3
        R_f=i
        break;
    end
end

%% Use corcondia to det rank of HNMR
LOF_H=[];CF_H=[];
for i=1:8
    [Factors_hnmr,it_h,lof_h,cf_h]=parafac(Z_h,i);
    %M_f = nmodel(Factors_ftir);
    LOF_H=[LOF_H;lof_h];
    %cf_h = corcond(Z_h,Factors_hnmr);
    CF_H=[CF_H;cf_h];
end

figure()
subplot(1,2,1)
plot(1:length(LOF_H),LOF_H,'-BX')
axis tight
xlabel('Number of components','fontweight','bold','FontSize',20)
ylabel('Lack of fit (LOF)','fontweight','bold','FontSize',20)
set(gca,'FontSize',20,'fontweight','bold')
subplot(1,2,2)
plot(1:length(CF_H),CF_H,'-BX')
axis tight
xlabel('Number of components','fontweight','bold','FontSize',20)
ylabel('Core consistency','fontweight','bold','FontSize',20)
set(gca,'FontSize',20,'fontweight','bold')

R_h=4; %rank of hnmr matrix from core consistency


%Define rank for joint factorization
R=min(R_h,R_f);



%% Fit a PARAFAC model for Z_ftir

% Set up optimization parameters
% Get the defaults
ncg_opts = ncg('defaults');
% Tighten the stop tolerance (norm of gradient). This is often too large.
ncg_opts.StopTol = 1.0e-6;
% Tighten relative change in function value tolearnce. This is often too large.
ncg_opts.RelFuncTol = 1.0e-20;
% Increase the number of iterations.
ncg_opts.MaxIters = 10^4;
% Only display every 10th iteration
ncg_opts.DisplayIters = 10;
% Display the final set of options
ncg_opts

% Initial guess can either be random or done via nonnegative double
% singular value decomposition
%[Ainit_ftir,~]=NNDSVD(double(tenmat(Z_ftir,1)*tenmat(Z_ftir,1)'),R_f,[]);
%[Binit_ftir,~]=NNDSVD(double(tenmat(Z_ftir,2)*tenmat(Z_ftir,2)'),R_f,[]);
%Cinit_ftir=double(tenmat(Z_ftir,3)*khatrirao(Binit_ftir,Ainit_ftir))*pinv((Binit_ftir'*Binit_ftir).*(Ainit_ftir'*Ainit_ftir));
%Minit_ftir={Ainit_ftir;Binit_ftir;Cinit_ftir};
%Minit_ftir = create_guess('Data',Z_ftir, 'Num_Factors', R_f, ...
  %    'Factor_Generator', 'nvecs');
%Call optimizer

%but below I choose to run with 100 random initial guesses and take the
%solution that results in least SSE.

SSE_ftir=[];M_FTIR = struct('tensrs',[]);

for i=1:100
[M_ftir,~,output_ftir] = cp_wopt(Z_ftir, P_ftir, R_f,'lower',0);%,'init', Minit_ftir);
SSE_ftir=[SSE_ftir;sum(sum(sum((double(tensor(M_ftir)-Z_ftir)).^2)))];
M_FTIR(i).tensrs=M_ftir;
end

csvwrite('SSE_ftir.csv',SSE_ftir);

[sse_minf,sse_locf]=min(SSE_ftir)

M_ftir=M_FTIR(sse_locf).tensrs;


csvwrite('M_ftirU1.csv',M_ftir.U{1});
csvwrite('M_ftirU2.csv',M_ftir.U{2});
csvwrite('M_ftirU3.csv',M_ftir.U{3});

%convert to a ttensor
R1 = length(M_ftir.lambda);  %<-- Number of factors in X.
core1 = tendiag(M_ftir.lambda, repmat(R1,1,ndims(M_ftir))); %<-- Create a diagonal core.
Y_ftir = ttensor(core1, M_ftir.U) ;%<-- Assemble the ttensor.

%% Fit PARAFAC for Z_hnmr
%Minit_hnmr = create_guess('Data',Z_hnmr, 'Num_Factors', R_h, ...
 %    'Factor_Generator', 'nvecs');
%[Ainit_hnmr,~]=NNDSVD(double(tenmat(Z_hnmr,1)*tenmat(Z_hnmr,1)'),R_h,[]);
%[Binit_hnmr,~]=NNDSVD(double(tenmat(Z_hnmr,2)*tenmat(Z_hnmr,2)'),R_h,[]);
%Cinit_hnmr=double(tenmat(Z_hnmr,3)*khatrirao(Binit_hnmr,Ainit_hnmr))*pinv((Binit_hnmr'*Binit_hnmr).*(Ainit_hnmr'*Ainit_hnmr));
%Minit_hnmr={Ainit_hnmr;Binit_hnmr;Cinit_hnmr};

SSE_hnmr=[];M_HNMR = struct('tensrs',[]);

for i=1:100
[M_hnmr,~,output_hnmr] = cp_wopt(Z_hnmr, P_hnmr, R_h,'lower',0);%,'init', Minit_ftir);
SSE_hnmr=[SSE_hnmr;sum(sum(sum((double(tensor(M_hnmr)-Z_hnmr)).^2)))];
M_HNMR(i).tensrs=M_hnmr;
end

csvwrite('SSE_hnmr.csv',SSE_hnmr);

[sse_minh,sse_loch]=min(SSE_hnmr)

M_hnmr=M_HNMR(sse_loch).tensrs;


csvwrite('M_hnmrU1.csv',M_hnmr.U{1});
csvwrite('M_hnmrU2.csv',M_hnmr.U{2});
csvwrite('M_hnmrU3.csv',M_hnmr.U{3});

%convert to a ttensor
R2 = length(M_hnmr.lambda);  %<-- Number of factors in X.
core2 = tendiag(M_hnmr.lambda, repmat(R2,1,ndims(M_hnmr))); %<-- Create a diagonal core.
Y_hnmr = ttensor(core2, M_hnmr.U) ;%<-- Assemble the ttensor.


%% Joint decompositon of FTIR and HNMR data tensors
%defining the objective function
%F=@(A,B,H1,H2)norm(P_ftir.*(Z_ftir-ktensor({A,B,H1})))+norm(P_hnmr.*(Z_hnmr-ktensor({A,B,H2})));

%initialize decision variables
A=(M_ftir.U{1}+M_hnmr.U{1})./2;
B=(M_ftir.U{2}+M_hnmr.U{2})./2;
H1=M_ftir.U{3};
H2=M_hnmr.U{3};

%%Defining grad of obj fun in case you want to implement your own optimizer
%  gradA=@(A,B,H1,H2) double(2*(tenmat((P_ftir.*ktensor({A,B,H1})),1)-tenmat((P_ftir.*Z_ftir),1))*khatrirao(H1,B))+double(2*(tenmat((P_hnmr.*ktensor({A,B,H2})),1)-tenmat((P_hnmr.*Z_hnmr),1))*khatrirao(H2,B));
%  gradB=@(A,B,H1,H2) double(2*(tenmat((P_ftir.*ktensor({A,B,H1})),2)-tenmat((P_ftir.*Z_ftir),2))*khatrirao(H1,A))+double(2*(tenmat((P_hnmr.*ktensor({A,B,H2})),2)-tenmat((P_hnmr.*Z_hnmr),2))*khatrirao(H2,A));
%  gradH1=@(A,B,H1,H2)double(2*(tenmat((P_ftir.*ktensor({A,B,H1})),3)-tenmat((P_ftir.*Z_ftir),3))*khatrirao(B,A));
%  gradH2=@(A,B,H1,H2)double(2*(tenmat((P_hnmr.*ktensor({A,B,H2})),3)-tenmat((P_hnmr.*Z_hnmr),3))*khatrirao(B,A));
% % 
% %Defining stepsizes
% stepA=@(A,B,H1,H2)A./(double(-tenmat(Z_ftir,1)*khatrirao(H1,B))+double(-tenmat(Z_hnmr,1)*khatrirao(H2,B)));
% stepB=@(A,B,H1,H2)B./(double(-tenmat(Z_ftir,2)*khatrirao(H1,A))+double(-tenmat(Z_hnmr,2)*khatrirao(H2,A)));
% stepH1=@(A,B,H1,H2)H1./(double(-tenmat(Z_ftir,3)*khatrirao(B,A)));
% stepH2=@(A,B,H1,H2)H2./(double(-tenmat(Z_hnmr,3)*khatrirao(B,A)));

% maxiter=3000;
% F0=F(A,B,H1,H2); Ft=F0;
% for i=1:maxiter
%     
%     A=double(tenmat(Z_ftir,1)*khatrirao(H1,B)+tenmat(Z_hnmr,1)*khatrirao(H2,B))*pinv(((H1'*H1).*(B'*B))+((H1'*H1).*(B'*B)));
%     B=double(tenmat(Z_ftir,2)*khatrirao(H1,A)+tenmat(Z_hnmr,2)*khatrirao(H2,A))*pinv(((H1'*H1).*(A'*A))+((H2'*H2).*(A'*A)));   
%     H1=double(tenmat(Z_ftir,3)*khatrirao(B,A))*pinv((B'*B).*(A'*A));
%     H2=double(tenmat(Z_hnmr,3)*khatrirao(B,A))*pinv((B'*B).*(A'*A));
%     Ft1=F(A,B,H1,H2)
%     ERE=abs((Ft-Ft1)/(F0-Ft1))
%     if ERE<=10^-6||Ft1>Ft
%         break;
%     end
%     Ft=Ft1;
%     
% end 
% 




%% Using LBFGSB optimizer to solve the joint decomposition
opts.x0=x0;
%opts    = struct( 'factr', 1e4, 'pgtol', 1e-4);
l=zeros(length(x0),1);
u=inf.*ones(length(x0),1);
fn=@(x)func(x,Z_ftir,Z_hnmr,P_ftir,P_hnmr);
gr=@(x)grad(x,Z_ftir,Z_hnmr,P_ftir,P_hnmr);
fabc= @(x)fminunc_wrapper(x,fn,gr); 

%[x,fout,info] = lbfgsb(@(x)pobobj(x,Z_ftir,Z_hnmr,P_ftir,P_hnmr),l,u,opts);
[x,fout,info] = lbfgsb(@(x)fabc(x),l,u,opts);

Aoptb=[x(1:8) x(9:16) x(17:24) x(25:32)];
Boptb=[x(33:74) x(75:116) x(117:158) x(159:200)];
H1optb=[x(201:1964) x(1965:3728) x(3729:5492) x(5493:7256)];
H2optb=[x(7257:7939) x(7940:8622) x(8623:9305) x(9306:9988)];

csvwrite('Aoptb.csv',Aoptb);
csvwrite('Boptb.csv',Boptb);
csvwrite('H1optb.csv',H1optb);
csvwrite('H2optb.csv',H2optb);

