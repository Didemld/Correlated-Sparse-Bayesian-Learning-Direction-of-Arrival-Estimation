
%% Initialization
clear;
close all
clc

n=100;                                          % signal dimension
m= 40;
cp = 1;% number of measurements
K=25;                                           % total number of nonzero coefficients
L=4;                                            % number of nonzero blocks
% rng(5000)
SNR = 20;    % Signal-to-noise ratio
% sigma2=1e-12;
Nsnapshot = 10;
rng(5)
% generate the block-sparse signal
x=zeros(n,1);
r=abs(randn(L,1)); r=r+1; r=round(r*K/sum(r));
r(L)=K-sum(r(1:L-1));                           % number of non-zero coefficients in each block
g=round(r*n/K);
g(L)=n-sum(g(1:L-1));
g_cum=cumsum(g);
for i=1:L
    % generate i-th block
    seg=rand(r(i),1);                % generate the non-zero block
    loc=randperm(g(i)-r(i));        % the starting position of non-zero block
    x_tmp=zeros(g(i), 1);
    x_tmp(loc(1):loc(1)-1+r(i))= seg;
    x(g_cum(i)-g(i)+1:g_cum(i), 1)=x_tmp;
end

% generate the measurement matrix
Phi=randn(m,n);
%A=Phi./(ones(m,1)*sqrt(sum(Phi.^2)));
A=Phi;
% noiseless measurements
x(x~=0) = 1;
% cp=0.8;
X = repmat(x,[1 Nsnapshot]);
% X(70,:) = 1;
% x(70) =1;
X = X.*exp(1i*2*pi*randn(size(X)));

% cp = 0.8;
% cp2 = 0;
% R = eye(n) + diag(cp*ones(n-1,1),1)+diag(cp*ones(n-1,1),-1)+diag(cp2*ones(n-2,1),2)+diag(cp2*ones(n-2,1),-2);
% X=sqrtm(R)*X;
% 

measure=A*X;

% % Observation noise, stdnoise = std(measure)*10^(-SNR/20);
% stdnoise=sqrt(sigma2);
% noise=randn(m,1)*stdnoise;
rnl = 10^(-SNR/20)*norm(x);
nwhite = complex(randn(m,1),randn(m,1))/sqrt(2*m);
noise = nwhite * rnl;	% error vector

% Noisy measurements
y=measure+noise;


%% Revoery via PC-SBL
eta=0;
sigma2 =1;
x_new=MPCSBL(y,A,sigma2,eta);
% x_new=mu_new;

org_x = norm(mean(abs(X),2))^2;
subplot(1,6,1)

nmse=norm(mean(abs(x_new),2)-mean(abs(X),2))^2/org_x

stem(mean(abs(X),2),'b')
hold on
stem(mean(abs(x_new),2),'r')
legend('Ground Truth','Reconstructed');
title(['EM SBL NMSE: ', num2str(nmse)])


%% Revoery via PC-SBL
eta=1;
sigma2 =1;
x_new=MPCSBL(y,A,sigma2,eta);
% x_new=mu_new;
subplot(1,6,2)
nmse=norm(mean(abs(x_new),2)-mean(abs(X),2))^2/org_x

stem(mean(abs(X),2),'b')
hold on
stem(mean(abs(x_new),2),'r')
title(['EM PCSBL NMSE: ', num2str(nmse)])



%% Revoery via FP-SBL
options = SBLSet();
options.Nsource = K;
options.convergence.error   = 10^(-2);

[x_new,report]=SBL( A , y, options );
x_new = reshape(x_new,size(X));

subplot(1,6,3)

nmse=norm(mean(abs(x_new),2)/max(mean(abs(x_new),2))-mean(abs(X),2))^2/org_x

stem(mean(abs(X),2),'b')
hold on
stem(mean(abs(x_new),2)/max(mean(abs(x_new),2)),'r')
% figure
% plot(report.results.error)
title(['FP SBL NMSE: ', num2str(nmse)])


 %% Revoery via FP-SBL
% options = SBLSet();
options.Nsource = K;
options.beta=0.1;
 [x_new,report]=SBL_PC( A , y, options );
x_new = reshape(x_new,size(X));

subplot(1,6,4)

nmse=norm(mean(abs(x_new),2)/max(mean(abs(x_new),2))-mean(abs(X),2))^2/org_x

stem(mean(abs(X),2),'b')
hold on
stem(mean(abs(x_new),2)/max(mean(abs(x_new),2)),'r')
title(['FP PCSBL NMSE: ', num2str(nmse)])


 %% Revoery via FP-SBL
% options = SBLSet();
options.Nsource = K;
options.beta=1;
options.fixedpoint = 2;
[x_new,report]=SBL_combined_tridiagonal_mv( A , y, options );
x_new = reshape(x_new,size(X));

% x_new(abs(x_new)>1) = 1;
subplot(1,6,5)
nmse=norm(mean(abs(x_new),2)/max(mean(abs(x_new),2))-mean(abs(X),2))^2/org_x

stem(mean(abs(X),2),'b')
hold on
stem(mean(abs(x_new),2)/max(mean(abs(x_new),2)),'r')
title(['FP TriSBL NMSE: ', num2str(nmse)])


%%

%%
x_real = X;
% For SR using pFISTA_diag_US
cs_params.Beta       = 1;
cs_params.L0         = [];         % Lipschitz constant
cs_params.LambdaBar  = 1e-7;
cs_params.Lambda     = 10^7;        % 1 - l1 regularization parameters. 0.0001 for TV, 0.009 for analysis
cs_params.IterMax    = 150;
cs_params.NonNegOrth = 0;
cs_params.MaxTimeLag = 1;
cs_params.sizeY = size(x_real);
cs_params.Lambda= 8*10^1;        % 1 - l1 regularization parameters. 0.0001 for TV, 0.009 for analysis
X_out_all = pFISTA_diag_US_ABZ_alternativeI_noconvA( y, A, cs_params );
if size(X_out_all,3) > 1
    X_out = squeeze(X_out_all);
    X_output = reshape(X_out,size(x_real,1),size(x_real,2),size(X_out,2));
else
    X_out = squeeze(X_out_all);
    X_output = reshape(X_out,size(x_real,1),size(x_real,2),size(X_out,3));
    
end
% x_ref = mean(abs(x_real)/max(max(abs(x_real))),3);

x_mmv = double(mean(abs(X_output),3));

subplot(1,6,6)
nmse=norm(mean(abs(x_mmv),2)-mean(abs(X),2))^2/org_x

stem(mean(abs(X),2),'b')
hold on
stem(mean(abs(x_mmv),2),'r')

title(['MMV FISTA NMSE: ', num2str(nmse)])

sgtitle(['m/n = ', num2str(m/n), ', K = ', num2str(K),  ', L = ', num2str(L), ', SNR = ', num2str(SNR)])


% %%
% x_new = A'*y;
% x_new = reshape(x_new,size(X));
% 
% % x_new(abs(x_new)>1) = 1;
% figure
% nmse=norm(mean(abs(x_new),2)-mean(abs(X),2))^2/org_x
% 
% stem(abs(x),'b')
% hold on
% stem(mean(abs(x_new),2),'r')
% title(['FP TriSBL NMSE: ', num2str(nmse)])
figure
imagesc(abs((X*X')))

