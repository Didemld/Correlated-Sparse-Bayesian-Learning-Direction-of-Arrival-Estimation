

%% Initialization
clear;
close all
clc
set=[15 20 25 30 35];
iter = 100;
rat = 0.1;
NMSE_set_sbl1 = zeros(length(set),iter);
NMSE_set_pcsbl = zeros(length(set),iter);
NMSE_set_sbl2 = zeros(length(set),iter);
NMSE_set_sbl3 = zeros(length(set),iter);
NMSE_set_sbl4 = zeros(length(set),iter);
NMSE_set_sbl5 = zeros(length(set),iter);

sup_set_sbl1 = zeros(length(set),iter);
sup_set_pcsbl = zeros(length(set),iter);
sup_set_sbl2 = zeros(length(set),iter);
sup_set_sbl3 = zeros(length(set),iter);
sup_set_sbl4 = zeros(length(set),iter);
sup_set_sbl5 = zeros(length(set),iter);
jj = 1;
for SNR = set

    for MC = 1:iter

        n=100;  
        m = 40;
        % signal dimension
        cp = 1;% number of measurements
        K=25;                                           % total number of nonzero coefficients
        L=3;                                            % number of nonzero blocks
        rng(MC)
        % siga2=1e-12;
        Nsnapshot = 1;
        % generate the block-sparse signal
        x=zeros(n,1);
        r=abs(randn(L,1)); r=r+1; r=round(r*K/sum(r));
        r(L)=K-sum(r(1:L-1));                           % number of non-zero coefficients in each block
        g=round(r*n/K);
        g(L)=n-sum(g(1:L-1));
        g_cum=cumsum(g);
        for i=1:L
            % generate i-th block
            seg = 1/sqrt(2)*(randn(r(i), 1) +1i*randn(r(i),1));              % generate the non-zero block
            cp = 0;
            cp2 = 0;
            R = eye(r(i)) + diag(cp*ones(r(i)-1,1),1)+diag(cp*ones(r(i)-1,1),-1)+diag(cp2*ones(r(i)-2,1),2)+diag(cp2*ones(r(i)-2,1),-2);
            seg=sqrtm(R)*seg;
            loc=randperm(g(i)-r(i));        % the starting position of non-zero block
            x_tmp=zeros(g(i), 1);
            x_tmp(loc(1):loc(1)-1+r(i))= seg;
            x(g_cum(i)-g(i)+1:g_cum(i), 1)=x_tmp;
        end

        % generate the measurement matrix
        % 
        % 
        Phi=randn(m,n);
         % Phi = randn(m,n) + sqrt(-1)*randn(m,n);
         A=Phi./(ones(m,1)*sqrt(sum(Phi.^2)));
        % noiseless measurements
        % x(x~=0) = 1;
        % cp=0.8;
        X = repmat(x,[1 Nsnapshot]);
        % X = X.*exp(1i*2*pi*randn(size(X)));
        % %

        % X(abs(X) < 0.3) = 0;
        % X(80,:) = 1;
        % x(80) =1;





        measure=A*X;

        % % Observation noise, stdnoise = std(measure)*10^(-SNR/20);
        % stdnoise=sqrt(sigma2);
        % noise=randn(m,1)*stdnoise;
        rnl = 10^(-SNR/20)*norm(x);
        nwhite = complex(randn(m,1),randn(m,1))/sqrt(2*m);
        noise = nwhite * rnl;	% error vector

        % Noisy measurements
        y=measure+noise;

        org_x = norm(mean(abs(X),2))^2;

        %% Revoery via PC-SBL
        eta=1;
        sigma2 =1;
        a=0.25;
        x_sbl1=MPCSBLa(y,A,sigma2,eta,a);
        % x_new=mu_new;

        % org_x = norm(mean(abs(X),2))^2;
        nmse_sbl1(MC)=norm(mean(abs(x_sbl1),2)-mean(abs(X),2))^2/org_x;

        % subplot(1,6,1)
        %
        %
        % stem(mean(abs(X),2),'b')
        % hold on
        % stem(mean(abs(x_new),2),'r')
        % legend('Ground Truth','Reconstructed');
        % title(['EM SBL NMSE: ', num2str(nmse)])


        %% Revoery via PC-SBL
        eta=1;
        sigma2 =1;
        a=0.5;
        x_pcsbl=MPCSBLa(y,A,sigma2,eta,a);
        % x_new=mu_new;
        % subplot(1,6,2)
        nmse_pcsbl(MC)=norm(mean(abs(x_pcsbl),2)-mean(abs(X),2))^2/org_x;

        % stem(mean(abs(X),2),'b')
        % hold on
        % stem(mean(abs(x_new),2),'r')
        % title(['EM PCSBL NMSE: ', num2str(nmse)])


        eta=1;
        sigma2 =1;
        a=0.75;
        x_sbl2=MPCSBLa(y,A,sigma2,eta,a);
        % x_new=mu_new;
        % subplot(1,5,5)
        nmse_sbl2(MC)=norm(mean(abs(x_sbl2),2)-mean(abs(X),2))^2/org_x;

        %%
        eta=1;
        sigma2 =1;
        a=1;
        x_sbl3=MPCSBLa(y,A,sigma2,eta,a);
        % x_new=mu_new;
        % subplot(1,5,5)
        nmse_sbl3(MC)=norm(mean(abs(x_sbl3),2)-mean(abs(X),2))^2/org_x;


        %%
        eta=1;
        sigma2 =1;
        a=1.25;
        x_sbl4=MPCSBLa(y,A,sigma2,eta,a);
        % x_new=mu_new;
        % subplot(1,5,5)
        nmse_sbl4(MC)=norm(mean(abs(x_sbl4),2)-mean(abs(X),2))^2/org_x;

        sup_X = zeros(size(X));
        sup_x_sbl1 = zeros(size(X));
        sup_x_sbl2 = zeros(size(X));
        sup_x_pcsbl = zeros(size(X));
        sup_x_sbl3 = zeros(size(X));
        sup_x_sbl4 = zeros(size(X));
       

        sup_X(abs(X)>rat) = 1;
        sup_x_sbl1(abs(x_sbl1) > rat)=1;
        sup_x_sbl2(abs(x_sbl2) > rat)=1;
        sup_x_pcsbl(abs(x_pcsbl) > rat)=1;
        sup_x_sbl3(abs(x_sbl3) > rat)=1;
        sup_x_sbl4(abs(x_sbl4) > rat)=1;

        sum_temp = zeros(size(X));
        sum_temp((sup_X+sup_x_sbl1)==2) =1;
        srate_sbl1(MC) = sum(sum_temp)/(sum(abs(sup_x_sbl1-sup_X)) + sum(abs(sup_X))); 

        sum_temp = zeros(size(X));
        sum_temp((sup_X+sup_x_sbl2)==2) =1;
        srate_sbl2(MC) = sum(sum_temp)/(sum(abs(sup_x_sbl2-sup_X)) + sum(abs(sup_X))); 

        sum_temp = zeros(size(X));
        sum_temp((sup_X+sup_x_pcsbl)==2) =1;
        srate_pcsbl(MC) = sum(sum_temp)/(sum(abs(sup_x_pcsbl-sup_X)) + sum(abs(sup_X))); 

        sum_temp = zeros(size(X));
        sum_temp((sup_X+sup_x_sbl3)==2) =1;
        srate_sbl3(MC) = sum(sum_temp)/(sum(abs(sup_x_sbl3-sup_X)) + sum(abs(sup_X))); 
        

        sum_temp = zeros(size(X));
        sum_temp((sup_X+sup_x_sbl4)==2) =1;
        srate_sbl4(MC)= sum(sum_temp)/(sum(abs(sup_x_sbl4-sup_X)) + sum(abs(sup_X))); 
      

    end
    NMSE_set_sbl1(jj,:) = nmse_sbl1;
    NMSE_set_pcsbl(jj,:) = nmse_pcsbl;
    NMSE_set_sbl2(jj,:) = nmse_sbl2;
    NMSE_set_sbl3(jj,:) = nmse_sbl3;
    NMSE_set_sbl4(jj,:) = nmse_sbl4;

    sup_set_sbl1(jj,:) = srate_sbl1;
    sup_set_pcsbl(jj,:) = srate_pcsbl;
    sup_set_sbl2(jj,:) = srate_sbl2;
    sup_set_sbl3(jj,:) = srate_sbl3;
    sup_set_sbl4(jj,:) = srate_sbl4;

    jj = jj + 1;
end

%%
    sup_set_sbl1(sup_set_sbl1<1) =0;
    sup_set_pcsbl(sup_set_pcsbl<1) = 0;
    sup_set_sbl2(sup_set_sbl2<1)=0;
    sup_set_sbl3(sup_set_sbl3<1)=0;
    sup_set_sbl4(sup_set_sbl4<1)=0;

data = set;
ll=2;
plot(data,(mean(NMSE_set_sbl1,2)),'LineWidth',ll)
hold on
plot(data,(mean(NMSE_set_pcsbl,2)),'LineWidth',ll)
hold on
plot(data,(mean(NMSE_set_sbl2,2)),'LineWidth',ll)
hold on
plot(data,(mean(NMSE_set_sbl3,2)),'LineWidth',ll)
hold on
plot(data,(mean(NMSE_set_sbl4,2)),'LineWidth',ll)


ylabel('NMSE')
xlabel('SNR')
legend('a=0.25', 'a=0.5','a=0.75','a=1','a=1.25')
grid on

figure

data = set;
ll=2;
plot(data,(mean(sup_set_sbl1,2)),'LineWidth',ll)
hold on
plot(data,(mean(sup_set_pcsbl,2)),'LineWidth',ll)
hold on
plot(data,(mean(sup_set_sbl2,2)),'LineWidth',ll)
hold on
plot(data,(mean(sup_set_sbl3,2)),'LineWidth',ll)
hold on
plot(data,(mean(sup_set_sbl4,2)),'LineWidth',ll)


ylabel('Support Recovery Rate')
xlabel('SNR')
legend('a=0.25', 'a=0.5','a=0.75','a=1','a=1.25')
grid on



% 
% NMSE_sr_sbl1 =zeros(size(NMSE_set_sbl1));
% NMSE_sr_sbl2 =zeros(size(NMSE_set_sbl2));
% NMSE_sr_sbl3 =zeros(size(NMSE_set_sbl3));
% NMSE_sr_sbl4 =zeros(size(NMSE_set_sbl4));
% NMSE_sr_pcsbl =zeros(size(NMSE_set_pcsbl));





% %%
%
% %%
% x_real = X;
% % For SR using pFISTA_diag_US
% cs_params.Beta       = 1;
% cs_params.L0         = [];         % Lipschitz constant
% cs_params.LambdaBar  = 1e-7;
% cs_params.Lambda     = 10^7;        % 1 - l1 regularization parameters. 0.0001 for TV, 0.009 for analysis
% cs_params.IterMax    = 150;
% cs_params.NonNegOrth = 0;
% cs_params.MaxTimeLag = 1;
% cs_params.sizeY = size(x_real);
% cs_params.Lambda= 8*10^1;        % 1 - l1 regularization parameters. 0.0001 for TV, 0.009 for analysis
% X_out_all = pFISTA_diag_US_ABZ_alternativeI_noconvA( y, A, cs_params );
% if size(X_out_all,3) > 1
%     X_out = squeeze(X_out_all);
%     X_output = reshape(X_out,size(x_real,1),size(x_real,2),size(X_out,2));
% else
%     X_out = squeeze(X_out_all);
%     X_output = reshape(X_out,size(x_real,1),size(x_real,2),size(X_out,3));
%
% end
% % x_ref = mean(abs(x_real)/max(max(abs(x_real))),3);
%
% x_mmv = double(mean(abs(X_output),3));
%
% subplot(1,6,6)
% nmse=norm(mean(abs(x_mmv),2)-mean(abs(X),2))^2/org_x
%
% stem(mean(abs(X),2),'b')
% hold on
% stem(mean(abs(x_mmv),2),'r')
%
% title(['MMV FISTA NMSE: ', num2str(nmse)])
%
% sgtitle(['m/n = ', num2str(m/n), ', K = ', num2str(K),  ', L = ', num2str(L), ', SNR = ', num2str(SNR)])


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
% figure
% imagesc(abs((X*X')))

