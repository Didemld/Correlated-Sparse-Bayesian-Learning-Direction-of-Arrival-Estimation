function [x_new,D]=MSBL_alternative_source(Y,A,sigma,eta,Nsource)

%%
% Input:
% y: measurements;
% A: sensing matrix;
% sigma: square root of the noise covariance;
% beta: parameter controling the relevance between the elements;
% Output:
% x_new: estimated sparse signal

%%
eta =0;
[m,n]=size(A);





%% Revoery via PC-SBL

% Initialization of parameters
a=0.5;
b=1e-10;
c=1e-10;
d=1e-10;
sy=size(Y,2);
iter=0;
iter_mx=100;
D=eye(n);
sigma2=1;
sigma2=sigma^2;
alpha_new=ones(n,1);
var_new=inv(A'*A/sigma2+D);
mu_old=ones(n,sy);
mu_new=1./sigma2*var_new*A'*Y;
gamma_new=1/sigma2;
SCM = (Y * Y') / size(Y,2);
while iter<iter_mx& norm(mu_new-mu_old)>1e-3
    iter=iter+1;
    mu_old=mu_new;
    mul=[mu_new(2:n,:);zeros(1,size(mu_new,2))];
    mur=[zeros(1,size(mu_new,2));mu_new(1:n-1,:)];
    var=diag(var_new);
    E=sum(abs(mu_new).^2,2)+var;
    alpha_new=E;
    idx1=find(alpha_new>1e10);
    alpha_new(idx1)=1e10;
    D=diag(1./(alpha_new));
    %=============================================
    % %  estimate the variance of noise
    % num=trace((Y-A*mu_old)'*(Y-A*mu_old))+trace(var_new*A'*A);
    % den=m;
    % sigma2=num/den;
    %==============================================
    %mse=norm(x_new-x)^2 % locate same peaks for all frequencies
    [ ~ , Ilocs] = findpeaks(alpha_new,'SORTSTR','descend','NPEAKS',Nsource);
    Apeak      = A(:,Ilocs);


    % only active replicas
    Am     = Apeak;

    % noise estimate
    sigma2 = real(trace( (eye(m)-Am*pinv(Am)) * SCM ) / ( m- Nsource ) );






    var_new=inv(A'*A/sigma2+D);

    mu_new=1/sigma2*var_new*A'*Y;

end
x_new=mu_new;

%