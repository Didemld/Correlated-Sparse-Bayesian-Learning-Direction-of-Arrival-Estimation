function [x_new,D]=MSBL_correlated_enhanced(Y,A,sigma,eta,inner_it)

%%
% Input:
% y: measurements;
% A: sensing matrix;
% sigma: square root of the noise covariance;
% beta: parameter controling the relevance between the elements;
% Output:
% x_new: estimated sparse signal

%%
% eta =0;
[m,n]=size(A);





%% Revoery via Corr-SBL
% figure
sy=size(Y,2);
iter=0;
iter_mx=100;
D=eye(n);
sigma2=1;
sigma2=sigma^2;
alpha_all=ones(n,1);
var_new=inv(A'*A/sigma2+D);
mu_old=ones(n,sy);
mu_new=1./sigma2*var_new*A'*Y;
gamma_new=1/sigma2;
alpha_all=var_new + mu_new*mu_new';
alpha_up = diag(diag(alpha_all));
while iter<iter_mx && norm(mu_new-mu_old)>1e-3
    iter=iter+1;
    mu_old=mu_new;
    % mul=[mu_new(2:n,:);zeros(1,size(mu_new,2))];
    % mur=[zeros(1,size(mu_new,2));mu_new(1:n-1,:)];
    % var=diag(var_new);
    % varl=[var(2:n);0];
    % varr=[0;var(1:n-1)];
    % E=sum(abs(mu_new).^2,2)+eta*sum(abs(mul).^2,2)+eta*sum(abs(mur).^2,2)+var+eta*varl+eta*varr;
    alpha_all=var_new + mu_new*mu_new';
    % idx1=find(alpha_new>1e10);
    % alpha_new(idx1)=1e10;
    % alf=[alpha_new(2:n,:); zeros(1,size(alpha_new,2))];                                %   left-shifted version
    % arf=[ zeros(1,size(alpha_new,2)); alpha_new(1:n-1,:)];                              %   right-shifted version
    alpha_up = diag(diag(alpha_all)) + eta*diag(diag(alpha_all,-1),-1)+eta*diag(diag(alpha_all,1),1);

    % diagup= diag(inv(alpha_up)*alpha_up*inv(alpha_up)-inv(alpha_up))

    eta2=1;
    % if iter==iter_mx && norm(mu_new-mu_old)<1e-3
    lrate = 0.001;
    for i = 1: inner_it
        alpha_inv = inv(alpha_up);
        update_part = alpha_inv*alpha_all*alpha_inv-alpha_inv;

        update_part = eta2*diag(diag(update_part,1),1)+eta2*diag(diag(update_part,-1),-1);

        update_part((update_part)<10^-6)=0;
        update_part(abs(update_part)<10^-6)=0;


        % diagup= diag(update_part)

        alpha_up = alpha_up + lrate * update_part;

        
    end
    [V,D]= eig(alpha_up);
        D = real(D);
        D((D) < 0) = 10^-6;
        alpha_up = V*D*V';
        (eig(D))
        i
    % end
    D=inv(alpha_up);
    %=============================================
    %  estimate the variance of noise
    num=trace((Y-A*mu_old)'*(Y-A*mu_old))+trace(var_new*A'*A);
    den=m;
    sigma2=num/den;
    %==============================================

    var_new=inv(A'*A/sigma2+D);

    mu_new=1/sigma2*var_new*A'*Y;

    if det(alpha_up) < 0
        'Undefined'
        break;
    end
    EM_variable(iter) = real(-0.5*log(det(alpha_up))-0.5*(trace(D*alpha_all)));
    plot(EM_variable)
    title('-0.5log(|\Gamma|)-0.5tr(\Gamma^{-1}(\Sigma_s + \mu_s \mu_s^H))')
    drawnow

end
x_new=mu_new;
% x_new(abs(x_new) < 0.3) =0;

%mse=norm(x_new-x)^2