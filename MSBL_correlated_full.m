function [x_new,EM_variable,alpha_new]=MSBL_correlated_full(Y,A,sigma,eta)

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





%% Revoery via PC-SBL
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
% figure
while iter<iter_mx& norm(mu_new-mu_old)>1e-3
    iter=iter+1;
    mu_old=mu_new;
    mul=[mu_new(2:n,:);zeros(1,size(mu_new,2))];
    mur=[zeros(1,size(mu_new,2));mu_new(1:n-1,:)];
    var=diag(var_new);
    % varl=[var(2:n);0];
    % varr=[0;var(1:n-1)];
    % E=sum(abs(mu_new).^2,2)+eta*sum(abs(mul).^2,2)+eta*sum(abs(mur).^2,2)+var+eta*varl+eta*varr;
    alpha_new=var_new + mu_new*mu_new';
    alpha_first = alpha_new;

    % idx1=find(alpha_new>1e10);
    % alpha_new(idx1)=1e10;
    % alf=[alpha_new(2:n,:); zeros(1,size(alpha_new,2))];                                %   left-shifted version
    % arf=[ zeros(1,size(alpha_new,2)); alpha_new(1:n-1,:)];                              %   right-shifted version
    % alpha_new = diag(diag(alpha_new)) + eta*diag(diag(mu_new*mu_new',-1),-1)+eta*diag(diag(mu_new*mu_new',1),1)+...
    % eta*diag(diag(alpha_new,-2),-2)+eta*diag(diag(alpha_new,2),2);
    % 
    % alpha_new = diag(diag(alpha_new)) + eta*diag(diag(alpha_new,-1),-1)+eta*diag(diag(alpha_new,1),1);
    % alpha_new = diag(ones(n,1)) + eta*diag(ones(n-1,1),-1)+eta*diag(ones(n-1,1),1);
    % alpha_new = 2*alpha_new;
    D=inv(alpha_new);

    eigs = eig(alpha_new)
    % min(eig(alpha_new))

    objective = D*(alpha_first)*D - D;
    %=============================================
    %  estimate the variance of noise
    num=trace((Y-A*mu_old)'*(Y-A*mu_old))+trace(var_new*A'*A);
    den=m;
    sigma2=num/den;
    %==============================================

    var_new=inv(A'*A/sigma2+D);

    mu_new=1/sigma2*var_new*A'*Y;

    EM_variable(iter) = real(-0.5*log(det(alpha_new))-0.5*(trace(D*alpha_first)));
        % EM_variable(iter) = real(-0.5*log(det(alpha_new)));
    EM_variable(iter)
    
    % % sum(sum(abs(diag(objective)))) 
    % % sum(sum(abs(diag(objective,1)))) 
    % % sum(sum(abs(diag(objective,-1))))
    % 
    % subplot(1,6,1)
    % imagesc(abs(var_new))
    % title('\Sigma_s + \mu_s \mu_s^H')
    % subplot(1,6,2)
    % imagesc(abs(alpha_new))
    % title('\Gamma = tridiag(\Sigma_s + \mu_s \mu_s^H)')
    % 
    % subplot(1,6,3)
    % imagesc(abs((objective)))
    % % imagesc(abs(diag(diag(objective))))
    % % imagesc(abs(diag(diag(objective))+diag(diag(objective,1),1)+diag(diag(objective,-1),-1)))
    % title('\Gamma^{-1}(\Sigma_s + \mu_s \mu_s^H)\Gamma^{-1} -\Gamma^{-1}')
    % 
    % % subplot(1,5,4)
    % % imagesc(abs(var_new))
    % subplot(1,6,4)
    % plot(abs(diag(alpha_new)))
    % hold on
    % plot(abs(diag(alpha_new,1)))
    % legend('Diagonal elements of \Gamma','Subdiagonal elements of \Gamma' )
    % drawnow
    % hold off
    % % subplot(1,6,6)
    % % plot(abs(diag(D)))
    % % hold on
    % % plot(abs(diag(alpha_new,1)))
    % % hold on
    % % plot(abs(diag(alpha_new,-1)))
    % % sgtitle(['Diagonal: ',num2str(sum(sum(abs(diag(objective)))))])
    % 
    % 
    % 
    % % EM_variable2 =-0.5*trace(D*alpha_first)
    % % EM_variable(iter)
    % subplot(1,6,6)
    % plot(EM_variable)
    % title('-0.5log(|\Gamma|)-0.5tr(\Gamma^{-1}(\Sigma_s + \mu_s \mu_s^H))')

end

% figure
% subplot(1,4,1)
% imagesc(abs(alpha_first))
% subplot(1,4,2)
% imagesc(abs(alpha_new))
% subplot(1,4,3)
% imagesc(abs(D))
% subplot(1,4,4)
% imagesc(abs(objective))
%  drawnow


x_new=mu_new;
%mse=norm(x_new-x)^2