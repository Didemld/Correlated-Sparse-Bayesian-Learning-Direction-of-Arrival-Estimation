function [x_new,EM_variable,alpha_new]=MSBL_correlated_all(Y,A,sigma,eta)

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
alpha_new=diag(ones(n,1));
var_new=inv(A'*A/sigma2+D);
mu_old=ones(n,sy);
mu_new=1./sigma2*var_new*A'*Y;
gamma_new=1/sigma2;
% figure
while iter<iter_mx& norm(mu_new-mu_old)>1e-3

    alpha_old = alpha_new;
    iter=iter+1;
    mu_old=mu_new;
    
    alpha_new=var_new + mu_new*mu_new';
    alpha_first = alpha_new;

    alpha_new = diag(diag(alpha_new)) + eta*diag(diag(alpha_new,-1),-1)+eta*diag(diag(alpha_new,1),1);

    % if eta == 0 && iter < 10
    % % alpha_new = diag(diag(alpha_first)) + 0.5*diag(diag(alpha_first,-1),-1)+0.5*diag(diag(alpha_first,1),1);
    % alpha_new = diag(diag(alpha_first));
    % end
    % alpha_new = diag(ones(n,1)) + eta*diag(ones(n-1,1),-1)+eta*diag(ones(n-1,1),1);
    % alpha_new = 2*alpha_new;

    alpha_diag = diag(diag(alpha_new));
    D=inv(alpha_new);

    % eigs = eig(alpha_first);
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

    value = A'*inv(sigma2*eye(m)+A*alpha_new*A')*(Y*Y'-sigma2*eye(m)-A*alpha_new*A')*inv(sigma2*eye(m)+A*alpha_new*A')*A;


    EM_variable(iter) = real(-0.5*log(det(alpha_new))-0.5*(trace(D*alpha_first)));
    % sigma2
    % %
    %      deger1 = trace(A*alpha_new*A')
    %
    %  deger2 = trace(A*alpha_new*A'+sigma2*eye(size(A,1)))
    % EM_variable(iter) = real(-0.5*(trace(D*(alpha_first-alpha_new))));
    % eig(D)
    % eig(alpha_first-alpha_new)
    % (eig(real(D*(alpha_first-alpha_new))))

    % EM_variable(iter) = real(0.5*(trace(inv(alpha_diag)*inv(inv(alpha_diag)+inv(alpha_new-alpha_diag))*inv(alpha_diag)*alpha_first)));
    % EM_variable(iter)
    EM_variable1 = real(-0.5*log(det(alpha_new)));
    EM_variabled1 = real(-0.5*log(det(alpha_diag)));

    EM_variable2 = real(-0.5*(trace(D*alpha_first)));
    EM_variabled2 = real(-0.5*(trace(inv(alpha_diag)*alpha_first)));


    sigma_z = sigma2*eye(m)+ A*alpha_old*A';
    trig = alpha_old*A'*inv(sigma_z)*(Y*Y')*inv(sigma_z)*A*alpha_old'- alpha_old*A'*inv(sigma_z)*A*alpha_old;
    det_val_ex=log(det(alpha_old))
    alpha_old_trig = diag(diag(alpha_old)) + eta*diag(diag(alpha_old,-1),-1)+eta*diag(diag(alpha_old,1),1);
    det_val_new = log(det(alpha_old_trig+diag(diag(trig)) + eta*diag(diag(trig,-1),-1)+eta*diag(diag(trig,1),1)))

    val_ex=log(det(alpha_old));
    val_ex2 = log(det(alpha_old_trig));
    det_val_new;
    det_val_new2 = det(alpha_old_trig) +det(diag(diag(trig)) + eta*diag(diag(trig,-1),-1)+eta*diag(diag(trig,1),1));
    trig2 = alpha_old*A'*inv(sigma_z)*A*alpha_old-alpha_old*A'*inv(sigma_z)*(Y*Y')*inv(sigma_z)*A*alpha_old';

    determin= det(trig2);
    fark = det(diag(diag(trig2)) + eta*diag(diag(trig2,-1),-1)+eta*diag(diag(trig2,1),1));

    % EM_variable(iter)
    % sum(sum(abs(diag(objective))))
    % sum(sum(abs(diag(objective,1))))
    % sum(sum(abs(diag(objective,-1))))

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
    % %
    % % % subplot(1,5,4)
    % % % imagesc(abs(var_new))
    % subplot(1,6,4)
    % plot(abs(diag(alpha_new)))
    % hold on
    % plot(abs(diag(alpha_new,1)))
    % legend('Diagonal elements of \Gamma','Subdiagonal elements of \Gamma' )
    % drawnow
    % hold off
    % subplot(1,6,5)
    % imagesc(abs((D*alpha_first)))
    % % hold on
    % % plot(abs(diag(alpha_new,1)))
    % % hold on
    % % plot(abs(diag(alpha_new,-1)))
    % % sgtitle(['Diagonal: ',num2str(sum(sum(abs(diag(objective)))))])
    %
    % % EM_variable2 =-0.5*trace(D*alpha_first)
    % % EM_variable(iter)
    % subplot(1,6,6)
    % plot(real(EM_variable))
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