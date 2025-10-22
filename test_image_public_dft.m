
clear;close all;
img=imread('lena256.bmp');
% img=mean(img,3);%for color image
ss=[128,128]; % size of image
X0=imresize(img,ss);
X0=double(X0);
X=X0;
[a,b]=size(X);
% Discrete Wavelet Transform
load DWTM.mat
ww = dftmtx(128);
M=64;
R=randn(M,a);
SNR=120;
measure=R*X;
% Observation noise
stdnoise = std(measure)*10^(-SNR/20);
noise = randn(M,1) * stdnoise;

Y=measure+noise;


A=R*ww';

figure(1);
X=reshape(X,ss);
imshow(uint8(X));
title('Original Image')




%=========================================================
%                Proposed PC-SBL algorithm
%=========================================================
% x_new=MSBL_correlated(y,A,sigma2,eta);

X3=zeros(a,b);
eta=1;
for i=1:128
    rec2=PCSBL(Y(:,i),A,stdnoise(i),eta); % recover the image column by column

    rec3=MSBL_correlated(Y(:,i),A,stdnoise(i),0.5); % recover the image column by column

    % blkLen =2;
    % 
    % blkStartLoc = 1:blkLen:128;
    % learnLambda = 0;
    % 
    % Result = BSBL_EM(A,Y(:,i),blkStartLoc,learnLambda);
    % rec4 = Result.x;

    blkLen =2;

    blkStartLoc = 1:blkLen:128;
    learnLambda = 0;

    blkLen =2;
    % blkStartLoc = 1:blkLen:n;
    % = 2: small noise;
    % = 0 : no noise
    Result5 = EBSBL_BO(A, Y(:,i), blkLen, 0);
    rec5 = Result5.x;

    X2dwt(:,i)=rec2;
    X3dwt(:,i)=rec3;
    % X4dwt(:,i)=rec4;
    X5dwt(:,i)=rec5;
end

%%
figure(2);
X2=ww'*X2dwt;  %  inverse-DWT transform
X22=reshape(X2,ss);
ERR=sqrt(sum(sum((X22-X0).^2,1),2)/sum(sum(X.^2,1),2))
imshow(uint8(X22));
title('PC-SBL');

figure(3);
X3=ww'*X3dwt;  %  inverse-DWT transform
X33=reshape(X3,ss);
ERR=sqrt(sum(sum((X33-X0).^2,1),2)/sum(sum(X.^2,1),2))
imshow(uint8(X33));
title('CorrSBL');

% figure(4);
% X4=ww'*X4dwt;  %  inverse-DWT transform
% X44=reshape(X4,ss);
% ERR=sqrt(sum(sum((X44-X0).^2,1),2)/sum(sum(X.^2,1),2))
% imshow(uint8(X44));
% title('BSBL');

figure(5);
X5=ww'*X5dwt;  %  inverse-DWT transform
X55=reshape(X5,ss);
ERR=sqrt(sum(sum((X55-X0).^2,1),2)/sum(sum(X.^2,1),2))
imshow(uint8(X55));
title('EBSBL');

