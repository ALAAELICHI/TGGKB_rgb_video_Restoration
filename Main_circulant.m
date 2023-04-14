clc
clear all
Normalize = @(x) (x - min(x(:)))/(max(x(:)) - min(x(:)));
vid = VideoReader('xylophone.mp4'); %read the video
numImgs = get(vid, 'NumberOfFrames');% get the number of frames
frames = read(vid);
N=240;
nbrframes=10; %number of frames studied
nrhs=3*nbrframes;
X =frames(1:N,1:N,:,1:nbrframes);
X=reshape(X,[240,240,nrhs]);
%A=blur(N,4,2);
[PSF, center] = psfGauss(N, 2);
[A, AA] = kronDecomp(PSF, center, 'periodic');
AA=zeros(N,N,nrhs);
BB=zeros(N,N,nrhs);
AA(:,:,1)=A;
BB(:,:,1)=A;
AA=fft(AA,[],3);
BB=fft(BB,[],3);
AA=ndsparse(AA);
BB=ndsparse(BB);
%-------------------------------------------------------------------------
B_exact= operator(AA,BB,X); 
E = randn(size(B_exact));
E = E/tnorm(E) ; 
nu=0.001*tnorm(B_exact)*E;
B = B_exact+nu;
%-------------------------------------------------------------------------
s=tnorm(nu);
eta=1.1;
mu0=10;
%-------------------------------------------------------------------------
tic
[Xmu,mu,m]=Tensor_GGKB(AA,BB,B,mu0,eta,s);
toc
Xmu=reshape(Xmu,[N,N,3,nbrframes]);
B=reshape(B,[N,N,3,nbrframes]);
X=reshape(X,[N,N,3,nbrframes]);
% RE=norm(X-Xmu)/norm(X);
% fprintf("RE=\n");disp(RE)
% SNR=snr(double(X), double(Xmu));
% fprintf("SNR=\n");disp(SNR)
% fprintf("Tensor GGKB steps=\n");disp(m)
% fprintf("Regularization parameter=\n");disp(1/mu)
% subplot(131); imshow(Normalize(double(X(:,:,:,5))),[]); title('Original Image');
% subplot(132);imshow(Normalize(double(B(:,:,:,5))),[]);title('Blurred and Noisy Image');
% subplot(133);imshow(Normalize(double(Xmu(:,:,:,5))),[]);title('Restored Image');
RE=tnorm(double(X)-double(Xmu))/tnorm(double(X))
 fprintf("RE=\n");disp(RE)
 SNR=snr(double(X), double(Xmu));
 fprintf("SNR=\n");disp(SNR)
 fprintf("Tensor GGKB steps=\n");disp(m)
 fprintf("Regularization parameter=\n");disp(1/mu)
% % subplot(131); imshow(double(X),[]); title('Original Image');
% % subplot(132);imshow(double(B),[]);title('Blurred and Noisy Image');
% % subplot(133);imshow(double(Xmu),[]);title('Restored Image');
% h1=figure;
%  imshow(Normalize(double(X(:,:,:,5))),[]);
% set(h1,'PaperSize',[6.5 6]); %set the paper size to what you want  
% print(h1,'papavo','-dpdf') % then print it
% h2=figure;
% imshow(Normalize(double(B(:,:,:,5))),[]);
% set(h2,'PaperSize',[6.5 6]); %set the paper size to what you want  
% print(h2,'papavb','-dpdf') % then print it
% h3=figure;
% imshow(Normalize(double(Xmu(:,:,:,5))),[]);
% set(h3,'PaperSize',[6.5 6]); %set the paper size to what you want  
% print(h3,'papavr','-dpdf') % then print it