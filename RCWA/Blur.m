function B = Blur(M,index)
% This function blur the device using index
% Input: M--the matrix which need to blur; index--blur index
% Output: blured matrix
% Caculate the kernel
ksize = index;
kernel = ones(ksize)/ksize^2;
[h,w] = size(M);

% Generate pad
p = ksize;
Ipad = zeros(h+2*p,w+2*p);
%Middle 
Ipad(p+1:p+h, p+1:p+w) = M;
%Top and Bottom
Ipad(1:p, p+1:p+w) = repmat(M(1,1:end), p, 1);
Ipad(p+h+1:end, p+1:p+w) = repmat(M(end,1:end), p, 1); 
%Left and Right
Ipad(p+1:p+h, 1:p) = repmat(M(1:end,1), 1, p);
Ipad(p+1:p+h, p+w+1:end) = repmat(M(1:end,end), 1, p); 
%Corners
Ipad(1:p, 1:p) = M(1,1); %Top-left
Ipad(1:p, p+w+1:end) = M(1,end); %Top-right
Ipad(p+h+1:end, 1:p) = M(end,1); %Bottom-left
Ipad(p+h+1:end,p+w+1:end) = M(end,end); %Bottom-right

%Generate image
[h1, w1] = size(Ipad);
kernelimagepad = zeros(h1,w1);
kernelimagepad(1:ksize, 1:ksize) = kernel;

% Perform 2D FFTs
fftimage = fft2(Ipad);
fftkernel = fft2(kernelimagepad);

%Set all zero values to minimum value
% fftkernel(fftkernel == 0) = 1e-6;

%Multiply FFTs
fftblurimagepad = fftimage.*fftkernel;

%Perform Inverse 2D FFT
blurimagepad  = ifft2(fftblurimagepad);

%Remove Padding and deplace back pattern
deplace = floor(ksize/2);
B = blurimagepad(ksize+1+deplace:ksize+h+deplace,ksize+1+deplace:ksize+w+deplace);








