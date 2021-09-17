function F = fft2g(A)
[W,H]=size(A);
x = (0:W-1);
y = (0:H-1);
at = 0:W-1;
bt = 0:H-1;
k = 2*pi*1i;

eyt = (exp((-k*y')/H)).^bt; 
F0yt = A*eyt;

ext = exp((-k*x')/W).^at; 
F = ext*F0yt;