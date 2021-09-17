function Ai = ifft2g(F)
[W,H]=size(F);
x = (0:W-1);
y = (0:H-1);
at = 0:W-1;
bt = 0:H-1;
k = 2*pi*1i;

eiyt = (exp((k*bt')/H)).^y; 
F0iyt = F*eiyt;

eixt = exp((k*at')/W).^x; 
Ai = (eixt*F0iyt).*1/(W*H);