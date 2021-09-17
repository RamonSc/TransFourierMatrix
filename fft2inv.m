function Ai=fft2inv(F)
[W,H]=size(F);
Ai = zeros(W,H);
c2=0;
f = 0;
for x = 0:W-1
    for y = 0:H-1
        for a = 0:W-1
            for b = 0:H-1
                f = f + F(a+1,b+1)*exp(((2*pi*x*a)*1i)/W)*exp(((2*pi*y*b)*1i)/H);
            end
        end
        Ai(x+1,y+1) = f*(1/(W*H));
        f = 0;
        c2 = c2 + (100/(W*H))
    end
end