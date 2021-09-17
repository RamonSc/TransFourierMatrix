function F=fft22(A)
[W,H]=size(A);
F = zeros(W,H);
f = 0;
c1=0;
k = -2*pi*1i;

for n1= 0:W-1
    for n2= 0:H-1
        for w1=0:W-1
            for w2=0:H-1
                %f = f+(A(w1+1,w2+1)*exp((-1*(2*pi*w1*n1)*1i)/W)*exp((-1*(2*pi*w2*n2)*1i)/H));
                f = f + (A(w1+1,w2+1)*exp(k*((w1*n1)/W + (w2*n2)/H))); % se demora menos
            end
        end
        F(n1+1,n2+1)=f;
        f=0;
        c1 = c1 + (100/(W*H))
    end
end