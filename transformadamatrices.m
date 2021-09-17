%% La transformada de fourier pero con matrices y dos for
clear; close all; clc;

AB1 = imread('D:\Documents\MATLAB\open.png');

%datos iniciales

% El plan de ataque consiste en volver la matriz
% de la imagen un unico valor y multiplicarlo con
% la matriz hecha con todos los valores de las
% exponenciales, eso columna por columna sería la
% suma de cada valor, al final da un vector y toca
% recontruir la matriz, teniendo en cuenta las
% transpuestas y en este caso hay una conjugada de
% por medio.

A=single(AB1(:,:,3));
%A = randi(10,3,3);
[W,H] = size(A);
k = -2*pi*1i;

A1 = reshape(A', [1,W*H]);

Exp = zeros(W*H, W*H);
c = 0;
b = 0;
a = 0;
r = 0;

tic
for x = 1:W*H
    for y = 1:W*H
        %a
        %b
        %r
        %c
        %Exp(x,y) = exp((k/W)*(a)*(r))*exp((k/H)*(b)*(c));
        Exp(x,y) = exp(k*((a*r)/W + (b*c)/H));
        
        b = b + 1;
        
        if b == H
            b = 0;
            a = a + 1;
        end
    end
    
    c = c + 1;
    a = 0;
    b = 0;
    
    if c == H
        c = 0;
        r = r + 1; 
    end
end

F1 = A1*Exp;

F = reshape(F1,H,W)';
F = conj(F);
toc

tic
fft22(A); % la funcion que hice para probar cual es mas rapida
toc

Fprueba = fft2(A); %la funcion de matlab para corroborar que este bien

%% mask que se llena de unos
clear; close all; clc;

AB1 = imread('D:\Documents\MATLAB\imagenn2.jpg');
A=single(AB1(:,:,3));
figure;imshow(uint8(A))
[W,H] = size(A);

if W>=H
    n=H;
else
    n=W;
end

Fp = fft2(A);
mask = zeros(W,H);
Ai = ifft2(Fp);

figure;imshow(uint8(Ai))
figure;
tic
for i = 1:n-1
    mask([1:i n-i:n],[1:i n-i:n]) = 1;
    
    Fmask = Fp.*mask;
    
    Amask = ifft2(Fmask);
    Amask = real(Amask);
    
    imshow(uint8((Amask)))
end
toc

%% transformada sin for
clear; close all; clc;

AB1 = imread('D:\Documents\MATLAB\open.png');
% 
A=single(AB1(:,:,3));
%A = randi(10,3,3);
[W,H] = size(A);
k = 2*pi*1i;
a = 0;
b = 0;
x = (0:W-1);
y = (0:H-1);

tic
Ftry1 = fft2(A);
toc

ey = exp((k*b*y')/H); %el vector de exponenciales cuando y varia
F0y = A(:,:)*ey;%un vector que queda de tal forma que 
% en cada valor, se multiplica el valor de A
% avanzando en x

ex = exp((k*a*x)/W); %el vector de exponenciales cuando avanza x varia
F00 = ex*F0y;%El F(0,0) ya que solo queda multiplicar 
%las constantes de exp con x variando a la
%sumatoria de y

%%
clear; close all; clc;

AB1 = imread('D:\Documents\MATLAB\open.png');
A=single(AB1(:,:,3));
Ftry = fft2(A);

%Ahora generalizando cuando(a,b) varía, en este
%caso se usa propiedades de las exponenciales para
%armar las matrices
tic
F_ab = fft2g(A);
toc

% Ahora la inversa general
Aitry = ifft2(Ftry);

tic
Ai_ab = ifft2g(F_ab);
toc

Ai_abr=real(Ai_ab);
a=Ai_abr + min(min(Ai_abr));
figure;imshow(uint8(a))

%% Recuperar la imagen limitando las frecuencias de F
clear; close all; clc;

AB1 = imread('D:\Documents\MATLAB\open.png');
A=single(AB1(:,:,3));
[W,H] = size(A);
ai = 1; %i limite inferior
as = 2; %s limie superior
bi = 1;
bs = 5; 

% limites con un recuadro exterior de 0
tic
FL1 = fft2g(A);
FL1([1:ai, W-as:end], :) = 0;
FL1(:, [1:bi, H-bs:end]) = 0;

AiL1 = ifft2g(FL1);
toc
ALr1=real(AiL1);
a1=ALr1 - min(min(ALr1));
figure;imshow(uint8(a1))

%Limites con un recuadro interior de 0

FL2 = fft2g(A);
FL2(ai:W-as, bi:H-bs) = 0;

AiL2 = ifft2g(FL2);

ALr2=real(AiL2);
a2=ALr2 - min(min(ALr2));
figure;imshow(uint8(a2))

















