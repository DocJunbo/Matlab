%% 例2-2
num = 4 * conv([1,2],[1,6,6]);
den = conv([1,0],conv([1,1],conv([1,1],[1 3 2 5])));
printsys(num,den)

%% 例2-3
num = [0 0 3 2;1 0 2 5];
den = [3 5 2 1];
printsys(num,den)

%% 例2-4
K = [3 4]';
Z = [-12 -1;inf -2];
P = [-3 -4 -5]';
P = [1 3 5 2];
R = roots(P)
P1 = poly(R)

%% 例2-6
clear;clc;
A = [0 0 1;-3/2 -2 -1/2;-3 0 -4];
B = [1 1;-1 -1;-1 -3];
C = [1 0 0;0 1 0];
D = zeros(2,2);
[num1,den1] = ss2tf(A,B,C,D,1)
[num2,den2] = ss2tf(A,B,C,D,2)

%% 例2-8
clear;clc;
num = [0 2 3; 1 2 1];
den = [1 0.4 1];
[A,B,C,D] = tf2ss(num,den)

%% 2-9
clear;clc;
num = [6 12 6 10];
den = [1 2 3 1 1];
[Z,P,K] = tf2zp(num,den)

%% 2-10
clear;clc;
K = 6;
Z = [-3];
P = [-1;-2;-5];
[A,B,C,D] = zp2ss(Z,P,K),[num,den] = zp2tf(Z,P,K)

%% 2-11
clear;clc;
num = [6 12 6 10];
den = [1 2 3 1 1];
[R,P,H] = residue(num,den)

%% 2-12
clear;clc;
A = [-5 8 0 0;-4 7 0 0;0 0 0 4;0 0 -2 6];
B = [4;-2;2;1];
C = [2 -2 -2 2];
D = 0;
[Am,Bm,Cm,Dm] = minreal(A,B,C,D)

%% 2-13
clear;clc;
A = [-5 8 0 0;-4 7 0 0;0 0 0 4;0 0 -2 6];
B = [4;-2;2;1];
C = [2 -2 -2 2];
D = 0;
[num,den] = ss2tf(A,B,C,D,1),[NUMm,DENn] = minreal(num,den)

%% 2-14
clear;clc;
A1 = [2 3;-1 4];
B1 = [1 0]';
C1 = [2 4];
D1 = 1;
A2 = [0 3;-3 -1];
B2 = [0 1]';
C2 = [1 3];
D2 = 2;
[A B C D] = series(A1,B1,C1,D1,A2,B2,C2,D2)

%% 2-15
clear;clc;
A1 = [2 3;-1 4];
B1 = [1 0]';
C1 = [2 4];
D1 = 1;
A2 = [0 3;-3 -1];
B2 = [0 1]';
C2 = [1 3];
D2 = 2;
[A B C D] = parallel(A1,B1,C1,D1,A2,B2,C2,D2)
num1 = 3;den1 = [1 4] ;num2 = [2 4];den2 = [1 2 3];
[num den] = parallel(num1,den1,num2,den2)

%% 2-16
clear;clc;
num1 = [2 1];den1 = [1 0];num2 = [5 2];den2 = [1 0];
g11 = tf(num1,den1);g22 = tf(num2,den2);
A = [0 3;-3 -1];B = [1 0;0 1];C = [2 1;0 1];Go = ss(A,B,C,0);
Gc = [g11,0;0,g22];
G = Go * Gc

%% 2-17
clear;clc;
numg = [2 5 1];deng = [1 2 3];numh = [5 10];denh = [1 10];
[num,den] = feedback(numg,deng,numh,denh);
printsys(num,den)

%% 2-18
clear;clc;
num1 = [10];den1 = [1 1];num2 = [1];den2 = [2 0.5];
num3 = [540];den3 = [1];num4 = [0.1];den4 = [1];
[na,da] = series(num1,den1,num2,den2);
[nb,db] = feedback(na,da,num4,den4,-1);
[nc,dc] = series(num3,den3,nb,db);
[num,den] = cloop(nc,dc,-1);
printsys(num,den)

%% 2-19
clear;clc;
s = tf('s');g11 = (2 * s + 1)/s;
g22 = (5 * s + 2)/s;
Gc = [g11,0;0,g22];
A = [0,3;-3,-1];
B = [1,0;0,1];
C = [2,1;0,1];
Go = ss(A,B,C,0);
H = eye(2);
GG = feedback(Go * Gc,H)

%% 2-20
clear;clc;
syms G1 G2 G3 G4 H1 H2 H3;
GG1 = feedbacksym(G4 * G3,H3);
GG2 = feedbacksym(GG1 * G2,H2/G4);
G = feedbacksym(GG2 * G1,H1);
pretty(G)

%% 2-24
clear;clc;
pade(0.1,3);

%% 2-25
clear;clc;
[A,B,C,D] = rmodel(3,2,2)

%% 2-26
clear;clc;
K = 6;Z = [-3];P = [-1 -2 -5]';
T = 0.1;
[A,B,C,D] = zp2ss(Z,P,K)

%% 2-27
clear;clc;
num = 1;den = conv([1,0],[1,1]);T = 1;
[numd1,dend1] = c2dm(num,den,T);
printsys(numd1,dend1,'z')

%% 2-28
clear;clc;
K = 1;Z = [0.7];P = [0.5];T = 0.1;
sys = zpk(Z,P,K,T),sys1 = d2d(sys,0.05)

%% 3-1
clear;clc;
r = 2;numo = 8;deno = [1 3 0];numh = 1;denh = 1;
[num,den] = feedback(numo,deno,numh,denh);[A,b,C,d] = tf2ss(num,den);
Tf = input('仿真时间Tf = ');h = input('计算步长h = ');
x = [zeros(length(A),1)];y = 0;t = 0;
for i = 1:Tf/h
    K1 = A * x + b * r;
    K2 = A * (x + h * K1/2) + b * r;
    K3 = A * (x + h * K2/2) + b * r;
    K4 = A * (x + h * K3) + b * r;
    x = x + h * (K1 + 2 * K2 + 2 * K3 + K4)/6;
    y = [y;C * x];t = [t;t(i)+h];
end
plot(t,y)

%% 3-2
clear;clc;
r=10;
P=[0.1 1 0.5 1;0 1 1 0;2 1 2 0;10 1 10 0];
W=[0 0 0 -1;1 0 0 0;0 1 0 0;0 0 1 0];
W0=[1;0;0;0];Wc=[0 0 0 1];
Tf=input('仿真时间Tf=');h=input('计算步长h=');
A1=diag(P(:,1));B1=diag(P(:,2));C1=diag(P(:,3));D1=diag(P(:,4));
H=B1-D1*W;Q=C1*W-A1;
A=inv(H)*Q;B=inv(H)*C1*W0;
x=[zeros(length(A),1)];y=[zeros(length(Wc(:,1)),1)];
t=0;
for i=1:Tf/h
    K1=A*x+B*r;
    K2=A*(x+h*K1/2)+B*r;
    K3=A*(x+h*K2/2)+B*r;
    K4=A*(x+h*K3)+B*r;
    x=x+h*(K1+2*K2+2*K3+K4)/6;
    y=[y,Wc*x];t=[t,t(i)+h];
end
plot(t,y)

%% 3-3
clear;clc;
r=2;
numo=8;deno=[1,3,0];[num,den]=cloop(numo,deno);[A,b,C,d]=tf2ss(num,den);                                            
Tf=input('仿真时间Tf=');h=input('计算步长h=');
x=[zeros(length(A),1)];y=0;t=0;
A=[A,b;zeros(1,length(A)),0];C=[C,0];x=[x;r];
eAt=eye(size(A))+A*h+A^2*h^2/2+A^3*h^3/(3*2)+A^4*h^4/(4*3*3);
for i=1:Tf/h
    x=eAt*x;y=[y;C*x];t=[t;t(i)+h];
end
plot(t,y)
