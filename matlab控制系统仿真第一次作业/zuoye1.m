%% 1-2
x = [-2 -3 4];
y(1) = 3*x(1)^2 + abs(x(2)) + sqrt(x(3))
y(2) = 3*x(1)^2 - x(2) - x(3)

%% 1-3
x = [-2 -3 4];
[y(1) y(2)] = fun(x)

%% 2
x = [4,12,10,3,1].*100/(4+12+10+3+1);
subplot(2,2,1);
pie(x);
title("二维饼图");
subplot(2,2,2);
pie3(x);
title("三维饼图");
subplot(2,2,3);
bar(x);
title("条形图");

%% 3
A = [1 2;3 4];
result{1} = A.^(0.5)
result{2} = A^(0.5)
result{3} = sqrt(A)
result{4} = sqrtm(A)

%% 4
clear;clc;
syms t;
x = t * sin(t);
y = t * (1- cos(t));
dxdt = diff(x,t);
dydt = diff(y,t);
dydx = dydt/dxdt;
simplify(dydx)

%% 5
clear;clc;
syms k T z lambda;
f = k*exp(-lambda*k*T);
F = ztrans(f,k,z);
simplify(F)
pretty(F)

%% 6
clear;clc;
syms t;
eqn = (cos(t)^2)*exp(-0.1*t) - 0.5*t == 0;
t = solve(eqn,t)

%% 7
clear;clc;
syms x y;
eqn1 = x^2 + y^2 == 1;
eqn2 = x*y == 2;
[x,y] = solve(eqn1,eqn2,x,y)

%% 8
clear;clc;
syms y(x)
eqn = diff(y,x,2) - 3*diff(y,x) + 2*y - x == 0;
Dy = diff(y,x);
cond = [y(0) == 0, Dy(0) == 1];
ySol(x) = dsolve(eqn,cond)
limit(ySol,0.5)

%% 9
clear;clc;
syms f(x) g(x);
eqns = [diff(f,x) == 3*f + 4*g, diff(g,x) == -4*f + 3*g];
Df = diff(f,x);
Dg = diff(g,x);
cond = [Df(0) == 0, Dg(0) == 1];
S = dsolve(eqns,cond);
simplify(S.f)
simplify(S.g)

%% 10
clear;clc;
t = 1930:10:2020;
y = [75.995,91.972,101.1111,123.203,131.669,150.697,179.323,203.212,226.505,249.693];
t_change = 1930:1:2020;
y1 = interp1(t,y,t_change);
y2 = interp1(t,y,t_change,'nearest');
y3 = interp1(t,y,t_change,'cubic');
y4 = interp1(t,y,t_change,'spline');
subplot(2,2,1);plot(t,y,'*',t_change,y1,'.');title('liner');
subplot(2,2,2);plot(t,y,'*',t_change,y2,'.');title('nearest');
subplot(2,2,3);plot(t,y,'*',t_change,y3,'.');title('cubic');
subplot(2,2,4);plot(t,y,'*',t_change,y4,'.');title('spline');

%% 11
clear;clc;
syms x y
eqns = [sin(x - y) == 0, cos(x + y) == 0];
vars = [x y];
[solv, solu] = solve(eqns,vars)

%% 12
clear;clc;
fun = @(t)exp(-t) * abs(sin(cos(t)));
x0 = 0;
[x,fovl] = fminsearch(fun,x0)

%% 13
clear;clc;
t = 0:0.01:5;
y = exp(-t).*cos(10*t);
plot(t,y);
n = length(y);
for i = n:-1:1
    if(abs(y(i) >= 0.05))
        m = i + 1;
        break;
    end
end
ts = (m - 1)*0.01
ys = y(m)
