clc; clear all; close all

%& Plotting for epigraphs
x = 0:0.01:2;
f1 = @(x) (x-1).^2;
figure, plot(x,f1(x), 'k'), hold on
fill(x,f1(x), [0    0.4470    0.7410]), xlabel('x'), ylabel('f(x)')
legend('f(x)', 'epi f')
set(gca,'Fontsize',20);
save_pdf_without_whitespace('f1.png')

x = -1:0.01:1;
f2 = @(x) abs(x);
figure, plot(x,f2(x), 'k'), hold on
fill(x,f2(x), [0    0.4470    0.7410]), xlabel('x'), ylabel('f(x)')
legend('f(x)', 'epi f')
set(gca,'Fontsize',20);
save_pdf_without_whitespace('f2.png')

x = -2:0.01:2;
f3 = @(x) 1/4*x.^4-x.^2/2;
figure, plot(x,f3(x), 'k'), hold on
fill(x,f3(x), [0    0.4470    0.7410]), xlabel('x'), ylabel('f(x)')
legend('f(x)', 'epi f')
set(gca,'Fontsize',20);
save_pdf_without_whitespace('f3.png')
%% Convec sets
%code to plot a heart shape in MATLAB
%set up mesh
clc; close all
t = linspace(0,2*pi,201);
r = sqrt(abs(2*sin(5*t)));
[x y]=pol2cart(t,r);
z=x;
figure, plot(x+30,y+30),fill(x+30,y+30, [0    0.4470    0.7410]), xlabel('x'), ylabel('f(x)')
axis off, axis equal
save_pdf_without_whitespace('c1.png')

t = linspace(0,2*pi,201);
r = (cos(t)+sin(t));
[x y]=pol2cart(t,r);
z=x;
figure, plot(x+30,y+30),fill(x+30,y+30, [0    0.4470    0.7410]), xlabel('x'), ylabel('f(x)')
axis off, axis equal
save_pdf_without_whitespace('c2.png')


x = 0:0.1:1;
y = x;
figure, plot([0,0,1,1,0],[0,1,1,0,0]),fill([0,0,1,1,0],[0,1,1,0,0], [0    0.4470    0.7410]), xlabel('x'), ylabel('f(x)')
axis off, axis equal
save_pdf_without_whitespace('c3.png')

%% Test proxy 1-D case
clc; close all
rng(1)
n = 5;
x = -1:0.1:1;
f1 = @(x) 1/2*x.^2;
prox_f1 = @(x,lambda) x./(lambda + 1);

lambda = 5;
x_k = 0.5;
x_k_1 = prox_f1(x_k,lambda);
x_k = (rand(n,1))*2-1;
x_k_1 = prox_f1(x_k,1);

plot(x,f1(x)), hold on,
quiver(x_k,f1(x_k),x_k_1-x_k, f1(x_k_1)-f1(x_k), '-r', 'linewidth', 2)
legend('f(x)', 'proximal updates')
save_pdf_without_whitespace('prox_f1.png')

f2 = @(x) abs(x);
x_k_1 = proxy_f2(x_k,lambda);

figure, plot(x,f2(x)), hold on,
quiver(x_k, f2(x_k), x_k_1-x_k, f2(x_k_1)-f2(x_k), '-r', 'linewidth', 2)
legend('f(x)', 'proximal updates')
save_pdf_without_whitespace('prox_f2.png')
%% Show separability of proxy
rng(1)
close all; clc
x = -1:0.1:1;
y = 0.1:0.1:1;
lambda = 5;
[X,Y] = meshgrid(x,y);

f1 = @(x,y) 1/2*(x.^2+y.^2);
f2 = @(x,y) abs(x)+abs(y);

prox_f1 = @(x,lambda) x./(lambda + 1);
Z = f1(X,Y);
figure, contour(X,Y,Z.*(Z>=0)), hold on
x_k = rand(5,1)*(max(x)-min(x))+min(x)
y_k =rand(5,1)*(max(y)-min(y))+min(y)

quiver(x_k,y_k,prox_f1(x_k,lambda)-x_k, prox_f1(y_k,lambda)-y_k, '-r', 'linewidth', 0.8)
save_pdf_without_whitespace('prox_f1_2d.png')

Z = f2(X,Y);
figure, contour(X,Y,Z.*(Z>=0)), hold on
x_k = rand(5,1)*(max(x)-min(x))+min(x)
y_k =rand(5,1)*(max(y)-min(y))+min(y)

quiver(x_k,y_k,proxy_f2(x_k,lambda)-x_k, proxy_f2(y_k,lambda)-y_k, '-r', 'linewidth', 0.8)
save_pdf_without_whitespace('prox_f2_2d.png')
%%
x = -1:0.1:1;
y = 0.1:0.1:1;

[X,Y] = meshgrid(x,y);

f1 = @(x,y) -(1/2*x.^2-y);
f1_y = @(x) (1/2*x.^2);
Z = f1(X,Y)
prox_f1 = @(x,lambda) x./(lambda + 1);

contour(X,Y,Z.*(Z>=0)), hold on
n = 1;
lambda = 0;

x = 0.5;%rand(n,1).*[-1,1];
y = 0.5;%rand(n,1).*[0.1,1];

prox_f1(x, lambda)

plot([x,prox_f1(x, lambda)],[y,f1_y(prox_f1(x, lambda))])
%quiver(x ,y ,v(1) ,v(2) ,v(3) ,s) % Plots a vector v starting from
%the point a and rescales the sise by s

function x_k_1 = proxy_f2(z,lambda)

x_k_1 = (z-lambda).*(z>=lambda) + (z+lambda).*(z<=-lambda);
end
