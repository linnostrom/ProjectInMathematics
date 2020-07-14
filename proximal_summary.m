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

%% Convex envelope
close all

x = -1:0.5:1;
f2 = @(x) abs(x);
f2_deriv = @(x,y) -1*(x<0) + 1*(x>0) + y*(x==0);
f2_p = @(x,y,yy,x_i) f2_deriv(x,yy)*x_i + (y - f2_deriv(x,y)*x);

for x_i = x
    
    %if f2(x_i) <= min(f2(x)- f2_deriv(x_i,1).*(x-x_i))
    figure(1), plot([min(x),max(x)]',[f2_p(x_i,f2(x_i),x_i,min(x)),f2_p(x_i,f2(x_i),x_i,max(x))]', 'r'), hold on
    
    if x_i == 0
        for y = -1:0.05:1
            plot([min(x),max(x)]',[f2_p(x_i,f2(x_i),y,min(x)),f2_p(x_i,f2(x_i),y,max(x))]', 'r'), hold on
        end
        A = 1
    end
    
end
figure(1), plot(x,f2(x), 'bl'), hold on
axis([min(x), max(x), 0, max(x)])
legend('f(x)', '\partialf(x)')
set(gca,'Fontsize',20);
save_pdf_without_whitespace2('subdifferential_f2.png')

x = -3:0.01:3;
f3 = @(x) 1/4*x.^4-x.^2/2;
f3_deriv = @(x) (4/4*x.^3-2*x/2);
f3_p = @(x,y,x_i) (4/4*x.^3-2*x/2)*x_i + (y - (4/4*x.^3-2*x/2)*x);
x = -3:0.1:3;
figure(2), plot(x,f3(x), 'bl'), hold on
for x_i = x
    
    if f3(x_i) <= min(f3(x)- f3_deriv(x_i).*(x-x_i))
        figure(2), plot([min(x),max(x)]',[f3_p(x_i,f3(x_i),min(x)),f3_p(x_i,f3(x_i),max(x))]', 'r'), hold on
    end
end
figure(2), plot(x,f3(x), 'bl'), hold on

axis([min(x), max(x),-1,5])
legend('f(x)', '\partialf(x)')
set(gca,'Fontsize',20);
save_pdf_without_whitespace2('subdifferential_f3.png')

x = -4:0.01:4;
f1 = @(x) (x-1).^2;
f1_deriv = @(x) 2*(x-1);
f1_p = @(x,y,x_i) 2*(x-1)*x_i + (y - 2*(x-1)*x);
x = -2:0.1:3;
figure(3), plot(x,f1(x), 'bl'), hold on
for x_i = x
    
    if f1(x_i) <= min(f1(x)- f1_deriv(x_i).*(x-x_i))
        figure(3), plot([min(x),max(x)]',[f1_p(x_i,f1(x_i),min(x)),f1_p(x_i,f1(x_i),max(x))]', 'r'), hold on
    end
end
figure(3), plot(x,f1(x), 'bl'), hold on

axis([min(x), max(x),-1,5])
legend('f(x)', '\partialf(x)')
set(gca,'Fontsize',20);
save_pdf_without_whitespace2('subdifferential_f1.png')

%% Convex sets
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

quiver(x_k,y_k,prox_f1(x_k,lambda)-x_k, prox_f1(y_k,lambda)-y_k, '-r', 'linewidth', 0.8), axis off
save_pdf_without_whitespace('prox_f1_2d.png')

Z = f2(X,Y);
figure, contour(X,Y,Z.*(Z>=0)), hold on
x_k = rand(5,1)*(max(x)-min(x))+min(x)
y_k =rand(5,1)*(max(y)-min(y))+min(y)

quiver(x_k,y_k,proxy_f2(x_k,lambda)-x_k, proxy_f2(y_k,lambda)-y_k, '-r', 'linewidth', 0.8), axis off
save_pdf_without_whitespace('prox_f2_2d.png')

%% Show separability of proxy in 3D
rng(1)
close all; clc
x = -1:0.1:1;
y = -1:0.1:1;
z = x.^2+(y.^2);
lambda = 5;
[X,Y] = meshgrid(x,y);
Z = X.^2+Y.^2;

f1 = @(x,y,z) 1/2*(x.^2+y.^2);
f2 = @(x,y) abs(x)+abs(y);

prox_f1 = @(x,lambda) x./(lambda + 1);
C = f1(X,Y);
figure, contour3(X,Y,Z), hold on

x_k = rand(5,1)*(max(x)-min(x))+min(x)
y_k =rand(5,1)*(max(y)-min(y))+min(y)
z_k =rand(5,1)*(max(z)-min(z))+min(z)

quiver3(x_k,y_k,z_k,prox_f1(x_k,lambda)-x_k, prox_f1(y_k,lambda)-y_k, prox_f1(z_k,lambda)-z_k, '-r', 'linewidth', 0.8), axis off
save_pdf_without_whitespace('prox_f1_3d.png')

Z = f2(X,Y);
figure, contour3(X,Y,Z), hold on
x_k = rand(5,1)*(max(x)-min(x))+min(x)
y_k =rand(5,1)*(max(y)-min(y))+min(y)
z_k =rand(5,1)*(max(z)-min(z))+min(z)

quiver3(x_k,y_k,z_k,proxy_f2(x_k,lambda)-x_k, proxy_f2(y_k,lambda)-y_k, proxy_f2(z_k,lambda)-z_k, '-r', 'linewidth', 0.8), axis off
save_pdf_without_whitespace('prox_f2_3d.png')

figure, contour3(X,Y,Z+C), hold on
x_k = rand(5,1)*(max(x)-min(x))+min(x)
y_k =rand(5,1)*(max(y)-min(y))+min(y)
z_k =rand(5,1)*(max(z)-min(z))+min(z)
quiver3(x_k,y_k,z_k,proxy_f2(x_k,lambda)+prox_f1(x_k,lambda)-x_k, proxy_f2(y_k,lambda)+prox_f1(y_k,lambda)-y_k, proxy_f2(z_k,lambda)+prox_f1(z_k,lambda)-z_k, '-r', 'linewidth', 0.8), axis off
save_pdf_without_whitespace('prox_f12_3d.png')

%% Show separability of proxy in 3D- different lambda
rng(1)
close all; clc
x = -1:0.1:1;
y = -1:0.1:1;
z = x.^2+(y.^2);

x_k = rand(2,1)*(max(x)-min(x))+min(x);
y_k =rand(2,1)*(max(y)-min(y))+min(y);
z_k =rand(2,1)*(max(z)-min(z))+min(z);

c = {'-.dr','--dg','-.dk','--dbl'};
i = 1;

[X,Y] = meshgrid(x,y);
Z = X.^2+Y.^2;

f1 = @(x,y,z) 1/2*(x.^2+y.^2);
f2 = @(x,y) abs(x)+abs(y);

prox_f1 = @(x,lambda) x./(lambda + 1);
C = f1(X,Y);
figure(1), contour(X,Y,Z.*(Z>=0))
Z = f2(X,Y);
figure(2), contour(X,Y,Z.*(Z>=0))
figure(3), contour(X,Y,(Z+C).*((Z+C)>=0))
lambda_vec = [10,5,3,1];

for lambda = lambda_vec
    
    figure(1), hold on, plot([x_k,prox_f1(x_k,lambda)]',[y_k,prox_f1(y_k,lambda)]',c{i})
    axis off
    figure(1), save_pdf_without_whitespace('prox_f1_3d_lambda.png')
    
    Z = f2(X,Y);
    figure(2), hold on, plot([x_k,proxy_f2(x_k,lambda)]',[y_k,proxy_f2(y_k,lambda)]',c{i})
    axis off
    figure(2), save_pdf_without_whitespace('prox_f2_3d_lambda.png')
    
    figure(3), hold on,  plot([x_k,prox_f1(x_k,lambda)+proxy_f2(x_k,lambda)]',[y_k,prox_f1(y_k,lambda)+proxy_f2(y_k,lambda)]',c{i})
    axis off
    figure(3),save_pdf_without_whitespace('prox_f12_3d_lambda.png')
    i = i+1;
    
end

function x_k_1 = proxy_f2(z,lambda)
x_k_1 = (z-lambda).*(z>=lambda) + (z+lambda).*(z<=-lambda);
end
