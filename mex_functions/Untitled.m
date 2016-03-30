clc;

gymax = 20;
symax = gymax;

gyrs = sqrt(gradient(y).^2); 
ny = gyrs<=gymax; 
my = mean(y(ny)); 
sy = sqrt((y-my).^2); 
ny = sy<=gymax;

plot(x, y, x(ny), y(ny))

%%
rnge=2:10; 

xx = x(rnge);
yy = y(rnge);

p = polyfit(log(xx), log(yy), 1);
p(2) = exp(p(2));

px = 1:0.01:20;
py = p(2)*px.^p(1);

hold on;
plot(px, py);
hold off;