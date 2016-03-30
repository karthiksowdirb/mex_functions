beta = 0.05;
i = find(betas==beta);

%betas	 =	[0.05,      0.15,       0.25,		0.35,       0.4,        0.425,      0.45,       0.475,      0.5,        0.525,      0.55,       0.65,       0.75,       0.85,       0.95,       1.0];
binsizes =  [1,         1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          1,          2,          2,          2,          2];
rnges    =  [1,20;      1,40;       1,40;       1,40;       1,35;       1,35;       1,35;       1,49;       1,40;       1,55;       1,55;       1,35;       2,329;      1,53;       1,58;       1,38;];

binsize = binsizes(i);

[y, x] = hist(pop(i).lifetimes, 1:binsize:2000);
nzy = y>0;
y = y(nzy);
x = x(nzy);

y = y / sum(y);

xx = x(2:end);
yy = y(2:end);

myfun = @(a) a(1)*xx.^(-a(2)).*exp(-a(3)*xx);
sumofsquares = @(ab)sum((yy - myfun(ab)).^2);

ip = [0.462 0.427 0.451];
p = fminsearch(sumofsquares, ip);

n = length(xx);
dof = n-2;
sdr = sqrt( sum((yy - myfun(p)).^2) / dof);
J = jacobianest(myfun, p);
Sigma = sdr^2*inv(J'*J);
se = sqrt(diag(Sigma))';

loglog(x, y, 'b+');
hold on;
nx = 1:2000;
loglog(nx, p(1)*nx.^(-p(2)).*exp(-p(3)*nx), 'r-');
    hold off;

xlim([1 2000]);
ylim([1e-8 1]);

disp(p);
disp(se);
