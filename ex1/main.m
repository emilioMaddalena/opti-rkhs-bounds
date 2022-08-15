clear all
close all
clc

rng(1)

a = 0.8; b = 8; c = 0.8;
f = @(z1,z2) 1 - a.*(z1.^2) + z2 + b*sin(c.*z2);

z_min = [-10; -10];
z_max = [10; 10];
nz = 2;

lengthscale = 5;
kernel = @(x1,x2) exp(-dist(x1,x2').^2 / (2*lengthscale^2));

Gamma = 1200;

d = 100;
Z = collect_data(d, z_min, z_max, 'grid');

fZ = f(Z(:,1),Z(:,2));

delta_bar = 1;
del = rand(d,1) * 2*delta_bar - delta_bar;

y = fZ + del;

%%

%z_gran = [1 1];
z_gran = [0.8 0.8];
[ZZ1, ZZ2] = meshgrid(z_min(1):z_gran(1):z_max(2),z_min(2):z_gran(2):z_max(2));
ZZ = [ZZ1(:) ZZ2(:)];

fZZ = f(ZZ1(:), ZZ2(:));

for i = 1:size(ZZ,1)
    [ub, lb] = opt_bnd(ZZ(i,:), [Z y], kernel, Gamma, delta_bar);
    UB(i) = ub; LB(i) = lb;
    if mod(i,10) == 0, disp(i); end
end
disp('Done optimal!')

for i = 1:size(ZZ,1)
    [ub, lb] = subopt_bnd(ZZ(i,:), [Z y], kernel, Gamma, delta_bar);
    SUB(i) = ub; SLB(i) = lb;
    if mod(i,10) == 0, disp(i); end
end
disp('Done suboptimal!')

%%
close all

sett = 'grid';

if strcmp(sett, 'rnd')
    text = ['Random sampling, $\bar \delta = ' num2str(delta_bar) '$'];
elseif strcmp(sett, 'grid')
    text = ['Grid sampling, $\bar \delta = ' num2str(delta_bar) '$'];
end

%'#77AC30'
surf(ZZ1,ZZ2,reshape(UB,size(ZZ1)), 'FaceColor', '#0072BD' ,'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on 
surf(ZZ1,ZZ2,reshape(fZZ,size(ZZ1)), 'FaceColor', 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none')

xlim([-10 10]); ylim([-10 10]); zlim([-100 50]);
xticks(-10:5:10); yticks(-10:5:10); zticks(-100:50:50);
set(gcf,'Color','w','Position',[100 100 360 250]);
title(text,'Interpreter','latex','FontSize',16); 
xl = xlabel('$z_1$','Interpreter','latex'); yl = ylabel('$z_2$','Interpreter','latex'); zl = zlabel('$f(z_1,z_2)$','Interpreter','latex');
xl.Position(1) = xl.Position(1) + abs(xl.Position(1) * 0.1); 
xl.Position(2) = xl.Position(2) + abs(xl.Position(2) * 0.15);
yl.Position(1) = yl.Position(1) + abs(yl.Position(1) * 1.5); 
yl.Position(2) = yl.Position(2) + abs(yl.Position(2) * 6.5);
yl.Position(3) = yl.Position(3) - abs(yl.Position(3) * 0.1);
zl.Position(1) = zl.Position(1) - abs(zl.Position(1) * 1.4); 
zl.Position(2) = zl.Position(2) + abs(zl.Position(2) * 0.2);
zl.Position(3) = zl.Position(3) - abs(zl.Position(3) * 4);

set(gca,'FontWeight','normal',...
        'FontSize', 16,...
        'TickLabelInterpreter','latex',...
        'FontName','Helvetica',...
        'TitleFontWeight','normal',...
        'TitleFontSize',1,...
        'view',[32,23]);

set(gca, 'Position', [0.16,0.124,0.66,0.78])
    
exportgraphics(gcf, ['ex1_' sett '_' num2str(delta_bar) '.pdf'], 'ContentType', 'vector')

%%

slice = -7.6;

ZZ2 = linspace(z_min(2),z_max(2),100);
ZZ1 = slice*ones(size(ZZ2));
ZZ = [ZZ1(:) ZZ2(:)];

fZZ = f(ZZ1(:), ZZ2(:));

for i = 1:size(ZZ,1)
    [ub, lb] = opt_bnd(ZZ(i,:), [Z y], kernel, Gamma, delta_bar);
    UB(i) = ub; LB(i) = lb;
    if mod(i,10) == 0, disp(i); end
end
disp('Done optimal!')

for i = 1:size(ZZ,1)
    [ub, lb] = subopt_bnd(ZZ(i,:), [Z y], kernel, Gamma, delta_bar);
    SUB(i) = ub; SLB(i) = lb;
    if mod(i,10) == 0, disp(i); end
end
disp('Done suboptimal!')

%%
close all

xx = ZZ2(ZZ1(:) == slice);
gt = fZZ(ZZ1(:) == slice);
ub = UB(ZZ1(:) == slice);
lb = LB(ZZ1(:) == slice);
sub = SUB(ZZ1(:) == slice);
slb = SLB(ZZ1(:) == slice);

figure
plot(xx,gt, 'k--', 'linewidth', 2); hold on;
fill([xx fliplr(xx)], [lb fliplr(ub)], 'y', 'facealpha', 0.15, 'linewidth', 1.5);

figure
plot(xx,gt, 'k--', 'linewidth', 2); hold on;
fill([xx fliplr(xx)], [slb fliplr(sub)], 'y', 'facealpha', 0.15, 'linewidth', 1.5);


figure
fill([xx fliplr(xx)], [slb fliplr(sub)], 'y', 'facealpha', 0.15, 'linewidth', 1.5); hold on
fill([xx fliplr(xx)], [lb fliplr(ub)], 'g', 'facealpha', 0.15, 'linewidth', 1.5);
plot(xx,gt, 'k--', 'linewidth', 2); hold on;

%%
legend('model', 'ground-truth', 'AutoUpdate', 'off')


xlim([-10 10]); ylim([-10 10]); zlim([-100 50]);
xticks(-10:5:10); yticks(-10:5:10); zticks(-100:50:50);
set(gcf,'Color','w','Position',[100 100 360 250]);
title(text,'Interpreter','latex','FontSize',16); 
xl = xlabel('$z_1$','Interpreter','latex'); yl = ylabel('$z_2$','Interpreter','latex'); zl = zlabel('$f(z_1,z_2)$','Interpreter','latex');
xl.Position(1) = xl.Position(1) + abs(xl.Position(1) * 0.1); 
xl.Position(2) = xl.Position(2) + abs(xl.Position(2) * 0.15);
yl.Position(1) = yl.Position(1) + abs(yl.Position(1) * 1.5); 
yl.Position(2) = yl.Position(2) + abs(yl.Position(2) * 6.5);
yl.Position(3) = yl.Position(3) - abs(yl.Position(3) * 0.1);
zl.Position(1) = zl.Position(1) - abs(zl.Position(1) * 1.4); 
zl.Position(2) = zl.Position(2) + abs(zl.Position(2) * 0.2);
zl.Position(3) = zl.Position(3) - abs(zl.Position(3) * 4);

set(gca,'FontWeight','normal',...
        'FontSize', 16,...
        'TickLabelInterpreter','latex',...
        'FontName','Helvetica',...
        'TitleFontWeight','normal',...
        'TitleFontSize',1,...
        'view',[32,23]);

set(gca, 'Position', [0.16,0.124,0.66,0.78])
    
exportgraphics(gcf, ['ex1_' sett '_' num2str(delta_bar) '.pdf'], 'ContentType', 'vector')