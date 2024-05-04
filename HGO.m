% This code is to verfy the stress response of transverse isotropic 
% HGO model under uniform tensile.
% About HGO model, please refer to Holzapfel 2006 titiled by
% A new Constitutive Framework for Arterial Wall Mechanics and a Comparative Study of Material Models
clc;clear;close all;

I1 = @(lambda1, lambda2, lambda3) lambda1.*lambda1 + lambda2.*lambda2 + lambda3.*lambda3;

I4 = @(lambda1, lambda2, lambda3, theta) sin(theta).*sin(theta).*lambda1.*lambda1 + cos(theta).*cos(theta).*lambda3.*lambda3;

ff = @(k2, kappa, theta, lambda1, lambda2, lambda3) k2 .* ( kappa .* I1(lambda1, lambda2, lambda3)...
    +(1.0 - 3.0 .* kappa) .* I4(lambda1, lambda2, lambda3, theta) - 1.0 );

dff_dlambda1 = @(k2, kappa, theta, lambda1, lambda2, lambda3) 2.0 .* k2 .* ( kappa .* I1(lambda1, lambda2, lambda3)...
    + (1.0 - 3.0 .* kappa) .* I4(lambda1, lambda2, lambda3, theta) - 1.0 ) .* 2.0 .* (kappa .* lambda1 + (1.0 - 3.0 .* kappa)...
    .* sin(theta).*sin(theta) .* lambda1 );

dff_dlambda2 = @(k2, kappa, theta, lambda1, lambda2, lambda3) 2.0 .* k2 .* ( kappa .* I1(lambda1, lambda2, lambda3)...
    + (1.0 - 3.0 .* kappa) .* I4(lambda1, lambda2, lambda3, theta) - 1.0 ) .* 2.0 .* kappa .* lambda2;

dff_dlambda3 = @(k2, kappa, theta, lambda1, lambda2, lambda3) 2.0 .* k2 .* ( kappa .* I1(lambda1, lambda2, lambda3)...
    + (1.0 - 3.0 .* kappa) .* I4(lambda1, lambda2, lambda3, theta) - 1.0 ) .* 2.0 .* (kappa .* lambda3 + (1.0 - 3.0 .* kappa)...
    .* cos(theta).*cos(theta) .* lambda3 );

Psi = @(c, k1, k2, theta, kappa, lambda1, lambda2, lambda3) 0.5 .* c .* I1(lambda1, lambda2, lambda3)...
    + k1 ./ k2 .* ( exp( ff(k2, kappa, theta, lambda1, lambda2, lambda3).^2 ) -1.0 );

dPsi_dlambda1 = @(c, k1, k2, theta, kappa, lambda1, lambda2, lambda3) c .* lambda1 + 2.0 * k1 ./ k2...
    .* ( exp( ff(k2, kappa, theta, lambda1, lambda2, lambda3).^2 ) - 1.0 ) .* ff(k2, kappa, theta, lambda1, lambda2, lambda3)...
    .* dff_dlambda1(k2, kappa, theta, lambda1, lambda2, lambda3);

dPsi_dlambda2 = @(c, k1, k2, theta, kappa, lambda1, lambda2, lambda3) c .* lambda2 + 2.0 * k1 ./ k2...
    .* ( exp( ff(k2, kappa, theta, lambda1, lambda2, lambda3).^2 ) - 1.0 ) .* ff(k2, kappa, theta, lambda1, lambda2, lambda3)...
    .* dff_dlambda2(k2, kappa, theta, lambda1, lambda2, lambda3);

dPsi_dlambda3 = @(c, k1, k2, theta, kappa, lambda1, lambda2, lambda3) c .* lambda3 + 2.0 * k1 ./ k2...
    .* ( exp( ff(k2, kappa, theta, lambda1, lambda2, lambda3).^2 ) - 1.0 ) .* ff(k2, kappa, theta, lambda1, lambda2, lambda3)...
    .* dff_dlambda3(k2, kappa, theta, lambda1, lambda2, lambda3);

lambda1 = linspace(1.0, 1.35, 50);
P1 = zeros(length(lambda1), 1);
c = 21.4; %kPa
k1 = 1018.8; %kPa
k2 = 20.0;
theta_group = [31, 32, 33, 34, 35, 36, 37] .* pi ./ 180;
kappa_group = [0.111, 0.222, 0.333, 0.444];

theta_labels = arrayfun(@(th) sprintf('$\\theta = %.0f^\\circ$', th / (pi / 180)), theta_group, 'UniformOutput', false);
kappa_labels = arrayfun(@(k) sprintf('$\\kappa = %0.3f$', k), kappa_group, 'UniformOutput', false);
legendItems = {};

colors = ['r', 'b', 'g', 'm', 'c', 'y', 'k'];
lineStyles = {'--', ':', '-.', '-'};
figure;
hold on;
for ii = 1:length(theta_group)
    theta = theta_group(ii);
    for jj = 1:length(kappa_group)
        kappa = kappa_group(jj);
        for kk = 1:length(lambda1)
            eq1 = @(p, lambda2, lambda3) -p/lambda2 + dPsi_dlambda2(c, k1, k2, theta, kappa, lambda1(kk), lambda2, lambda3);

            eq2 = @(p, lambda2, lambda3) -p/lambda3 + dPsi_dlambda3(c, k1, k2, theta, kappa, lambda1(kk), lambda2, lambda3);

            initial_guess = [1e6; sqrt(1/lambda1(kk))];

            options = optimoptions('fsolve', 'TolFun', 1e-6, 'MaxIter', 400, 'Display', 'iter');
            sol = fsolve(@(x) [eq1(x(1), x(2), 1.0/(lambda1(kk)*x(2))), eq2(x(1), x(2), 1.0/(lambda1(kk)*x(2)))], initial_guess);

            p1 = sol(1);
            lambda2 = sol(2);
            lambda3 = 1.0/(lambda1(kk)*lambda2);
            p2 = lambda2 * dPsi_dlambda2(c, k1, k2, theta, kappa, lambda1(kk), lambda2, lambda3);
            p3 = lambda3 * dPsi_dlambda3(c, k1, k2, theta, kappa, lambda1(kk), lambda2, lambda3);
            if abs(p1 - p2) < 1e-6 && abs(p2 - p3) < 1e-6 && abs(p1 - p3) < 1e-6 && abs(lambda1(kk)*lambda2*lambda3-1) < 1e-6
                P1(kk) = -p1 / lambda1(kk) + dPsi_dlambda1(c, k1, k2, theta, kappa, lambda1(kk), lambda2, lambda3);
            else
                disp('ERROR: can not calculate correct pressure p!');
                return
            end
        end
        color = colors(mod(ii-1, length(colors)) + 1);
        lineStyle = lineStyles(mod(jj-1, length(lineStyles)) + 1);

        plot(lambda1, P1, 'Color', color, 'LineStyle', lineStyle, 'DisplayName',...
            sprintf('%s, %s', theta_labels{ii}, kappa_labels{jj}), 'LineWidth', 2.5);
        legendItem = strcat(theta_labels{ii}, ', ', kappa_labels{jj});
        legendItems{end+1} = legendItem;
        hold on
    end
end

hXLabel = xlabel('$\lambda_1$', 'interpreter', 'latex');
hYLabel = ylabel('$P_1$', 'interpreter', 'latex');

set( gca, 'Box', 'on', 'TickDir'     , 'out', ...
    'TickLength'  , [.02 .02], ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'  , ...
    'YGrid'       , 'on' , ...
    'XGrid'       , 'on' , ...
    'XColor'      , [0 0 0 ], ...
    'YColor'      , [0 0 0 ], ...
    'LineWidth'   , 2 );
set(gca,'FontSize', 25,'fontWeight','bold');
set([hXLabel, hYLabel], 'FontName', 'Helvetica', 'FontSize', 30, 'FontWeight', 'bold');

l = legend(legendItems);
set(l, 'interpreter', 'latex', 'fontsize', 25, 'box', 'off', 'location', 'SouthOutside', ...
    'Orientation', 'horizontal', 'FontWeight', 'bold', 'FontName', 'Helvetica', 'NumColumns', 3);
l.ItemTokenSize = [30, 10];

X = 60.0;
Y = X;
xMargin = 3;
yMargin = 3;
xSize = X - 2 * xMargin;
ySize = Y - 2 * yMargin;
set(gcf, 'Units','centimeters', 'Position',[5 5 xSize ySize]);
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize',[X Y]);
set(gcf, 'PaperPosition',[xMargin yMargin xSize ySize]);
set(gcf, 'PaperOrientation','portrait');

saveas(gcf, 'HGO.jpg');