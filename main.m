clc,clear;

% piece-wise constant parameters
Theta = @(t) [3.63; -8.58; 20.2; -21.9; -51.88] ...
             + [-22.22; 23.74; -82.66; 31.45; 73.33] * (t>=4 ) ...
             - [-22.22; 23.74; -82.66; 31.45; 73.33] * (t>=13);
         
% reference signal
t_r = [8, 17];
ref = @(t) 1.0*(t>=0 & t<=t_r(1)) + 0.0*(t>t_r(1) & t<=t_r(2)) ...
            + 1.0*(t>t_r(2));

AWR_CMRAC_Results = AWR_CMRAC;

Glush_CMRAC_Results = Glush_CMRAC;

t = AWR_CMRAC_Results.t;      
Thetas = Theta(t);



figure('name', 'state and control', 'position', [200,150,1100,1100])
subplot(311)
p1=plot(t, AWR_CMRAC_Results.x(1,:), 'color', [1.0, 0.49, 0.25], 'linewidth', 2);grid on;hold on
plot(t, Glush_CMRAC_Results.x(1,:), 'color', [0.2, 0.63, 0.79], 'linestyle','--','linewidth', 2)
plot(t, ref(t), 'k--', 'linewidth', 2)
lgd = legend({'\, AWR-CMRAC', '\, Ref [29]', '\, reference'}, ...
               'interpreter','latex','fontsize', 22, 'box', 'off');  
lgd.ItemTokenSize = [60,20];  % set the legend length
set(get(p1,'parent'),'linewidth',1.9)  % set linewidth of axes
set(gca,'Position',[0.1,0.69,0.87,0.27],'fontsize',22)
set(gca, 'GridColor', 'k', 'GridLineStyle', ':')
set(gca,'XTicklabel',[]);
ylabel('$x_1$','interpreter','latex','fontsize',22)
ylim([-0.1, 1.1])
subplot(312)
p1=plot(t, AWR_CMRAC_Results.x(2,:), 'color', [1.0, 0.49, 0.25], 'linewidth', 2);grid on;hold on
plot(t, Glush_CMRAC_Results.x(2,:), 'color', [0.2, 0.63, 0.79], 'linestyle','--','linewidth', 2)
set(get(p1,'parent'),'linewidth',1.9)
set(gca,'Position',[0.1,0.38,0.87,0.27],'fontsize',22)
set(gca, 'GridColor', 'k', 'GridLineStyle', ':')
set(gca,'XTicklabel',[]);
ylabel('$x_2$','interpreter','latex','fontsize',22)
ylim([-1.1,1.5])
subplot(313)
p1=plot(t, AWR_CMRAC_Results.u, 'color', [1.0, 0.49, 0.25], 'linewidth', 2);grid on;hold on
plot(t, Glush_CMRAC_Results.u, 'color', [0.2, 0.63, 0.79], 'linestyle','--', 'linewidth', 2)
set(get(p1,'parent'),'linewidth',1.9)
set(gca,'Position',[0.1,0.07,0.87,0.27],'fontsize',22)
set(gca, 'GridColor', 'k', 'GridLineStyle', ':')
ylabel('$u$','interpreter','latex','fontsize',22)
xlabel('$t$','interpreter','latex','fontsize',22)
ylim([-40, 80])

figure('name', 'parameters', 'position', [1350,150,1100,1100])
subplot(311)
p1=plot(t, Thetas(1,:), 'color', [1.0, 0.27, 0.0], 'linewidth', 2, 'linestyle','--');grid on;hold on
set(get(p1,'parent'),'linewidth',1.9)
plot(t, AWR_CMRAC_Results.theta_est(1,:), 'color', [1.0, 0.27, 0.0], 'linewidth', 2);hold on
plot(t, Thetas(2,:), 'color', [0.6, 0.19, 0.8], 'linewidth', 2, 'linestyle','--');hold on
plot(t, AWR_CMRAC_Results.theta_est(2,:), 'color', [0.6, 0.19, 0.8], 'linewidth', 2);hold on
plot(t, Thetas(3,:), 'color', [0.85, 0.65, 0.13], 'linewidth', 2, 'linestyle','--');hold on
plot(t, AWR_CMRAC_Results.theta_est(3,:), 'color', [0.85, 0.65, 0.13], 'linewidth', 2);hold on
plot(t, Thetas(4,:), 'color', [0.2, 0.63, 0.79], 'linewidth', 2, 'linestyle','--');hold on
plot(t, AWR_CMRAC_Results.theta_est(4,:), 'color', [0.2, 0.63, 0.79], 'linewidth', 2);hold on
plot(t, Thetas(5,:), 'color', [0.24, 0.71, 0.44], 'linewidth', 2, 'linestyle','--');
plot(t, AWR_CMRAC_Results.theta_est(5,:), 'color', [0.24, 0.71, 0.44], 'linewidth', 2);hold on
set(gca,'Position',[0.1,0.69,0.87,0.27],'fontsize',22)
set(gca, 'GridColor', 'k', 'GridLineStyle', ':')
set(gca,'XTicklabel',[]);
ylabel('$\hat\Theta$','interpreter','latex','fontsize',22)
ylim([-75, 100])
subplot(312)
p1=plot(t, Thetas(1,:), 'color', [1.0, 0.27, 0.0], 'linewidth', 2, 'linestyle','--');grid on;hold on
set(get(p1,'parent'),'linewidth',1.9)
plot(t, Glush_CMRAC_Results.theta_est(1,:), 'color', [1.0, 0.27, 0.0], 'linewidth', 2);hold on
plot(t, Thetas(2,:), 'color', [0.6, 0.19, 0.8], 'linewidth', 2, 'linestyle','--');hold on
plot(t, Glush_CMRAC_Results.theta_est(2,:), 'color', [0.6, 0.19, 0.8], 'linewidth', 2);hold on
plot(t, Thetas(3,:), 'color', [0.85, 0.65, 0.13], 'linewidth', 2, 'linestyle','--');hold on
plot(t, Glush_CMRAC_Results.theta_est(3,:), 'color', [0.85, 0.65, 0.13], 'linewidth', 2);hold on
plot(t, Thetas(4,:), 'color', [0.2, 0.63, 0.79], 'linewidth', 2, 'linestyle','--');hold on
plot(t, Glush_CMRAC_Results.theta_est(4,:), 'color', [0.2, 0.63, 0.79], 'linewidth', 2);hold on
plot(t, Thetas(5,:), 'color', [0.24, 0.71, 0.44], 'linewidth', 2, 'linestyle','--');
plot(t, Glush_CMRAC_Results.theta_est(5,:), 'color', [0.24, 0.71, 0.44], 'linewidth', 2);hold on
set(gca,'Position',[0.1,0.38,0.87,0.27],'fontsize',22)
set(gca, 'GridColor', 'k', 'GridLineStyle', ':')
set(gca,'XTicklabel',[]);
ylabel('$\hat\Theta$','interpreter','latex','fontsize',22)
ylim([-75, 100])
subplot(313)
p1=plot(t, sum(AWR_CMRAC_Results.errs.^2).^(1/2), 'color', [1.0, 0.49, 0.25], 'linewidth', 2);grid on;hold on
plot(t, sum(Glush_CMRAC_Results.errs.^2).^(1/2), 'color', [0.2, 0.63, 0.79], 'linestyle','--', 'linewidth', 2)
lgd = legend({'\, AWR-CMRAC', '\, Ref [29]'}, ...
               'interpreter','latex','fontsize', 22, 'box', 'off');  
lgd.ItemTokenSize = [60,20];  % set the legend length
set(get(p1,'parent'),'linewidth',1.9)
set(gca,'Position',[0.1,0.07,0.87,0.27],'fontsize',22)
set(gca, 'GridColor', 'k', 'GridLineStyle', ':')
ylabel('$||\xi||$','interpreter','latex','fontsize',22)
xlabel('$t$','interpreter','latex','fontsize',22)
ylim([-20, 200])
