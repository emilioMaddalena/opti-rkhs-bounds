%
% Setting for the phase portrait figures (temporary.mat):
% N = 8
% KRR data = 400
% Bounds data = 800
%
% delta_bar = 0.01
% delta_aug = 1.2; 
% gamma_aug = 1.2;
%

figure(1)

title('Sub-optimal bounds','Interpreter','latex','FontSize',16); sgtitle('');
legend('Position',[0.5,0.79,0.12,0.05],'Interpreter','latex','NumColumns',2)
xlabel('$c_A$','Interpreter','latex'); 
ylabel('$c_B$','Interpreter','latex'); 
set(gca,'FontWeight','normal',...
        'FontSize', 16,...
        'TickLabelInterpreter','latex',...
        'FontName','Helvetica',...
        'TitleFontWeight','normal',...
        'TitleFontSize',1);
set(gcf,'Color','w');
set(gcf,'Position',[100 100 350 250]);
h = findall(gcf,'Type','Line');
for i = 1:numel(h), uistack(h(i),'top'); end

exportgraphics(gcf, 'cstr_phase_subopt.pdf', 'ContentType', 'vector')

% set(gca, 'FontName', 'Latin Modern Roman')
% xlabel('$x$', 'FontName', 'Helvetica', 'Interpreter', 'latex', 'FontSize', 15); 
% ylabel('$f(x)$', 'Interpreter', 'latex', 'FontSize', 15); ax = gca; ax.FontSize = 14;
%xticks([1 1.5 2 2.5 3])
%yticks([0.5 1 1.5 2])



figure(1)
title('Optimal bounds','Interpreter','latex','FontSize',16); sgtitle('');
legend('Position',[0.5,0.79,0.12,0.05],'Interpreter','latex','NumColumns',2)
xlabel('$c_A$','Interpreter','latex'); 
ylabel('$c_B$','Interpreter','latex'); 
set(gca,'FontWeight','normal','FontSize',16,'FontName','Helvetica','TickLabelInterpreter','latex','TitleFontWeight','normal','TitleFontSize',1);
set(gcf,'color','w');
% originally: set(gcf,'Position',[100 100 700 500])
set(gcf,'Position',[100 100 350 250])
h = findall(gcf,'Type','Line');
for i = 1:numel(h), uistack(h(i),'top'); end

exportgraphics(gcf, 'cstr_phase_opt.pdf', 'ContentType', 'vector')


% set(gca, 'FontName', 'Latin Modern Roman')
% xlabel('$x$', 'FontName', 'Helvetica', 'Interpreter', 'latex', 'FontSize', 15); 
% ylabel('$f(x)$', 'Interpreter', 'latex', 'FontSize', 15); ax = gca; ax.FontSize = 14;
%xticks([1 1.5 2 2.5 3])
%yticks([0.5 1 1.5 2])


%
% settings for the time-domin plots
%
% gammas =
%   162.1770  173.3966  185.7223  183.9826  173.2991  149.6120  122.2163   32.0720
%    81.2464   74.7939   88.0815   73.7187   64.7745   61.6786   44.3390   19.0168
% delta_bar = 0.01
%    
% delta_aug = 4
% gamma_aug = 1.2
% N = 8
% KRR data = 400
% Bounds data = 200
%

figure(1)
sgtitle('Sub-optimal bounds','Interpreter','latex','FontSize',16); 
subplot(1,2,1)
title(''); 
xlabel('Time step','Interpreter','latex'); ylabel('$c_A$','Interpreter','latex'); 
ylim([1.3 3.1]); xticks(0:8); yticks([1.5 2 2.5 3]); set(gca, 'TickLabelInterpreter', 'latex');
set(gca,'FontWeight','normal','FontSize',16,'FontName','Helvetica'); set(gca,'TitleFontWeight','normal','TitleFontSize',1);
legend('Location','southwest', 'NumColumns',2); grid on
subplot(1,2,2)
title(''); 
xlabel('Time step','Interpreter','latex'); ylabel('$c_B$','Interpreter','latex'); 
ylim([0.4 1.6]); xticks(0:8); yticks([0.4 0.7 1.1 1.4]); set(gca, 'TickLabelInterpreter', 'latex');
set(gca,'FontWeight','normal','FontSize',16,'FontName','Helvetica'); set(gca,'TitleFontWeight','normal','TitleFontSize',1);
legend('Location','northeast','NumColumns',2)
set(gcf,'color','w'); grid on
set(gcf,'Position',[10 10 780 250])
exportgraphics(gcf, 'cstr_time_subopt.pdf', 'ContentType', 'vector')


figure(2)
sgtitle('Optimal bounds','Interpreter','latex','FontSize',16); 
subplot(1,2,1)
title(''); 
xlabel('Time step','Interpreter','latex'); ylabel('$c_A$','Interpreter','latex'); 
ylim([1.3 3.1]); xticks(0:8); yticks([1.5 2 2.5 3]); set(gca, 'TickLabelInterpreter', 'latex');
set(gca,'FontWeight','normal','FontSize',16,'FontName','Helvetica'); set(gca,'TitleFontWeight','normal','TitleFontSize',1);
legend('Location','southwest', 'NumColumns',2); grid on
subplot(1,2,2)
title(''); 
xlabel('Time step','Interpreter','latex'); ylabel('$c_B$','Interpreter','latex'); 
ylim([0.4 1.6]); xticks(0:8); yticks([0.4 0.7 1.1 1.4]); set(gca, 'TickLabelInterpreter', 'latex');
set(gca,'FontWeight','normal','FontSize',16,'FontName','Helvetica'); set(gca,'TitleFontWeight','normal','TitleFontSize',1);
legend('Location','northeast', 'NumColumns',2)
set(gcf,'color','w'); grid on
set(gcf,'Position',[100 100 780 250])
exportgraphics(gcf, 'cstr_time_opt.pdf', 'ContentType', 'vector')
