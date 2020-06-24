% xlim([0.01 200])
% ylim([1 10*10^3])
%    set(gca,'YScale','log');
%    set(gca,'XScale','log');
%    set(gca,'XMinorTick','on');
%    set(gca,'YMinorTick','on');
    % set(gca,'YTick', [20 50 80]*10^3);
%        set(gca,'XTick',[0.1 1 10 100]);
         set(gca,'ticklength',0.01*get(gca,'ticklength'))
 %axis tight;
 %  xlabel({'\it{\tau_{\it i}}'},'FontSize',36)
%     ylabel({'P(\tau)'},'FontSize',36,'FontWeight','bold')
    
set(findall(gcf,'-property','FontSize'),'FontName','Times New Roman','FontSize',18,'linewidth',6)   
xlabel('\it \bf Node index','Interpreter','tex','FontSize',24)
ylabel('\it \bf t','Interpreter','tex','FontSize',24) 
%clear all
set(gcf, 'PaperPositionMode', 'auto','position', [0, 0, 500, 500]);
%axis square 
%minor_ticks_red
set(legend,'color','none');
set(legend, 'Box', 'off');
hold on;
