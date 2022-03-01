
run('../EssayCodesV8/fig4n6_SetParameter_Pon_Reac_bar_V7.m');

version = '202202';

if ndb
    fig_Reac_bar_filename=strcat('./',version,'figsupp_ndb_Reac_bar_cellcycle_',num2str(cell_cycle),'hours_dTat_',d_Tat_str);
    YL=[-0.5,15];
else
    fig_Reac_bar_filename=strcat('./',version,'fig4_db_Reac_bar_cellcycle_',num2str(cell_cycle),'hours_dTat_',d_Tat_str);
    YL=[-0.5,15];
end

c = categorical(["Untreated" "AC" "NE"...
    "AC+NE" "AC+NS"]);
c = reordercats(c,{'Untreated' 'AC' 'NE' ...
    'AC+NE' 'AC+NS' });

figure(1)
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
ax2=axes;
Reac_bar_100=Reac_bar*100;
b=bar(c,Reac_bar_100 ,'BaseValue',0,'EdgeColor',[0 0 0],'LineWidth',1);
%xlabel()
b.FaceColor = 'flat';
b.CData(1,:) = [0.75 0.75 0.75];
b.CData(2,:) = [0 1 0];
b.CData(3,:) = [1 0 1];
b.CData(4,:) = [1 0 0];
b.CData(5,:) = [0 0 1];

ytickformat(ax2, '%g%%');
set(gca,'YLim',YL,'LineWidth',3)%,'XTick',1:5,...
grid on
set(gca,'fontsize',ax_ticks_fontsize,'fontweight','b','fontname',ax_ticks_fontname);
ylabel('Reactivation','fontsize',xylabel_fontsize,'fontname',xylabel_text_fontname)
saveas(gcf,fig_Reac_bar_filename,'epsc');
close