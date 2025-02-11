% Script for creating all individual graphs for "A dysregulated methylome modulates the radiosensitivity of cancer" article
% Panel figures are then created in Affinity Designer from the individual graphs.
% Anders Berglund 2024
%%

% Please change this acordingle to where you have downloaded the data
BaseDir= "/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/";
ResultDir "RESULTS"

PanelFigDir =

%% Figure 1

% Figure 1a

% load pre calculated PCA model and Groups 
load(fullfile(BaseDir,ResultDir,"PanCan_Results","PCA_ALL_50.mat"))
load(fullfile(BaseDir,ResultDir,"PanCan_Results","Group.mat"))

%Calculate t-SNE model
rng('default')
Ytsne = tsne(PCA_ALL_50.T,'Algorithm','exact','NumPCAComponents',0,'Perplexity',30);
CMap = [GetPalette('Tab20',[2 1]);[ 0 0 0];GetPalette('Tab20',[4 3 6 5 8 7 10 9 14 13 16 15 18 17 20 19])];
nx=1;ny=2;
fh=figure('Name','Scatter Plot','Color','w','Tag','Bar Plot','Units','inches');
fh.Position(3:4) = [3.8 2.7];
ah = axes(fh,'NextPlot','add','tag','Scatter Plot','box','on','Layer','top','FontSize',7);
axis square
ah.LineWidth = 0.5;
gh=gscatter(ah,Ytsne(:,nx),Ytsne(:,ny),Group,CMap,'o',4);
nudgeX=range(Ytsne(:,nx))/50;
ah.XLim = [min(Ytsne(:,nx))-nudgeX  max(Ytsne(:,nx)) + nudgeX];
nudgeY=range(Ytsne(:,ny))/50;
ah.YLim = [min(Ytsne(:,ny))-nudgeY  max(Ytsne(:,ny)) + nudgeY];
legend('off')
% create legend with thicker lines
hCopy = copyobj(fh.Children(1).Children, ah); 
for i=1:19;set(hCopy(i),'XData', NaN', 'YData', NaN);end
for i=1:19;hCopy(i).LineWidth = 1;end
legend(flipud(hCopy))
fh.Children(1).Location='eastoutside';
fh.Children(1).Box='off';
ah.XTick=[];
ah.YTick=[];
ah.XLabel.String = 't-SNE 1';
ah.YLabel.String = 't-SNE 2';
ah.Position=[0.0340    0.0869    0.7222    0.9035];
fh.Children(1).Position = [0.6702    0.1997    0.3412    0.8015];
exportgraphics(gcf,'tSNE_All_Samples.pdf')
exportgraphics(gcf,'tSNE_All_Samples.png','Resolution',600)

% Figure 1b Volcano Plots
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_HNSC_Radiation/RESULTS_2023/RESULTS_HNSC_M450_RT_HPV_Neg_187_DSS')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_HNSC_Radiation/RESULTS_2023/RESULTS_HNSC_M450_NoRT_HPV_Neg_167_DSS.mat')

load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_PRAD_Radiation/RESULTS_2023/RESULTS_PRAD_M450_NoRT_411_PFI.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_PRAD_Radiation/RESULTS_2023/RESULTS_PRAD_M450_RT_59_PFI.mat')

load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_SKCM_Radiation/RESULTS_2023/RESULTS_SKCM_M450_NoRT_351_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_SKCM_Radiation/RESULTS_2023/RESULTS_SKCM_M450_RT_67_DSS.mat')

load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_BRCA_Radiation/RESULTS_2023/RESULTS_BRCA_M450_NoRT_362_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_BRCA_Radiation/RESULTS_2023/RESULTS_BRCA_M450_RT_282_DSS.mat')

PlotSize = [1.8 1.9];

% Create Volcano Plots, highligthing top 1% of the CpG-probes

fh = VolcanoPlotResults(RESULTS_HNSC_M450_RT_HPV_Neg_187_DSS,'HR coxreg DSS',0,'p coxreg DSS',2,'TopPrctile',99,'FigureSize',PlotSize,'FontSize',7,'EqualXLim');
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(fh.Children,'HNSC_M450_RT_HPV_Neg_187_DSS_Volcano_Plot.pdf')
exportgraphics(fh.Children,'HNSC_M450_RT_HPV_Neg_187_DSS_Volcano_Plot.png','Resolution',600)

fh = VolcanoPlotResults(RESULTS_HNSC_M450_NoRT_HPV_Neg_167_DSS,'HR coxreg DSS',0,'p coxreg DSS',2,'TopPrctile',99,'FigureSize',PlotSize,PlotSize,'FontSize',7,'EqualXLim' );
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(fh.Children,'HNSC_M450_NoRT_HPV_Neg_167_DSS_Volcano_Plot.pdf')
exportgraphics(fh.Children,'HNSC_M450_NoRT_HPV_Neg_167_DSS_Volcano_Plot.png','Resolution',600)

fh = VolcanoPlotResults(RESULTS_PRAD_M450_RT_59_PFI,'HR coxreg PFI',0,'p coxreg PFI',2,'TopPrctile',99,'FigureSize',PlotSize,'FontSize',7,'EqualXLim');
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(fh.Children,'PRAD_M450_RT_59_PFI_Volcano_Plot.pdf')
exportgraphics(fh.Children,'PRAD_M450_RT_59_PFI_Volcano_Plot.png','Resolution',600)

fh = VolcanoPlotResults(RESULTS_PRAD_M450_NoRT_411_PFI,'HR coxreg PFI',0,'p coxreg PFI',2,'TopPrctile',99,'FigureSize',PlotSize,PlotSize,'FontSize',7,'EqualXLim');
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(fh.Children,'PRAD_M450_NoRT_411_PFI_Volcano_Plot.pdf')
exportgraphics(fh.Children,'PRAD_M450_NoRT_411_PFI_Volcano_Plot.png','Resolution',600)


fh = VolcanoPlotResults(RESULTS_SKCM_M450_RT_67_DSS,'HR coxreg DSS',0,'p coxreg DSS',2,'TopPrctile',99,'FigureSize',PlotSize,'FontSize',7,'EqualXLim');
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(fh.Children,'SKCM_M450_RT_67_DSS_Volcano_Plot.pdf')
exportgraphics(fh.Children,'SKCM_M450_RT_67_DSS_Volcano_Plot.png','Resolution',600)


fh = VolcanoPlotResults(RESULTS_SKCM_M450_NoRT_351_DSS,'HR coxreg DSS',0,'p coxreg DSS',2,'TopPrctile',99,'FigureSize',PlotSize,'FontSize',7,'EqualXLim');
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(fh.Children,'SKCM_M450_NoRT_351_DSS_Volcano_Plot.pdf')
exportgraphics(fh.Children,'SKCM_M450_NoRT_351_DSS_Volcano_Plot.png','Resolution',600)


fh = VolcanoPlotResults(RESULTS_BRCA_M450_RT_282_DSS,'HR coxreg DSS',0,'p coxreg DSS',2,'TopPrctile',99,'FigureSize',PlotSize,'FontSize',7,'EqualXLim');
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(fh.Children,'BRCA_M450_RT_282_DSS_Volcano_Plot.pdf')
exportgraphics(fh.Children,'BRCA_M450_RT_282_DSS_Volcano_Plot.png','Resolution',600)

fh = VolcanoPlotResults(RESULTS_BRCA_M450_NoRT_362_DSS,'HR coxreg DSS',0,'p coxreg DSS',2,'TopPrctile',99,'FigureSize',PlotSize,'FontSize',7,'EqualXLim');
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(fh.Children,'BRCA_M450_NoRT_362_DSS_Volcano_Plot.pdf')
exportgraphics(fh.Children,'BRCA_M450_NoRT_362_DSS_Volcano_Plot.png','Resolution',600)

% Figure 1c 
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_HNSC_Radiation/RESULTS_2023/RESULTS_HNSC_M450_RT_HPV_Pos_51_DSS.mat')
PlotSize = [1.8 1.9];
fh = VolcanoPlotResults(RESULTS_HNSC_M450_RT_HPV_Pos_51_DSS,'HR coxreg DSS',0,'p coxreg DSS',2,'TopPrctile',99,'FigureSize',PlotSize,'FontSize',7,'EqualXLim' );
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(fh.Children,'HNSC_M450_RT_HPV_Pos_51_DSS_Volcano_Plot.pdf')
exportgraphics(fh.Children,'HNSC_M450_RT_HPV_Pos_51_DSS_Volcano_Plot.png','Resolution',600)

% Figure 1d Density scatter plots comparing results across different cohorts
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_HNSC_Radiation/RESULTS_2023/RESULTS_HNSC_M450_RT_HPV_Neg_187_DSS')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_HNSC_Radiation/RESULTS_2023/RESULTS_HNSC_M450_NoRT_HPV_Neg_167_DSS.mat')

x= -log10(RESULTS_HNSC_M450_RT_HPV_Neg_187_DSS.X(:,5));
y = -log10(RESULTS_HNSC_M450_NoRT_HPV_Neg_167_DSS.X(:,5));
fh=figure('Name','Scatter Plot','Color','w','Tag','Bar Plot','Units','inches');
fh.Position(3:4) = [2 2.1];
ah = axes(fh,'NextPlot','add','tag','Scatter Plot','box','on','Layer','top','FontSize',6);
ah.LineWidth = 0.5;
DensScat(x,y,'ColorMap',colorcet('L08'),'TargetAxes',ah,'mSize',10);
ah.XLim=[0 max(x)+0.2];
ah.YLim=[0 max(y)+0.2];
xlabel('HPV(−)HNSCC RT -log_1_0(p)')
ylabel('HPV(−)HNSCC NoRT -log_1_0(p)')
[r p]=corr(x,y,'Type','Pearson','Rows','pairwise');
text(max(x),max(y),sprintf('r=%.2f',r),'HorizontalAlignment','right','VerticalAlignment','top','FontSize',6)
fh.Renderer='painters';
exportgraphics(fh,'HNSC_RT_vs_NoRT_p_values.pdf')
exportgraphics(fh,'HNSC_RT_vs_NoRT_p_values.png','Resolution',300)

% Figure 1e LGG
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_LGG_Radiation/RESULTS_2023/RESULTS_LGG_M450_NoRT_185_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_LGG_Radiation/RESULTS_2023/RESULTS_LGG_M450_RT_269_DSS.mat')

x = -log10(RESULTS_LGG_M450_RT_269_DSS.X(:,5));
y = -log10(RESULTS_LGG_M450_NoRT_185_DSS.X(:,5));
fh=figure('Name','Scatter Plot','Color','w','Tag','Bar Plot','Units','inches');
fh.Position(3:4) = [2 2];
ah = axes(fh,'NextPlot','add','tag','Scatter Plot','box','on','Layer','top','FontSize',6);
ah.LineWidth = 0.5;
DensScat(x,y,'ColorMap',colorcet('L08'),'TargetAxes',ah,'mSize',10);
ah.XLim=[0 max(x)+0.2];
ah.YLim=[0 max(y)+0.2];
xlabel('LGG RT -log_1_0(p)')
ylabel('LGG NoRT -log_1_0(p)')
[r p]=corr(x,y,'Type','Pearson','Rows','pairwise');
text(max(x),max(y),sprintf('r=%.2f',r),'HorizontalAlignment','right','VerticalAlignment','top','FontSize',6)
fh.Renderer='painters';
exportgraphics(fh,'LGG_RT_vs_NoRT_p_values.pdf')
exportgraphics(fh,'LGG_RT_vs_NoRT_p_values.png','Resolution',300)



% Figure 1f Comparing degree of methylation across all cohorts
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_PRAD_Radiation/RESULTS_2023/RESULTS_PRAD_M450_NoRT_411_PFI.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_PRAD_Radiation/RESULTS_2023/RESULTS_PRAD_M450_RT_59_PFI.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_HNSC_Radiation/RESULTS_2023/RESULTS_HNSC_M450_RT_HPV_Neg_187_DSS')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_HNSC_Radiation/RESULTS_2023/RESULTS_HNSC_M450_NoRT_HPV_Neg_167_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_SKCM_Radiation/RESULTS_2023/RESULTS_SKCM_M450_NoRT_351_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_SKCM_Radiation/RESULTS_2023/RESULTS_SKCM_M450_RT_67_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_BRCA_Radiation/RESULTS_2023/RESULTS_BRCA_M450_NoRT_362_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_BRCA_Radiation/RESULTS_2023/RESULTS_BRCA_M450_RT_282_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_SARC_Radiation/RESULTS_2023/RESULTS_SARC_M450_NoRT_178_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_SARC_Radiation/RESULTS_2023/RESULTS_SARC_M450_RT_59_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_STAD_Radiation/RESULTS_2023/RESULTS_STAD_M450_NoRT_312_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_STAD_Radiation/RESULTS_2023/RESULTS_STAD_M450_RT_60_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_CESC_Radiation/RESULTS_2023/RESULTS_CESC_M450_NoRT_120_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_CESC_Radiation/RESULTS_2023/RESULTS_CESC_M450_RT_134_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_LGG_Radiation/RESULTS_2023/RESULTS_LGG_M450_NoRT_185_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_LGG_Radiation/RESULTS_2023/RESULTS_LGG_M450_RT_269_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_GBM_Radiation/RESULTS_2023/RESULTS_GBM_M450_NoRT_32_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_GBM_Radiation/RESULTS_2023/RESULTS_GBM_M450_RT_80_DSS.mat')

X = zeros(9,2);
[X(1,1)] = CalcVolcanoStatResults(RESULTS_HNSC_M450_RT_HPV_Neg_187_DSS,'HR coxreg DSS','p coxreg DSS');
[X(1,2)] = CalcVolcanoStatResults(RESULTS_HNSC_M450_NoRT_HPV_Neg_167_DSS,'HR coxreg DSS','p coxreg DSS');
[X(2,1)] = CalcVolcanoStatResults(RESULTS_PRAD_M450_RT_59_PFI,'HR coxreg PFI','p coxreg PFI');
[X(2,2)] = CalcVolcanoStatResults(RESULTS_PRAD_M450_NoRT_411_PFI,'HR coxreg PFI','p coxreg PFI');
[X(3,1)] = CalcVolcanoStatResults(RESULTS_SKCM_M450_RT_67_DSS,'HR coxreg DSS','p coxreg DSS');
[X(3,2)] = CalcVolcanoStatResults(RESULTS_SKCM_M450_NoRT_351_DSS,'HR coxreg DSS','p coxreg DSS');
[X(4,1)] = CalcVolcanoStatResults(RESULTS_BRCA_M450_RT_282_DSS,'HR coxreg DSS','p coxreg DSS');
[X(4,2)] = CalcVolcanoStatResults(RESULTS_BRCA_M450_NoRT_362_DSS,'HR coxreg DSS','p coxreg DSS');
[X(5,1)] = CalcVolcanoStatResults(RESULTS_SARC_M450_RT_59_DSS,'HR coxreg DSS','p coxreg DSS');
[X(5,2)] = CalcVolcanoStatResults(RESULTS_SARC_M450_NoRT_178_DSS,'HR coxreg DSS','p coxreg DSS');
[X(6,1)] = CalcVolcanoStatResults(RESULTS_STAD_M450_RT_60_DSS,'HR coxreg DSS','p coxreg DSS');
[X(6,2)] = CalcVolcanoStatResults(RESULTS_STAD_M450_NoRT_312_DSS,'HR coxreg DSS','p coxreg DSS');
[X(7,1)] = CalcVolcanoStatResults(RESULTS_CESC_M450_RT_134_DSS,'HR coxreg DSS','p coxreg DSS');
[X(7,2)] = CalcVolcanoStatResults(RESULTS_CESC_M450_NoRT_120_DSS,'HR coxreg DSS','p coxreg DSS');
[X(8,1)] = CalcVolcanoStatResults(RESULTS_GBM_M450_RT_80_DSS,'HR coxreg DSS','p coxreg DSS');
[X(8,2)] = CalcVolcanoStatResults(RESULTS_GBM_M450_NoRT_32_DSS,'HR coxreg DSS','p coxreg DSS');
[X(9,1)] = CalcVolcanoStatResults(RESULTS_LGG_M450_RT_269_DSS,'HR coxreg DSS','p coxreg DSS');
[X(9,2)] = CalcVolcanoStatResults(RESULTS_LGG_M450_NoRT_185_DSS,'HR coxreg DSS','p coxreg DSS');

TCGA_Ids={'HPV(−)HNSCC','PRAD','SKCM','BRCA','SARC','STAD','CESC','GBM','LGG'};
fh=figure('Name','Bar Plot','Color','w','Tag','Bar Plot','Units','inches');
fh.Position(3:4) = [4,2];
ah = axes(fh,'NextPlot','add','tag','Volcano Plot','box','off','Layer','top','FontSize',7);
ah.LineWidth = 0.5;
sh1=scatter(1:9,X(:,1),100,GetPalette('Tab20',[1 3 5 7 9 13 15 17 19]),'MarkerFaceColor','flat','MarkerEdgeColor',[0.1 0.1 0.1]);
sh2=scatter(1:9,X(:,2),100,GetPalette('Tab20',[2 4 6 8 10 14 16 18 20]),'Linewidth',1);
ylabel('(SumHyper - SumHypo) /SumTotal'); % Greek sum symbol is added in Affinity due to a bug in MATLAB while exporting fonts
ah.XLim=[0.5 9.5];
ah.XTick=1:9;
ah.XTickLabel=TCGA_Ids;
ah.XTickLabelRotation=-45;
line(ah.XLim,[0 0],'Linewidth',1,'Color','k','LineStyle','-.')
line(1:9,X(:,1),'Linewidth',1,'Color','k')
scatter(8.2,0.9,100,[0.3 0.3 0.3],'MarkerFaceColor','flat','MarkerEdgeColor',[0.1 0.1 0.1]);
text(8.5,0.9,'RT','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',7)
scatter(8.2,0.65,100,[0.3 0.3 0.3],'MarkerEdgeColor',[0.1 0.1 0.1],'Linewidth',1);
text(8.5,0.65,'NoRT','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',7)

exportgraphics(fh.Children,'DotPlot_RT_vs_NoRT.pdf')
exportgraphics(fh.Children,'DotPlot_RT_vs_NoRT.png','Resolution',600)


% KM plots
% Figure 1g
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_HNSC_Radiation/RESULTS_2023/AverageMethylation/HyperHypo_HNSC_M450_RT_HPV_Neg_187.mat')
[p_LR,fH,stats] = MatSurv(HyperHypo_HNSC_M450_RT_HPV_Neg_187.SURVIVAL.SurvTime(:,3),HyperHypo_HNSC_M450_RT_HPV_Neg_187.SURVIVAL.SurvEvent(:,3),HyperHypo_HNSC_M450_RT_HPV_Neg_187.X(:,2),'cutpoint','median',...
'Timeunit','Months','Print',1,'RT_KMplot',1,'BaseFontSize',6,'XStep',24,'LineWidth',0.75,'XLim',[0 120],'legend',false,'ylabel','Survival Probability',...
'XTickFontSize',0,'YTickFontSize',0,'LegendFontSize',0,'PvalFontSize',0,'CensorLineWidth',0.5,'LineColor',GetPalette('Lancet',[2 1]),'xlabel','DSS by Hyper-methylation (Months)');
%fH.Children(1).String = {'High Hyper','Low hyper'};
fH.Children(1).Children(1).String = {'Low'};
fH.Children(1).Children(2).String = {'High'};
fH.Children(1).Children(1).FontWeight='normal';
fH.Children(1).Children(2).FontWeight='normal';
fH.Children(1).LineWidth=0.5;
fH.Units = 'inches';
fH.Position(3:4)= [2.1 2.2];
fH.Children(1).Position=[0.1747 0.1549 0.7892 0.8210];
exportgraphics(fH,'Hyper_HNSC_M450_RT_HPV_Neg_187_DSS_KM_Plot.pdf')
exportgraphics(fH,'Hyper_HNSC_M450_RT_HPV_Neg_187_DSS_KM_Plot.png','Resolution',600)

% Figure 1h
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_PRAD_Radiation/RESULTS_2023/AverageMethylation/HyperHypo_RT_59.mat')
[p_LR,fH,stats] = MatSurv(HyperHypo_RT_59.SURVIVAL.SurvTime(:,2),HyperHypo_RT_59.SURVIVAL.SurvEvent(:,2),HyperHypo_RT_59.X(:,2),'cutpoint','median',...
'Timeunit','Months','Print',1,'RT_KMplot',1,'BaseFontSize',6,'XStep',12,'LineWidth',0.75,'XLim',[0 60],'legend',false,'ylabel','Progression Free Probability',...
'XTickFontSize',0,'YTickFontSize',0,'LegendFontSize',0,'PvalFontSize',0,'CensorLineWidth',0.5,'LineColor',GetPalette('Lancet',[2 1]),'xlabel','PFI by Hyper-methylation (Months)');
%fH.Children(1).String = {'High Hyper','Low hyper'};
fH.Children(1).Children(1).String = {'Low'};
fH.Children(1).Children(2).String = {'High'};
fH.Children(1).Children(1).FontWeight='normal';
fH.Children(1).Children(2).FontWeight='normal';
fH.Children(1).LineWidth=0.5;
fH.Units = 'inches';
fH.Position(3:4)= [2.1 2.2];
fH.Children(1).Position=[0.1747 0.1549 0.7892 0.8210];
exportgraphics(fH,'PRAD_M450_RT_59_Hyper_PFI_KM_Plot.pdf')
exportgraphics(fH,'PRAD_M450_RT_59_Hyper_PFI_KM_Plot.png','Resolution',600)




%% Figure 2

% HNSC HPV-
% Figure 2a

load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_HNSC_Radiation/RESULTS_2023/RESULTS_HNSC_M450_RT_HPV_Neg_187_DSS.mat')

fh = ChrPlotDiff(RESULTS_HNSC_M450_RT_HPV_Neg_187_DSS,'chr1',[],'HR coxreg DSS','HR coxreg DSS','p coxreg DSS',[6 2],[1 80],'cytoband','REGION',{[203000000 226100000],4,''});
fh.Children(4).YLim=[-4.5 4.5];
fh.Children(4).CLim = [0 6];
fh.Renderer='painters';
exportgraphics(fh,'Chr_01_Cox_HNSC_RT_HPV_Neg.pdf')
exportgraphics(fh,'Chr_01_Cox_HNSC_RT_HPV_Neg.png','Resolution',600)


%Figure 2b
genes={'TMEM183A','RBBP5','ELK4','RASSF5','DYRK3','MAPKAPK2','CD55','TRAF5','SLC30A1','ATF3','SDE2'};
fh = ChrPlotDiff(RESULTS_HNSC_M450_RT_HPV_Neg_187_DSS,'chr1',[203000000 226100000],'HR coxreg DSS','HR coxreg DSS','p coxreg DSS',[4 1.7],[1 80],'GENES',genes);
fh.Children(2).CLim = [0 6];

% Adjust position of the legend for highlighted genes
h_SDE2 = findobj(fh,'String','SDE2');
h_SDE2.Position = [224.5 2.543 0];

h_TMEM183A = findobj(fh,'String','TMEM183A');
h_TMEM183A.Position = [202.7 3.55 0];

h_RBBP5 = findobj(fh,'String','RBBP5');
h_RBBP5.Position = [205.6 3.35 0];

h_RASSF5 = findobj(fh,'String','RASSF5');
h_RASSF5.Position = [206.9 1.3 0];

h_CD55 = findobj(fh,'String','CD55');
h_CD55.Position = [207.7 1.92 0];

h_MAPKAPK2 = findobj(fh,'String','MAPKAPK2');
h_MAPKAPK2.Position = [206.7092 2.13 0];

h_TRAF5 = findobj(fh,'String','TRAF5');
h_TRAF5.Position = [209.7000 1.3100 0];
exportgraphics(fh,'Chr_01_Cox_HNSC_RT_HPV_Neg_Zoom.pdf')
exportgraphics(fh,'Chr_01_Cox_HNSC_RT_HPV_Neg_Zoom.png','Resolution',600)



%Figure 2c
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_HNSC_Radiation/RESULTS_2023/RESULTS_HNSC_M450_NoRT_HPV_Neg_167_DSS.mat')
fh = ChrPlotDiff(RESULTS_HNSC_M450_NoRT_HPV_Neg_167_DSS,'chr1',[],'HR coxreg DSS','HR coxreg DSS','p coxreg DSS',[6 2],[1 80],'YValCutOff',1,'ColorValCutOff',2,'cytoband');
fh.Children(4).YLim=[-4.5 4.5];
fh.Children(4).CLim = [0 6];
fh.Renderer='painters';
exportgraphics(fh,'Chr_01_Cox_HNSC_NoRT_HPV_Neg.pdf')
exportgraphics(fh,'Chr_01_Cox_HNSC_NoRT_HPV_Neg.png','Resolution',600)


%Figure 2d
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_HNSC_Radiation/RESULTS_2023/AverageMethylation/SURVIVAL_HyperHypo_HNSC_M450_RT_HPV_Neg_187.mat')
samples=SURVIVAL_HyperHypo_HNSC_M450_RT_HPV_Neg_187.RowId(6:end)
fh=Chr_Survival_Dot_Plot(SURVIVAL_HyperHypo_HNSC_M450_RT_HPV_Neg_187,'HR logrank DSS',samples,'p logrank DSS',[1 0.05 0.01],[5 20 50],GetPalette('Tab10',1),3,1.8);
ylabel('Chromosome')
exportgraphics(fh,'Chr_Survival_HNSC_M450_RT_HPV_Neg_187_logrank.pdf')
exportgraphics(fh,'Chr_Survival_HNSC_M450_RT_HPV_Neg_187_logrank.png','Resolution',600)



% PRAD
% Figure 2e
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_PRAD_Radiation/RESULTS_2023/RESULTS_PRAD_M450_RT_59_PFI.mat')
fh = ChrPlotDiff(RESULTS_PRAD_M450_RT_59_PFI,'chr19',[],'HR coxreg PFI','HR coxreg PFI','p coxreg PFI',[6 1.95],[1 80],'cytoband','REGION',{[43500000 47230900],8,''});
fh.Children(4).YLim=[-8.5 8.5];
fh.Children(4).CLim = [0 5];
fh.Renderer='painters';
exportgraphics(fh,'Chr_19_Cox_PRAD_RT.pdf')
exportgraphics(fh,'Chr_19_Cox_PRAD_RT.png','Resolution',600)


% Figure 2f
genes={'ZNF575','CLPTM1','ERCC1','IRF2BP1','BBC3'}
fh = ChrPlotDiff(RESULTS_PRAD_M450_RT_59_PFI,'chr19',[43500000 47230900],'HR coxreg PFI','HR coxreg PFI','p coxreg PFI',[4 1.65],[1 80],'GENES',genes);
fh.Children(2).CLim = [0 5];
% Adjust position of the legend for highlighted genes
fh_BBC3 = findobj(fh,'String','BBC3');
fh_BBC3.Position = [46.83 6.65 0];
fh.Renderer='painters';
exportgraphics(fh,'Chr_19_Cox_PRAD_RT_Zoom.pdf')
exportgraphics(fh,'Chr_19_Cox_PRAD_RT_Zoom.png','Resolution',600)

% Figure 2g
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_PRAD_Radiation/RESULTS_2023/RESULTS_PRAD_M450_NoRT_411_PFI.mat')
fh = ChrPlotDiff(RESULTS_PRAD_M450_NoRT_411_PFI,'chr19',[],'HR coxreg PFI','HR coxreg PFI','p coxreg PFI',[6 1.95],[1 80],'cytoband');
fh.Children(4).YLim=[-3.5 3.5];
fh.Children(4).CLim = [0 5];
fh.Renderer='painters';
exportgraphics(fh,'Chr_19_Cox_PRAD_NoRT.pdf')
exportgraphics(fh,'Chr_19_Cox_PRAD_NoRT.png','Resolution',600)

% Figure 2h
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_PRAD_Radiation/RESULTS_2023/AverageMethylation/SURVIVAL_HyperHypo_RT_59.mat')
samples=SURVIVAL_HyperHypo_RT_59.RowId(6:end)
fh=Chr_Survival_Dot_Plot(SURVIVAL_HyperHypo_RT_59,'HR logrank PFI',samples,'p logrank PFI',[1 0.05 0.01],[5 20 50],GetPalette('Tab10',2),3,1.8);
ylabel('Chromosome')
exportgraphics(fh,'Chr_Survival_PRAD_M450_RT_59_logrank.pdf')
exportgraphics(fh,'Chr_Survival_PRAD_M450_RT_59_logrank.png','Resolution',600)


%SKCM

% Figure 2i
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_SKCM_Radiation/RESULTS_2023/RESULTS_SKCM_M450_RT_67_DSS.mat')
fh = ChrPlotDiff(RESULTS_SKCM_M450_RT_67_DSS,'chr2',[],'HR coxreg DSS','HR coxreg DSS','p coxreg DSS',[6 1.95],[1 80],'cytoband','REGION',{[235500000 242088874],7.5,''});
fh.Children(4).YLim=[-8.3 8.3];
fh.Children(4).CLim = [0 5];
fh.Renderer='painters';
exportgraphics(fh,'Chr_02_Cox_SKCM_RT.pdf')
exportgraphics(fh,'Chr_02_Cox_SKCM_RT.png','Resolution',600)


% Figure 2j
genes={'AGAP1','RBM44','HDAC4','ING5'}
fh = ChrPlotDiff(RESULTS_SKCM_M450_RT_67_DSS,'chr2',[ 235500000 242088874],'HR coxreg DSS','HR coxreg DSS','p coxreg DSS',[4 1.65],[1 80],'GENES',genes);
fh.Children(2).CLim = [0 5];
fh.Renderer='painters';
exportgraphics(fh,'Chr_02_Cox_SKCM_RT_Zoom.pdf')
exportgraphics(fh,'Chr_02_Cox_SKCM_RT_Zoom.png','Resolution',600)

% Figure 2k
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_SKCM_Radiation/RESULTS_2023/RESULTS_SKCM_M450_NoRT_351_DSS.mat')
fh = ChrPlotDiff(RESULTS_SKCM_M450_NoRT_351_DSS,'chr2',[],'HR coxreg DSS','HR coxreg DSS','p coxreg DSS',[6 1.95],[1 80],'cytoband','REGION',{[235500000 242088874],7.5,''});
fh.Children(4).YLim=[-2.5 2.5];
fh.Children(4).CLim = [0 5];
fh.Renderer='painters';
exportgraphics(fh,'Chr_02_Cox_SKCM_NoRT.pdf')
exportgraphics(fh,'Chr_02_Cox_SKCM_NoRT.png','Resolution',600)




%% Figure 3

% Figure 3a
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_PRAD_Radiation/RESULTS_2023/RESULTS_PRAD_M450_RT_59_PFI.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_HNSC_Radiation/RESULTS_2023/RESULTS_HNSC_M450_RT_HPV_Neg_187_DSS')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_SKCM_Radiation/RESULTS_2023/RESULTS_SKCM_M450_RT_67_DSS')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_BRCA_Radiation/RESULTS_2023/RESULTS_BRCA_M450_RT_282_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_SARC_Radiation/RESULTS_2023/RESULTS_SARC_M450_RT_59_DSS')


X = [log2(RESULTS_HNSC_M450_RT_HPV_Neg_187_DSS.X(:,8))> 1 & -log10(RESULTS_HNSC_M450_RT_HPV_Neg_187_DSS.X(:,5))>2 ,...
    log2(RESULTS_PRAD_M450_RT_59_PFI.X(:,8))> 1 & -log10(RESULTS_PRAD_M450_RT_59_PFI.X(:,5))>2 ,...
    log2(RESULTS_SKCM_M450_RT_67_DSS.X(:,8))> 1 & -log10(RESULTS_SKCM_M450_RT_67_DSS.X(:,5))>2 ,...
    log2(RESULTS_BRCA_M450_RT_282_DSS.X(:,8))> 1 & -log10(RESULTS_BRCA_M450_RT_282_DSS.X(:,5))>2 ,...
    log2(RESULTS_SARC_M450_RT_59_DSS.X(:,8))> 1 & -log10(RESULTS_SARC_M450_RT_59_DSS.X(:,5))>2];

s=sum(X,2);
indx=find(s>=2);

% Pie chart
[IDs,numID] = GroupCount(RESULTS_HNSC_M450_RT_HPV_Neg_187_DSS.RowAnnotation(indx,15),1);
fh=figure('Name','Pie Chart','Color','w','Tag','Pie Chart','Units','inches');
fh.Colormap=flipud(fh.Colormap);
IDs=strrep(IDs,'Empty Cell','Missing');
IDs=strrep(IDs,'_',' ');
fh.Position(3:4) = [4 2];
p_indx= [2 6 7 4 5 3 1];
ph = pie(numID(p_indx),[1 0 0 0 0 0 0 ]);
set(findobj(ph,'FontSize',10),'FontSize',7);
lh = legend(IDs(p_indx));
lh.FontSize=7;
lh.Interpreter='none';
set( gca,'FontSize',7);
set(gcf, 'Color', 'w');
lh.Location='eastoutside';
 SavePDF_AEB('CPG_Probe_Type_Pie_Chart')


% Figure 3b GSEA analysis was done using the online tool and the results was then exported to .tsv files
GSEA_REACTOME_Hyper_216 = Read_GSEA_File('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/Combined_Analysis_2023/GSEA_REACTOME_Hyper/GSEA_RECTOME_216_Hyper.tsv')
fh1  = GSEA_Dot_Plot(GSEA_REACTOME_Hyper_216,[],5.2,2.5, [10 20 30],'Description');
SavePDF_AEB('GSEA_REACTOME_Hyper_216')

% Figure 3c
GSEA_TFT_GRD_Hyper_216 = Read_GSEA_File('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/Combined_Analysis_2023/GSEA_TF/TFT_GTRD.tsv')
fh2  = GSEA_Dot_Plot(GSEA_TFT_GRD_Hyper_216,[],4,2.5, [10 20 30],'Name');
SavePDF_AEB('GSEA_GTRD')


% KM plots using the 216 selected CpG-probes

load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_PRAD_Radiation/DATA/PRAD_M450_RT_59.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_HNSC_Radiation/DATA/HNSC_M450_RT_HPV_Neg_187')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_SKCM_Radiation/DATA/SKCM_M450_RT_67')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_BRCA_Radiation/DATA/BRCA_M450_RT_282.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_SARC_Radiation/DATA/SARC_M450_RT_59')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/Combined_Analysis_2023/CpG_probes_216.mat')

HNSC_M450_RT_HPV_Neg_187_216probes =  EditVariablesDATA(HNSC_M450_RT_HPV_Neg_187,CpG_probes_216,'Keep');
SKCM_M450_RT_67_216probes =  EditVariablesDATA(SKCM_M450_RT_67,CpG_probes_216,'Keep');
PRAD_M450_RT_59_216probes =  EditVariablesDATA(PRAD_M450_RT_59,CpG_probes_216,'Keep');
BRCA_M450_RT_282_216probes =  EditVariablesDATA(BRCA_M450_RT_282,CpG_probes_216,'Keep');
SARC_M450_RT_59_216probes =  EditVariablesDATA(SARC_M450_RT_59,CpG_probes_216,'Keep');



Width=1.68;
Hight = 1.6;

% Figure 3d
[p_LR,fH,stats] = MatSurv(HNSC_M450_RT_HPV_Neg_187.SURVIVAL.SurvTime(:,3),HNSC_M450_RT_HPV_Neg_187.SURVIVAL.SurvEvent(:,3),mean(HNSC_M450_RT_HPV_Neg_187_216probes.X,2,'omitnan'),'cutpoint','median',...
'Timeunit','Months','Print',1,'RT_KMplot',1,'BaseFontSize',5,'XStep',24,'LineWidth',0.75,'XLim',[0 120],'legend',false,'ylabel','Survival Probability',...
'XTickFontSize',0,'YTickFontSize',0,'LegendFontSize',0,'PvalFontSize',0,'CensorLineWidth',0.5,'LineColor',GetPalette('Lancet',[2 1]),'xlabel','DSS by RRMS (Months)');
%fH.Children(1).String = {'High Hyper','Low hyper'};
fH.Children(1).Children(1).String = {'Low'};
fH.Children(1).Children(2).String = {'High'};
fH.Children(1).Children(1).FontWeight='normal';
fH.Children(1).Children(2).FontWeight='normal';
fH.Children(1).LineWidth=0.5;
fH.Units = 'inches';
fH.Position(3)=Width;fH.Position(4)=Hight;
fH.Children(1).Position=[0.18 0.18 0.78 0.8];
exportgraphics(fH,'Methylation_Score_HNSC_M450_RT_HPV_Neg_187_DSS_KM_Plot.pdf')
exportgraphics(fH,'Methylation_Score_HNSC_M450_RT_HPV_Neg_187_DSS_KM_Plot.png','Resolution',600)

% Figure 3e
[p_LR,fH,stats] = MatSurv(PRAD_M450_RT_59.SURVIVAL.SurvTime(:,2),PRAD_M450_RT_59.SURVIVAL.SurvEvent(:,2),mean(PRAD_M450_RT_59_216probes.X,2,'omitnan'),'cutpoint','median',...
'Timeunit','Months','Print',1,'RT_KMplot',1,'BaseFontSize',5,'XStep',12,'LineWidth',0.75,'XLim',[0 60],'legend',false,'ylabel','Progression Free Probability',...
'XTickFontSize',0,'YTickFontSize',0,'LegendFontSize',0,'PvalFontSize',0,'CensorLineWidth',0.5,'LineColor',GetPalette('Lancet',[2 1]),'xlabel','PFI by RRMS (Months)');
%fH.Children(1).String = {'High Hyper','Low hyper'};
fH.Children(1).Children(1).String = {'Low'};
fH.Children(1).Children(2).String = {'High'};
fH.Children(1).Children(1).FontWeight='normal';
fH.Children(1).Children(2).FontWeight='normal';
fH.Children(1).LineWidth=0.5;
fH.Units = 'inches';
fH.Position(3)=Width;fH.Position(4)=Hight;
fH.Children(1).Position=[0.18 0.18 0.78 0.8];
exportgraphics(fH,'Methylation_Score_PRAD_M450_RT_59_Hyper_PFI_KM_Plot.pdf')
exportgraphics(fH,'Methylation_Score_PRAD_M450_RT_59_Hyper_PFI_KM_Plot.png','Resolution',600)

% Figure 3f
[p_LR,fH,stats] = MatSurv(SKCM_M450_RT_67_216probes.SURVIVAL.SurvTime(:,3),SKCM_M450_RT_67_216probes.SURVIVAL.SurvEvent(:,3),mean(SKCM_M450_RT_67_216probes.X,2,'omitnan'),'cutpoint','median',...
'Timeunit','Months','Print',1,'RT_KMplot',1,'BaseFontSize',5,'XStep',24,'LineWidth',0.75,'XLim',[0 120],'legend',false,'ylabel','Survival Probability',...
'XTickFontSize',0,'YTickFontSize',0,'LegendFontSize',0,'PvalFontSize',0,'CensorLineWidth',0.5,'LineColor',GetPalette('Lancet',[2 1]),'xlabel','DSS by RRMS (Months)');
%fH.Children(1).String = {'High Hyper','Low hyper'};
fH.Children(1).Children(1).String = {'Low'};
fH.Children(1).Children(2).String = {'High'};
fH.Children(1).Children(1).FontWeight='normal';
fH.Children(1).Children(2).FontWeight='normal';
fH.Children(1).LineWidth=0.5;
fH.Units = 'inches';
fH.Position(3)=Width;fH.Position(4)=Hight;
fH.Children(1).Position=[0.18 0.18 0.78 0.8];
exportgraphics(fH,'Methylation_Score_SKCM_M450_RT_67_DSS_KM_Plot.pdf')
exportgraphics(fH,'Methylation_Score_SKCM_M450_RT_67_DSS_KM_Plot.png','Resolution',600)


% Figure 3g
[p_LR,fH,stats] = MatSurv(BRCA_M450_RT_282_216probes.SURVIVAL.SurvTime(:,3),BRCA_M450_RT_282_216probes.SURVIVAL.SurvEvent(:,3),mean(BRCA_M450_RT_282_216probes.X,2,'omitnan'),'cutpoint','median','TimeMax',120,...
'Timeunit','Months','Print',1,'RT_KMplot',1,'BaseFontSize',5,'XStep',24,'LineWidth',0.75,'XLim',[0 120],'legend',false,'ylabel','Survival Probability',...
'XTickFontSize',0,'YTickFontSize',0,'LegendFontSize',0,'PvalFontSize',0,'CensorLineWidth',0.5,'LineColor',GetPalette('Lancet',[2 1]),'xlabel','DSS by RRMS (Months)');
%fH.Children(1).String = {'High Hyper','Low hyper'};
fH.Children(1).Children(1).String = {'Low'};
fH.Children(1).Children(2).String = {'High'};
fH.Children(1).Children(1).FontWeight='normal';
fH.Children(1).Children(2).FontWeight='normal';
fH.Children(1).LineWidth=0.5;
fH.Units = 'inches';
fH.Position(3)=Width;fH.Position(4)=Hight;
fH.Children(1).Position=[0.18 0.18 0.78 0.8];
exportgraphics(fH,'Methylation_Score_BRCA_M450_RT_282_DSS_KM_Plot.pdf')
exportgraphics(fH,'Methylation_Score_BRCA_M450_RT_282_DSS_KM_Plot.png','Resolution',600)

% Figure 3h
[p_LR,fH,stats] = MatSurv(SARC_M450_RT_59_216probes.SURVIVAL.SurvTime(:,3),SARC_M450_RT_59_216probes.SURVIVAL.SurvEvent(:,3),mean(SARC_M450_RT_59_216probes.X,2,'omitnan'),'cutpoint','median',...
'Timeunit','Months','Print',1,'RT_KMplot',1,'BaseFontSize',5,'XStep',12,'LineWidth',0.75,'XLim',[0 60],'legend',false,'ylabel','Survival Probability',...
'XTickFontSize',0,'YTickFontSize',0,'LegendFontSize',0,'PvalFontSize',0,'CensorLineWidth',0.5,'LineColor',GetPalette('Lancet',[2 1]),'xlabel','DSS by RRMS (Months)');
%fH.Children(1).String = {'High Hyper','Low hyper'};
fH.Children(1).Children(1).String = {'Low'};
fH.Children(1).Children(2).String = {'High'};
fH.Children(1).Children(1).FontWeight='normal';
fH.Children(1).Children(2).FontWeight='normal';
fH.Children(1).LineWidth=0.5;
fH.Units = 'inches';
fH.Position(3)=Width;fH.Position(4)=Hight;
fH.Children(1).Position=[0.18 0.18 0.78 0.8];
exportgraphics(fH,'Methylation_Score_SARC_M450_RT_59_DSS_KM_Plot.pdf')
exportgraphics(fH,'Methylation_Score_SARC_M450_RT_59_DSS_KM_Plot.png','Resolution',600)


%% Figure 4 CCLE data

% Figure 4a
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/CCLE/RESULTS/Melanoma/RESULTS_Melanoma_23_AUC_High_vs_Low.mat')
fh = VolcanoPlotResults(RESULTS_Melanoma_23_AUC_High_vs_Low,'Delta Average',0,'p t-test',2,'TopPrctile',99,'FigureSize',[2.5 2.5]  );

fh.Children.Position=[0.05    0.15    0.90    0.75];
fh.Renderer='painters';
exportgraphics(fh.Children,'Melanoma_23_AUC_High_vs_Low_VolcanoPlot.pdf')
exportgraphics(fh.Children,'Melanoma_23_AUC_High_vs_Low_VolcanoPlot.png','Resolution',600)

% Figure 4b
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/CCLE/RESULTS/Melanoma/RESULTS_Melanoma_23_AUC_Corr.mat')
fh=figure('Name','Volcano Plot','Color','w','Tag','Volcano Plot figure','GraphicsSmoothing','off','Units','inches','Renderer','painters');
fh.Position(3:4) = [3 2.5];
ah = axes(fh,'NextPlot','add','tag','Volcano Plot','box','off','Linewidth',0.5,'YAxisLocation','origin','Layer','top','FontSize',7);
ah.XGrid = 'on';
ah.YGrid = 'on';
uistack(ah,'top');
indx_NS = ~(RESULTS_Melanoma_23_AUC_Corr.X(:,6) < 0.05);
indx_S = (RESULTS_Melanoma_23_AUC_Corr.X(:,6) < 0.05);
h_NS = line(ah,RESULTS_Melanoma_23_AUC_Corr.X(indx_NS,5),RESULTS_Melanoma_23_AUC_Corr.X(indx_NS,9),'LineStyle','none','Marker','.','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerSize',2);
ah.XLim = [-1 1];
ah.YLim = [0 14];
ah.YAxisLocation='left';
box on;
fh = DensScat( RESULTS_Melanoma_23_AUC_Corr.X(indx_S,5),RESULTS_Melanoma_23_AUC_Corr.X(indx_S,9),'TargetAxes',ah);
xlabel('Spearman''s \rho');
ylabel('M-value range');
line([0 0],fh.Children(2).YLim,'LineWidth',1,'Color','k','LineStyle','--');
fh.Children(2).XTick=-1:0.5:1;
fh.Children(2).YTick=0:2:14;
fh.Children(2).Position = [0.10   0.10    0.67    0.87];
exportgraphics(fh,'M-Value_range_vs_Correlation_DensityScatterPlot_CCLE_Melanoma.pdf')
exportgraphics(fh,'M-Value_range_vs_Correlation_DensityScatterPlot_CCLE_Melanoma.png','Resolution',600)

%Figure 4c

load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/CCLE/DATA/CCLE_M450_Melanoma_23.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/CCLE/RESULTS/Melanoma/AUC_23_Mel.mat')

RESULTS_CCLE_M450_Melanoma_23_Average_Methylation = CalculateDiffMethylation(CCLE_M450_Melanoma_23,[]);

[rb pb] = corr(RESULTS_CCLE_M450_Melanoma_23_Average_Methylation.X, AUC_23_Mel,'Type','Spearman');

FontSize = 7;
LineWidth = 0.5;
GridLines = 'on';

fh=Dot_Plot(rb(2:end),pb(2:end),[1 0.05 0.01],[5 20 40],GetPalette('Tab10',3),3,2.5);
xlabel('Spearman''s \rho');
ylabel('Chromosome');
exportgraphics(fH,'CCLE_Melanoma_Chromosome_Average_Methylation_Corr_AUC.pdf')
SavePDF_AEB('CCLE_Melanoma_Chromosome_Average_Methylation_Corr_AUC')

%Figure 4d

load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/CCLE/RESULTS/Melanoma/RESULTS_Melanoma_23_AUC_Corr.mat')

fh = ChrPlotDiff(RESULTS_Melanoma_23_AUC_Corr,'chr8',[],'r Spearman','p Spearman','Range',[6 2.15],[1 80],'cytoband','SizeLegend',[1 0.05 0.01 0.001],[5 20 40 80],[0.25 0.4 0.55 0.7]);
fh.Children(4).YLim=[-0.88 0.88];
fh.Renderer='painters';
exportgraphics(fh,'CCLE_Melanoma_Chromosome_Average_Methylation_Corr_AUC.pdf')
exportgraphics(fh,'CCLE_Melanoma_Chromosome_Average_Methylation_Corr_AUC.png','Resolution',600)


%Figure 4e
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/CCLE/DATA/CCLE_M450_Melanoma_23.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/CCLE/RESULTS/Melanoma/AUC_23_Mel.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/Combined_Analysis_2023/CpG_probes_216.mat')
CCLE_M450_Melanoma_23_216_probes =  EditVariablesDATA(CCLE_M450_Melanoma_23,CpG_probes_216,'Keep');
x=mean(CCLE_M450_Melanoma_23_216_probes.X,2,'omitnan');
[r p_val] = corr(x, AUC_23_Mel,'Type','Spearman')

fh=figure('Name','Bar Plot','Color','w','Tag','Scatter Plot','Units','inches');
fh.Position(3:4) = [2.2,2.2];
ah = axes(fh,'NextPlot','add','tag','Volcano Plot','box','off','Layer','top','FontSize',7);
ah.XGrid = 'on';
ah.YGrid = 'on';
ah.LineWidth = 0.5;
sh=scatter(ah,x,AUC_23_Mel,50,'o','MarkerFaceColor',GetPalette('Tab10',3),'MarkerEdgeColor',[0.2 0.2 0.2]);
AlphaValue=0.9;
sh.MarkerEdgeAlpha = AlphaValue;
sh.MarkerFaceAlpha = AlphaValue;
box on
xlabel('RRMS')
ylabel({'AUC'});
p = polyfit(x,AUC_23_Mel,1) 
f = polyval(p,x);
hold on
line(x,f,'Color','black');
axis square
text(ah,0.161,5.5,{'\rho=0.57';'p=0.006'},'FontSize',7)
exportgraphics(fh,'AUC_vs_MethylationScore_CCLE_Melanoma.pdf')
exportgraphics(fh,'AUC_vs_MethylationScore_CCLE_Melanoma.png','Resolution',600)



%% Supplementary figure 1
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_SARC_Radiation/RESULTS_2023/RESULTS_SARC_M450_NoRT_178_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_SARC_Radiation/RESULTS_2023/RESULTS_SARC_M450_RT_59_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_STAD_Radiation/RESULTS_2023/RESULTS_STAD_M450_NoRT_312_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_STAD_Radiation/RESULTS_2023/RESULTS_STAD_M450_RT_60_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_CESC_Radiation/RESULTS_2023/RESULTS_CESC_M450_NoRT_120_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_CESC_Radiation/RESULTS_2023/RESULTS_CESC_M450_RT_134_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_LGG_Radiation/RESULTS_2023/RESULTS_LGG_M450_NoRT_185_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_LGG_Radiation/RESULTS_2023/RESULTS_LGG_M450_RT_269_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_GBM_Radiation/RESULTS_2023/RESULTS_GBM_M450_NoRT_32_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_GBM_Radiation/RESULTS_2023/RESULTS_GBM_M450_RT_80_DSS.mat')



PlotSize = [1.8 1.9];
%SARC
fh = VolcanoPlotResults(RESULTS_SARC_M450_RT_59_DSS,'HR coxreg DSS',0,'p coxreg DSS',2,'TopPrctile',99,'FigureSize',PlotSize,'FontSize',7,'EqualXLim');
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(fh.Children,'SARC_M450_RT_59_DSS_Volcano_Plot.pdf')
exportgraphics(fh.Children,'SARC_M450_RT_59_DSS_Volcano_Plot.png','Resolution',600)

fh = VolcanoPlotResults(RESULTS_SARC_M450_NoRT_178_DSS,'HR coxreg DSS',0,'p coxreg DSS',2,'TopPrctile',99,'FigureSize',PlotSize,'FontSize',7,'EqualXLim');
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(fh.Children,'SARC_M450_NoRT_178_DSS_Volcano_Plot.pdf')
exportgraphics(fh.Children,'SARC_M450_NoRT_178_DSS_Volcano_Plot.png','Resolution',600)


% STAD
fh = VolcanoPlotResults(RESULTS_STAD_M450_RT_60_DSS,'HR coxreg DSS',0,'p coxreg DSS',2,'TopPrctile',99,'FigureSize',PlotSize,'FontSize',7,'EqualXLim');
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(fh.Children,'STAD_M450_RT_60_DSS_Volcano_Plot.pdf')
exportgraphics(fh.Children,'STAD_M450_RT_60_DSS_Volcano_Plot.png','Resolution',600)

fh = VolcanoPlotResults(RESULTS_STAD_M450_NoRT_312_DSS,'HR coxreg DSS',0,'p coxreg DSS',2,'TopPrctile',99,'FigureSize',PlotSize,'FontSize',7,'EqualXLim');
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(fh.Children,'STAD_M450_NoRT_312_DSS_Volcano_Plot.pdf')
exportgraphics(fh.Children,'STAD_M450_NoRT_312_DSS_Volcano_Plot.png','Resolution',600)


% CESC
fh = VolcanoPlotResults(RESULTS_CESC_M450_RT_134_DSS,'HR coxreg DSS',0,'p coxreg DSS',2,'TopPrctile',99,'FigureSize',PlotSize,'FontSize',7,'EqualXLim');
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(fh.Children,'CESC_M450_RT_134_DSS_Volcano_Plot.pdf')
exportgraphics(fh.Children,'CESC_M450_RT_134_DSS_Volcano_Plot.png','Resolution',600)


fh = VolcanoPlotResults(RESULTS_CESC_M450_NoRT_120_DSS,'HR coxreg DSS',0,'p coxreg DSS',2,'TopPrctile',99,'FigureSize',PlotSize,'FontSize',7,'EqualXLim');
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(fh.Children,'CESC_M450_NoRT_120_DSS_Volcano_Plot.pdf')
exportgraphics(fh.Children,'CESC_M450_NoRT_120_DSS_Volcano_Plot.png','Resolution',600)


% GBM
fh = VolcanoPlotResults(RESULTS_GBM_M450_RT_80_DSS,'HR coxreg DSS',0,'p coxreg DSS',2,'TopPrctile',99,'FigureSize',PlotSize,'FontSize',7,'EqualXLim');
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(fh.Children,'GBM_M450_RT_80_DSS_Volcano_Plot.pdf')
exportgraphics(fh.Children,'GBM_M450_RT_80_DSS_Volcano_Plot.png','Resolution',600)


fh = VolcanoPlotResults(RESULTS_GBM_M450_NoRT_32_DSS,'HR coxreg DSS',0,'p coxreg DSS',2,'TopPrctile',99,'FigureSize',PlotSize,'FontSize',7,'EqualXLim');
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(fh.Children,'GBM_M450_NoRT_32_DSS_Volcano_Plot.pdf')
exportgraphics(fh.Children,'GBM_M450_NoRT_32_DSS_Volcano_Plot.png','Resolution',600)



% GBM
fh = VolcanoPlotResults(RESULTS_LGG_M450_RT_269_DSS,'HR coxreg DSS',0,'p coxreg DSS',2,'TopPrctile',99,'FigureSize',PlotSize,'FontSize',7,'EqualXLim');
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(fh.Children,'LGG_M450_RT_269_DSS_Volcano_Plot.pdf')
exportgraphics(fh.Children,'LGG_M450_RT_269_DSS_Volcano_Plot.png','Resolution',600)


fh = VolcanoPlotResults(RESULTS_LGG_M450_NoRT_185_DSS,'HR coxreg DSS',0,'p coxreg DSS',2,'TopPrctile',99,'FigureSize',PlotSize,'FontSize',7,'EqualXLim');
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(fh.Children,'LGG_M450_NoRT_185_DSS_Volcano_Plot.pdf')
exportgraphics(fh.Children,'LGG_M450_NoRT_185_DSS_Volcano_Plot.png','Resolution',600)

%% Tables
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_HNSC_Radiation/DATA/HNSC_M450_RT_HPV_Neg_187.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_HNSC_Radiation/DATA/HNSC_M450_NoRT_HPV_Neg_167.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_HNSC_Radiation/DATA/HNSC_M450_RT_HPV_Pos_51.mat')
WriteData(HNSC_M450_RT_HPV_Neg_187,'HNSC_M450_RT_HPV_Neg_187.txt','Survival','SampleAnnotationOnly')
WriteData(HNSC_M450_NoRT_HPV_Neg_167,'HNSC_M450_NoRT_HPV_Neg_167.txt','Survival','SampleAnnotationOnly')
WriteData(HNSC_M450_RT_HPV_Pos_51,'HNSC_M450_RT_HPV_Pos_51.txt','Survival','SampleAnnotationOnly')



load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_PRAD_Radiation/DATA/PRAD_M450_RT_59.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_PRAD_Radiation/DATA/PRAD_M450_NoRT_411.mat')
WriteData(PRAD_M450_RT_59,'PRAD_M450_RT_59.txt','Survival','SampleAnnotationOnly')
WriteData(PRAD_M450_NoRT_411,'PRAD_M450_NoRT_411.txt','Survival','SampleAnnotationOnly')


load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_SKCM_Radiation/DATA/SKCM_M450_RT_67.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_SKCM_Radiation/DATA/SKCM_M450_NoRT_351.mat')
WriteData(SKCM_M450_RT_67,'SKCM_M450_RT_67.txt','Survival','SampleAnnotationOnly')
WriteData(SKCM_M450_NoRT_351,'SKCM_M450_NoRT_351.txt','Survival','SampleAnnotationOnly')
clear

load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_BRCA_Radiation/DATA/BRCA_M450_RT_282.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_BRCA_Radiation/DATA/BRCA_M450_NoRT_362.mat')
WriteData(BRCA_M450_RT_282,'BRCA_M450_RT_282.txt','Survival','SampleAnnotationOnly')
WriteData(BRCA_M450_NoRT_362,'BRCA_M450_NoRT_362.txt','Survival','SampleAnnotationOnly')
clear

load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_SARC_Radiation/DATA/SARC_M450_RT_59.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_SARC_Radiation/DATA/SARC_M450_NoRT_178.mat')
WriteData(SARC_M450_RT_59,'SARC_M450_RT_59.txt','Survival','SampleAnnotationOnly')
WriteData(SARC_M450_NoRT_178,'SARC_M450_NoRT_178.txt','Survival','SampleAnnotationOnly')
clear

load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_STAD_Radiation/DATA/STAD_M450_RT_60.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_STAD_Radiation/DATA/STAD_M450_NoRT_312.mat')
WriteData(STAD_M450_RT_60,'STAD_M450_RT_60.txt','Survival','SampleAnnotationOnly')
WriteData(STAD_M450_NoRT_312,'STAD_M450_NoRT_312.txt','Survival','SampleAnnotationOnly')
clear

load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_CESC_Radiation/DATA/CESC_M450_RT_134.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_CESC_Radiation/DATA/CESC_M450_NoRT_120.mat')
WriteData(CESC_M450_RT_134,'CESC_M450_RT_134.txt','Survival','SampleAnnotationOnly')
WriteData(CESC_M450_NoRT_120,'CESC_M450_NoRT_120.txt','Survival','SampleAnnotationOnly')
clear

load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_GBM_Radiation/DATA/GBM_M450_RT_80.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_GBM_Radiation/DATA/GBM_M450_NoRT_32.mat')
WriteData(GBM_M450_RT_80,'GBM_M450_RT_80.txt','Survival','SampleAnnotationOnly')
WriteData(GBM_M450_NoRT_32,'GBM_M450_NoRT_32.txt','Survival','SampleAnnotationOnly')
clear

load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_LGG_Radiation/DATA/LGG_M450_RT_269.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_LGG_Radiation/DATA/LGG_M450_NoRT_185.mat')
WriteData(LGG_M450_RT_269,'LGG_M450_RT_269.txt','Survival','SampleAnnotationOnly')
WriteData(LGG_M450_NoRT_185,'LGG_M450_NoRT_185.txt','Survival','SampleAnnotationOnly')
clear


load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/CCLE/DATA/CCLE_M450_Melanoma_23.mat')
WriteData(CCLE_M450_Melanoma_23,'CCLE_M450_Melanoma_23.txt','SampleAnnotationOnly')
