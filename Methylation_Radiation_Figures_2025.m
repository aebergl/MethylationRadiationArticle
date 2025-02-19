% Script for creating all individual graphs for "A dysregulated methylome modulates the radiosensitivity of cancer" article
% Panel figures are then created in Affinity Designer from the individual graphs.
% Anders Berglund 2024
%% Directory paths

% Please change this accordingly to where you have downloaded the data
BaseDir                 = "/Users/berglund.anders/Documents/RadiationArticleData";
ResultDir               = "ResultFiles";
AverageMethylationDir   = "AverageMethylation";

PanelFigDir             = "/Users/berglund.anders/Documents/RadiationArticleData/PanelFigures";
%% Figure Settings

png_res = 600;


%% Check that necessary funcions are available

if ~exist('GetPalette','file')
    error("Please download GetPalette.m from https://github.com/aebergl/AEB_COLOR and add it to the path") 
end
if ~exist('VolcanoPlotResults','file')
    error("Please download VolcanoPlotResults.m from https://github.com/aebergl/BioinformaticsAEB and add it to the path") 
end

if ~exist('DensScat','file')
    error("Please download DensScat.m from https://github.com/aebergl/DensScat and add it to the path") 
end
if ~exist('MatSurv','file')
    error("Please download MatSurv.m from https://github.com/aebergl/MatSurv and add it to the path") 
end

%% Figure 1
FigureDir = "Figure_01";
%Create Figure_01 directory
[status, Msg] = mkdir(PanelFigDir,FigureDir);
if ~status
    error('Could not create %s, reason: %s',FigureDir,Msg)
end

% Figure 1a

% load pre calculated PCA model and Groups 
load(fullfile(BaseDir,ResultDir,"PCA_ALL_50.mat"))
load(fullfile(BaseDir,ResultDir,"Group.mat"))

% Update name for HNSC
Group = strrep(Group,'HNSC HPV- NoRT','HPV(-)HNSCC NoRT');
Group = strrep(Group,'HNSC HPV- RT','HPV(-)HNSCC RT');
Group = strrep(Group,'HNSC HPV+ RT','HPV(+)HNSCC RT');

%Calculate t-SNE model
rng('default')
Ytsne = tsne(PCA_ALL_50.T,'Algorithm','exact','NumPCAComponents',0,'Perplexity',30);
% Create figure
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
fh.Children(1).Position = [0.6702    0.175    0.3412    0.8015];
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1a_tSNE_All_Samples.pdf'));
fh.Renderer='painters';
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1a_tSNE_All_Samples.png'),'Resolution',png_res)
close(fh)
clear PCA_ALL_50 Group status Msg CMap Ytsne nx ny fh ah gh nudgeY nudgeX hCopy i 


% Figure 1b Volcano Plots
% Load Result files
load(fullfile(BaseDir,ResultDir,"RESULTS_HNSC_M450_RT_HPV_Neg_187_DSS.mat"))
load(fullfile(BaseDir,ResultDir,"RESULTS_HNSC_M450_NoRT_HPV_Neg_167_DSS.mat"))

load(fullfile(BaseDir,ResultDir,"RESULTS_PRAD_M450_NoRT_411_PFI.mat"))
load(fullfile(BaseDir,ResultDir,"RESULTS_PRAD_M450_RT_59_PFI.mat"))

load(fullfile(BaseDir,ResultDir,"RESULTS_SKCM_M450_NoRT_351_DSS.mat"))
load(fullfile(BaseDir,ResultDir,"RESULTS_SKCM_M450_RT_67_DSS.mat"))

load(fullfile(BaseDir,ResultDir,"RESULTS_BRCA_M450_NoRT_362_DSS.mat"))
load(fullfile(BaseDir,ResultDir,"RESULTS_BRCA_M450_RT_282_DSS.mat"))


PlotSize = [1.8 1.9];

% Create Volcano Plots, highligthing top 1% of the CpG-probes

fh = VolcanoPlotResults(RESULTS_HNSC_M450_RT_HPV_Neg_187_DSS,'HR coxreg DSS',0,'p coxreg DSS',2,'TopPrctile',99,'FigureSize',PlotSize,'FontSize',7,'EqualXLim','XlimCrop');
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1b_HNSC_M450_RT_HPV_Neg_187_DSS_Volcano_Plot.pdf'));
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1b_HNSC_M450_RT_HPV_Neg_187_DSS_Volcano_Plot.png'),'Resolution',png_res)
close(fh)

fh = VolcanoPlotResults(RESULTS_HNSC_M450_NoRT_HPV_Neg_167_DSS,'HR coxreg DSS',0,'p coxreg DSS',2,'TopPrctile',99,'FigureSize',PlotSize,'FontSize',7,'EqualXLim','XlimCrop' );
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1b_HNSC_M450_NoRT_HPV_Neg_167_DSS_Volcano_Plot.pdf'));
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1b_HNSC_M450_NoRT_HPV_Neg_167_DSS_Volcano_Plot.png'),'Resolution',png_res)
close(fh)

fh = VolcanoPlotResults(RESULTS_PRAD_M450_RT_59_PFI,'HR coxreg PFI',0,'p coxreg PFI',2,'TopPrctile',99,'FigureSize',PlotSize,'FontSize',7,'EqualXLim','XlimCrop');
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1b_PRAD_M450_RT_59_PFI_Volcano_Plot.pdf'));
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1b_PRAD_M450_RT_59_PFI_Volcano_Plot.png'),'Resolution',png_res)
close(fh)

fh = VolcanoPlotResults(RESULTS_PRAD_M450_NoRT_411_PFI,'HR coxreg PFI',0,'p coxreg PFI',2,'TopPrctile',99,'FigureSize',PlotSize,'FontSize',7,'EqualXLim','XlimCrop');
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1b_PRAD_M450_NoRT_411_PFI_Volcano_Plot.pdf'));
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1b_PRAD_M450_NoRT_411_PFI_Volcano_Plot.png'),'Resolution',png_res)
close(fh)

fh = VolcanoPlotResults(RESULTS_SKCM_M450_RT_67_DSS,'HR coxreg DSS',0,'p coxreg DSS',2,'TopPrctile',99,'FigureSize',PlotSize,'FontSize',7,'EqualXLim','XlimCrop');
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1b_SKCM_M450_RT_67_DSS_Volcano_Plot.pdf'));
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1b_SKCM_M450_RT_67_DSS_Volcano_Plot.png'),'Resolution',png_res)
close(fh)

fh = VolcanoPlotResults(RESULTS_SKCM_M450_NoRT_351_DSS,'HR coxreg DSS',0,'p coxreg DSS',2,'TopPrctile',99,'FigureSize',PlotSize,'FontSize',7,'EqualXLim','XlimCrop');
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1b_SKCM_M450_NoRT_351_DSS_Volcano_Plot.pdf'));
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1b_SKCM_M450_NoRT_351_DSS_Volcano_Plot.png'),'Resolution',png_res)
close(fh)

fh = VolcanoPlotResults(RESULTS_BRCA_M450_RT_282_DSS,'HR coxreg DSS',0,'p coxreg DSS',2,'TopPrctile',99,'FigureSize',PlotSize,'FontSize',7,'EqualXLim','XlimCrop');
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1b_BRCA_M450_RT_282_DSS_Volcano_Plot.pdf'));
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1b_BRCA_M450_RT_282_DSS_Volcano_Plot.png'),'Resolution',png_res)
close(fh)

fh = VolcanoPlotResults(RESULTS_BRCA_M450_NoRT_362_DSS,'HR coxreg DSS',0,'p coxreg DSS',2,'TopPrctile',99,'FigureSize',PlotSize,'FontSize',7,'EqualXLim','XlimCrop');
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1b_BRCA_M450_NoRT_362_DSS_Volcano_Plot.pdf'));
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1b_BRCA_M450_NoRT_362_DSS_Volcano_Plot.png'),'Resolution',png_res)
close(fh)

% Figure 1c 
load(fullfile(BaseDir,ResultDir,"RESULTS_HNSC_M450_RT_HPV_Pos_51_DSS.mat"))
fh = VolcanoPlotResults(RESULTS_HNSC_M450_RT_HPV_Pos_51_DSS,'HR coxreg DSS',0,'p coxreg DSS',2,'TopPrctile',99,'FigureSize',PlotSize,'FontSize',7,'EqualXLim','XlimCrop' );
fh.Children.Position=[0.02    0.15    0.96    0.78];
fh.Renderer='painters';
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1c_HNSC_M450_RT_HPV_Pos_51_DSS_Volcano_Plot.pdf'));
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1c_HNSC_M450_RT_HPV_Pos_51_DSS_Volcano_Plot.png'),'Resolution',png_res)
close(fh)

% Figure 1d Comparing degree of methylation across all cohorts

load(fullfile(BaseDir,ResultDir,"RESULTS_SARC_M450_NoRT_178_DSS.mat"))
load(fullfile(BaseDir,ResultDir,"RESULTS_SARC_M450_RT_59_DSS.mat"))

load(fullfile(BaseDir,ResultDir,"RESULTS_STAD_M450_NoRT_312_DSS.mat"))
load(fullfile(BaseDir,ResultDir,"RESULTS_STAD_M450_RT_60_DSS.mat"))

load(fullfile(BaseDir,ResultDir,"RESULTS_CESC_M450_NoRT_120_DSS.mat"))
load(fullfile(BaseDir,ResultDir,"RESULTS_CESC_M450_RT_134_DSS.mat"))

load(fullfile(BaseDir,ResultDir,"RESULTS_LGG_M450_NoRT_185_DSS.mat"))
load(fullfile(BaseDir,ResultDir,"RESULTS_LGG_M450_RT_269_DSS.mat"))

load(fullfile(BaseDir,ResultDir,"RESULTS_GBM_M450_NoRT_32_DSS.mat"))
load(fullfile(BaseDir,ResultDir,"RESULTS_GBM_M450_RT_80_DSS.mat"))


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

TCGA_Ids={'HPV(-)HNSCC','PRAD','SKCM','BRCA','SARC','STAD','CESC','GBM','LGG'};
fh=figure('Name','Bar Plot','Color','w','Tag','Bar Plot','Units','inches');
fh.Position(3:4) = [4,2];
ah = axes(fh,'NextPlot','add','tag','Volcano Plot','box','off','Layer','top','FontSize',7);
ah.LineWidth = 0.5;
sh1=scatter(1:9,X(:,1),100,GetPalette('Tab20',[1 3 5 7 9 13 15 17 19]),'MarkerFaceColor','flat','MarkerEdgeColor',[0.1 0.1 0.1]);
sh2=scatter(1:9,X(:,2),100,GetPalette('Tab20',[2 4 6 8 10 14 16 18 20]),'Linewidth',1);
ylabel('(\SigmaHyper - \SigmaHypo) /\SigmaTotal');
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
fh.Renderer='painters';
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1d_DotPlot_RT_vs_NoRT.pdf'));
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1d_DotPlot_RT_vs_NoRT.png'),'Resolution',png_res)
close(fh)

clear X TCGA_Ids fh ah sh1 sh2 

% Figure 1e Density scatter plots comparing HNSC RT vs NoRT
indx = strcmp('p coxreg DSS',RESULTS_HNSC_M450_RT_HPV_Neg_187_DSS.ColId);
x = -log10(RESULTS_HNSC_M450_RT_HPV_Neg_187_DSS.X(:,indx));
y = -log10(RESULTS_HNSC_M450_NoRT_HPV_Neg_167_DSS.X(:,indx));
fh=figure('Name','Scatter Plot','Color','w','Tag','Bar Plot','Units','inches');
fh.Position(3:4) = [2 2.1];
ah = axes(fh,'NextPlot','add','tag','Scatter Plot','box','on','Layer','top','FontSize',6);
ah.LineWidth = 0.5;
DensScat(x,y,'ColorMap',colorcet('L08'),'TargetAxes',ah,'mSize',10);
ah.XLim=[0 max(x)+0.2];
ah.YLim=[0 max(y)+0.2];
xlabel('HPV(−)HNSCC RT -log_1_0(p)')
ylabel('HPV(−)HNSCC NoRT -log_1_0(p)')
[r, ~]=corr(x,y,'Type','Pearson','Rows','pairwise');
text(max(x),max(y),sprintf('r=%.2f',r),'HorizontalAlignment','right','VerticalAlignment','top','FontSize',6)
fh.Renderer='painters';
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1e_HNSC_RT_vs_NoRT_p_values_Density_Scatter.pdf'));
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1e_HNSC_RT_vs_NoRT_p_values_Density_Scatter.png'),'Resolution',png_res)
close(fh)
clear indx x y fh ah r 

% Figure 1f Density scatter plots comparing LGG RT vs NoRT
indx = strcmp('p coxreg DSS',RESULTS_LGG_M450_RT_269_DSS.ColId);
x = -log10(RESULTS_LGG_M450_RT_269_DSS.X(:,indx));
y = -log10(RESULTS_LGG_M450_NoRT_185_DSS.X(:,indx));
fh=figure('Name','Scatter Plot','Color','w','Tag','Bar Plot','Units','inches');
fh.Position(3:4) = [2 2];
ah = axes(fh,'NextPlot','add','tag','Scatter Plot','box','on','Layer','top','FontSize',6);
ah.LineWidth = 0.5;
DensScat(x,y,'ColorMap',colorcet('L08'),'TargetAxes',ah,'mSize',10);
ah.XLim=[0 max(x)+0.2];
ah.YLim=[0 max(y)+0.2];
xlabel('LGG RT -log_1_0(p)')
ylabel('LGG NoRT -log_1_0(p)')
[r, ~]=corr(x,y,'Type','Pearson','Rows','pairwise');
text(max(x),max(y),sprintf('r=%.2f',r),'HorizontalAlignment','right','VerticalAlignment','top','FontSize',6)
fh.Renderer='painters';
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1f_LGG_RT_vs_NoRT_p_values_Density_Scatter.pdf'));
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1f_LGG_RT_vs_NoRT_p_values_Density_Scatter.png'),'Resolution',png_res)
close(fh)
clear indx x y fh ah r 




% KM plots
% Figure 1g
load(fullfile(BaseDir,AverageMethylationDir,"HyperHypo_HNSC_M450_RT_HPV_Neg_187.mat"))
indx_HypMet = strcmp('Hyper methylation',HyperHypo_HNSC_M450_RT_HPV_Neg_187.ColId);
indx_DSS = strcmp('DSS',HyperHypo_HNSC_M450_RT_HPV_Neg_187.SURVIVAL.SurvivalTypes);
[~,fH,~] = MatSurv(HyperHypo_HNSC_M450_RT_HPV_Neg_187.SURVIVAL.SurvTime(:,indx_DSS),HyperHypo_HNSC_M450_RT_HPV_Neg_187.SURVIVAL.SurvEvent(:,indx_DSS),HyperHypo_HNSC_M450_RT_HPV_Neg_187.X(:,indx_HypMet),'cutpoint','median',...
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
fh.Renderer='painters';
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1g_Hyper_HNSC_M450_RT_HPV_Neg_187_DSS_KM_Plot.pdf'));
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1g_Hyper_HNSC_M450_RT_HPV_Neg_187_DSS_KM_Plot.png'),'Resolution',png_res)
close(fH)
clear indx_HypMet indx_DSS fH


% Figure 1h
load(fullfile(BaseDir,AverageMethylationDir,"HyperHypo_RT_59.mat"))
indx_HypMet = strcmp('Hyper methylation',HyperHypo_RT_59.ColId);
indx_PFI = strcmp('PFI',HyperHypo_RT_59.SURVIVAL.SurvivalTypes);


[p_LR,fH,stats] = MatSurv(HyperHypo_RT_59.SURVIVAL.SurvTime(:,indx_PFI),HyperHypo_RT_59.SURVIVAL.SurvEvent(:,indx_PFI),HyperHypo_RT_59.X(:,indx_HypMet),'cutpoint','median',...
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
fh.Renderer='painters';
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1h_Hyper_PRAD_M450_RT_59_PFI_KM_Plot.pdf'));
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_1h_Hyper_PRAD_M450_RT_59_PFI_KM_Plot.png'),'Resolution',png_res)
close(fH)
clear indx_HypMet indx_PFI fH

% Clear all files related to Figure 1

clear FigureDir RESULTS_* HyperHypo_*


%% Figure 2

FigureDir = "Figure_02";
%Create Figure_02 directory
[status, Msg] = mkdir(PanelFigDir,FigureDir);
if ~status
    error('Could not create %s, reason: %s',FigureDir,Msg)
end


% HNSC HPV-
% Figure 2a
load(fullfile(BaseDir,ResultDir,"RESULTS_HNSC_M450_RT_HPV_Neg_187_DSS.mat"))
fh = ChrPlotDiff(RESULTS_HNSC_M450_RT_HPV_Neg_187_DSS,'CpG_beg','CpG_chrm','chr1','HR coxreg DSS','HR coxreg DSS','p coxreg DSS','cytoband','mb','REGION',{[203000000 226100000] [-2 4],''});
fh.Children(4).YLim=[-4.5 4.5];
fh.Children(4).CLim = [0 6];
fh.Renderer='painters';
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_2a_Chr_01_Cox_HNSC_RT_HPV_Neg.pdf'));
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_2a_Chr_01_Cox_HNSC_RT_HPV_Neg.png'),'Resolution',png_res)
close(fh)

%Figure 2b
genes={'TMEM183A','RBBP5','ELK4','RASSF5','DYRK3','MAPKAPK2','CD55','TRAF5','SLC30A1','ATF3','SDE2'};
fh = ChrPlotDiff(RESULTS_HNSC_M450_RT_HPV_Neg_187_DSS,'CpG_beg','CpG_chrm','chr1','HR coxreg DSS','HR coxreg DSS','p coxreg DSS','PlotRange',[203000000 226100000],'FigSize',[4 1.7],'GENES','gene_HGNC',genes);
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

fh.Renderer='painters';
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_2b_Chr_01_Cox_HNSC_RT_HPV_Neg_Zoom.pdf'));
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_2b_Chr_01_Cox_HNSC_RT_HPV_Neg_Zoom.png'),'Resolution',png_res)
close(fh)
clear genes h_* fh


%Figure 2c
load(fullfile(BaseDir,ResultDir,"RESULTS_HNSC_M450_NoRT_HPV_Neg_167_DSS.mat"))

fh = ChrPlotDiff(RESULTS_HNSC_M450_NoRT_HPV_Neg_167_DSS,'CpG_beg','CpG_chrm','chr1','HR coxreg DSS','HR coxreg DSS','p coxreg DSS','cytoband','mb');

fh.Children(4).YLim=[-4.5 4.5];
fh.Children(4).CLim = [0 6];
fh.Renderer='painters';
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_2c_Chr_01_Cox_HNSC_NoRT_HPV_Neg.pdf'));
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_2b_Chr_01_Cox_HNSC_RT_HPV_Neg_Zoom.png'),'Resolution',png_res)
close(fh)

%Figure 2d
load(fullfile(BaseDir,AverageMethylationDir,"SURVIVAL_HyperHypo_HNSC_M450_RT_HPV_Neg_187.mat"))
RowsToUse = SURVIVAL_HyperHypo_HNSC_M450_RT_HPV_Neg_187.RowId(6:end);
fh=Chr_Survival_Dot_Plot(SURVIVAL_HyperHypo_HNSC_M450_RT_HPV_Neg_187,'HR logrank DSS',RowsToUse,'p logrank DSS',[1 0.05 0.01],[5 20 50],GetPalette('Tab10',1),3,1.8);
ylabel('Chromosome')
fh.Renderer='painters';
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_2d_hr_Survival_HNSC_M450_RT_HPV_Neg_187_logrank.pdf'));
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_2d_hr_Survival_HNSC_M450_RT_HPV_Neg_187_logrank.png'),'Resolution',png_res)
close(fh)
clear RowsToUse RESULTS_HNSC_M450_NoRT_HPV_Neg_167_DSS SURVIVAL_HyperHypo_HNSC_M450_RT_HPV_Neg_187 RESULTS_HNSC_M450_RT_HPV_Neg_187_DSS fh



% PRAD
% Figure 2e
load(fullfile(BaseDir,ResultDir,"RESULTS_PRAD_M450_RT_59_PFI.mat"))
fh = ChrPlotDiff(RESULTS_PRAD_M450_RT_59_PFI,'CpG_beg','CpG_chrm','chr19','HR coxreg PFI','HR coxreg PFI','p coxreg PFI','cytoband','mb','FigSize',[6 1.95],'REGION',{[43500000 47230900] [-4 8],''});
fh.Children(4).YLim=[-8.5 8.5];
fh.Children(4).CLim = [0 5];
fh.Renderer='painters';
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_2e_Chr_19_Cox_PRAD_RT.pdf'));
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_2e_Chr_19_Cox_PRAD_RT.png'),'Resolution',png_res)
close(fh)



% Figure 2f
genes={'ZNF575','CLPTM1','ERCC1','IRF2BP1','BBC3'};
fh = ChrPlotDiff(RESULTS_PRAD_M450_RT_59_PFI,'CpG_beg','CpG_chrm','chr19','HR coxreg PFI','HR coxreg PFI','p coxreg PFI','PlotRange',[43500000 47230900],'FigSize',[4 1.65],'GENES','gene_HGNC',genes);

fh.Children(2).CLim = [0 5];
% Adjust position of the legend for highlighted genes
fh_BBC3 = findobj(fh,'String','BBC3');
fh_BBC3.Position = [46.83 6.65 0];
fh.Renderer='painters';
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_2f_Chr_19_Cox_PRAD_RT_Zoom.pdf'));
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_2f_Chr_19_Cox_PRAD_RT_Zoom.png'),'Resolution',png_res)
close(fh)

% Figure 2g
load(fullfile(BaseDir,ResultDir,"RESULTS_PRAD_M450_NoRT_411_PFI.mat"))
fh = ChrPlotDiff(RESULTS_PRAD_M450_NoRT_411_PFI,'CpG_beg','CpG_chrm','chr19','HR coxreg PFI','HR coxreg PFI','p coxreg PFI','FigSize',[6 1.95],'cytoband','mb');
fh.Children(4).YLim=[-3.5 3.5];
fh.Children(4).CLim = [0 5];
fh.Renderer='painters';
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_2g_Chr_19_Cox_PRAD_NoRT.pdf'));
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_2g_Chr_19_Cox_PRAD_NoRT.png'),'Resolution',png_res)
close(fh)


% Figure 2h
load(fullfile(BaseDir,AverageMethylationDir,"SURVIVAL_HyperHypo_RT_59.mat"))
RowsToUse = SURVIVAL_HyperHypo_RT_59.RowId(6:end);

fh=Chr_Survival_Dot_Plot(SURVIVAL_HyperHypo_RT_59,'HR logrank PFI',RowsToUse,'p logrank PFI',[1 0.05 0.01],[5 20 50],GetPalette('Tab10',2),3,1.8);
ylabel('Chromosome')
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_2h_Chr_Survival_PRAD_M450_RT_59_logrank.pdf'));
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_2h_Chr_Survival_PRAD_M450_RT_59_logrank.png'),'Resolution',png_res)
close(fh)
clear RowsToUse RESULTS_PRAD_M450_RT_59_PFI RESULTS_PRAD_M450_NoRT_411_PFI SURVIVAL_HyperHypo_RT_59 fh fh_BBC3



%SKCM

% Figure 2i
load(fullfile(BaseDir,ResultDir,"RESULTS_SKCM_M450_RT_67_DSS.mat"))
fh = ChrPlotDiff(RESULTS_SKCM_M450_RT_67_DSS,'CpG_beg','CpG_chrm','chr2','HR coxreg DSS','HR coxreg DSS','p coxreg DSS','cytoband','mb','FigSize',[6 1.95],'REGION',{[235500000 242088874] [-6 8],''});
fh.Children(4).YLim=[-8.3 8.3];
fh.Children(4).CLim = [0 5];
fh.Renderer='painters';
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_2i_Chr_02_Cox_SKCM_RT.pdf'));
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_2i_Chr_02_Cox_SKCM_RT.png'),'Resolution',png_res)
close(fh)


% Figure 2j
genes={'AGAP1','RBM44','HDAC4','ING5'};
fh = ChrPlotDiff(RESULTS_SKCM_M450_RT_67_DSS,'CpG_beg','CpG_chrm','chr2','HR coxreg DSS','HR coxreg DSS','p coxreg DSS','PlotRange',[235500000 242088874],'FigSize',[4 1.65],'GENES','gene_HGNC',genes);
fh.Children(2).YLim=[-6 8];
fh.Children(2).CLim = [0 5];
fh.Renderer='painters';
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_2j_Chr_02_Cox_SKCM_RT_Zoom.pdf'));
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_2j_Chr_02_Cox_SKCM_RT_Zoom.png'),'Resolution',png_res)
close(fh)


% Figure 2k
load(fullfile(BaseDir,ResultDir,"RESULTS_SKCM_M450_NoRT_351_DSS.mat"))
fh = ChrPlotDiff(RESULTS_SKCM_M450_NoRT_351_DSS,'CpG_beg','CpG_chrm','chr2','HR coxreg DSS','HR coxreg DSS','p coxreg DSS','cytoband','mb','FigSize',[6 1.95]);

fh.Children(4).YLim=[-2.5 2.5];
fh.Children(4).CLim = [0 5];
fh.Renderer='painters';
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_2k_Chr_02_Cox_SKCM_NoRT.pdf'));
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_2k_Chr_02_Cox_SKCM_NoRT.png'),'Resolution',png_res)
close(fh)

clear fh RESULTS_SKCM_M450_RT_67_DSS RESULTS_SKCM_M450_NoRT_351_DSS genes FigureDir


%% Figure 3
FigureDir = "Figure_03";
%Create Figure_03 directory
[status, Msg] = mkdir(PanelFigDir,FigureDir);
if ~status
    error('Could not create %s, reason: %s',FigureDir,Msg)
end

% Figure 3a
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_PRAD_Radiation/RESULTS_2023/RESULTS_PRAD_M450_RT_59_PFI.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_HNSC_Radiation/RESULTS_2023/RESULTS_HNSC_M450_RT_HPV_Neg_187_DSS')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_SKCM_Radiation/RESULTS_2023/RESULTS_SKCM_M450_RT_67_DSS')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_BRCA_Radiation/RESULTS_2023/RESULTS_BRCA_M450_RT_282_DSS.mat')
load('/Users/bergluae/AEBERGL/USR/SUNGJUNE/TCGA_Radiation/TCGA_SARC_Radiation/RESULTS_2023/RESULTS_SARC_M450_RT_59_DSS')

load(fullfile(BaseDir,ResultDir,"RESULTS_HNSC_M450_RT_HPV_Neg_187_DSS.mat"))
load(fullfile(BaseDir,ResultDir,"RESULTS_PRAD_M450_RT_59_PFI.mat"))
load(fullfile(BaseDir,ResultDir,"RESULTS_SKCM_M450_RT_67_DSS.mat"))
load(fullfile(BaseDir,ResultDir,"RESULTS_BRCA_M450_RT_282_DSS.mat"))
load(fullfile(BaseDir,ResultDir,"RESULTS_SARC_M450_RT_59_DSS.mat"))


X = [log2(RESULTS_HNSC_M450_RT_HPV_Neg_187_DSS.X(:,8))> 1 & -log10(RESULTS_HNSC_M450_RT_HPV_Neg_187_DSS.X(:,5))>2 ,...
    log2(RESULTS_PRAD_M450_RT_59_PFI.X(:,8))> 1 & -log10(RESULTS_PRAD_M450_RT_59_PFI.X(:,5))>2 ,...
    log2(RESULTS_SKCM_M450_RT_67_DSS.X(:,8))> 1 & -log10(RESULTS_SKCM_M450_RT_67_DSS.X(:,5))>2 ,...
    log2(RESULTS_BRCA_M450_RT_282_DSS.X(:,8))> 1 & -log10(RESULTS_BRCA_M450_RT_282_DSS.X(:,5))>2 ,...
    log2(RESULTS_SARC_M450_RT_59_DSS.X(:,8))> 1 & -log10(RESULTS_SARC_M450_RT_59_DSS.X(:,5))>2];

s = sum(X,2);
indx = find(s>=2);

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
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_3a_CPG_Probe_Type_Pie_Chart.pdf'));
exportgraphics(gcf,fullfile(PanelFigDir,FigureDir,'Figure_3a_CPG_Probe_Type_Pie_Chart.png'),'Resolution',png_res)
close(fh)

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
