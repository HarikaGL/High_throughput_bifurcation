% ===================== make_comparison_figures.m ==================================
% Re-creates, in MATLAB, every figure of the comparison report that was originally drawn in
% Python, from the v2 ensemble (data_gen/*_v2_1e6), the original ensemble (for the
% "MATCONT without the inverse method" comparisons) and the ML outputs
% (ML_label_comparison/labels_Exact_Class_v2).
%
% Output: Comparison_figures/<same file names as Comparison/figures/>  (.fig, .pdf and .png)
%         Comparison_figures/figure_numbers.txt   (every number quoted in captions)
% The original Python figures in Comparison/figures are not touched.
% Compile Comparison/comparison_matlab.tex, which reads its figures from figures_matlab/.
%
% Needs base MATLAB R2020a or newer (readtable, exportgraphics, ode15s, fzero).
% No toolbox is required. The exact steady-state analysis and the normal forms are
% implemented below as local functions (same algorithm as ss_structure.m: every
% steady state from the explicit parametrisation in u = K/P).
% Run time: about 20-40 min the first time (exact analysis of ~80,000 parameter sets, the
% normal forms and 500 ode15s integrations); results are cached in figures_matlab/cache_*.mat,
% so later runs (e.g. after changing a plot) take about a minute. Delete the cache files
% (or set USE_CACHE = false) after changing the data.
% =======================================================================================
clear; close all; clc

%% ------------------------------- SETTINGS ---------------------------------------------
ROOT   = fileparts(fileparts(mfilename('fullpath')));   % .../Ensemble_ML (script lives in Comparison/)
if isempty(ROOT) || ~exist(fullfile(ROOT,'data_gen'),'dir')
    ROOT = '/Users/harika/Desktop/Ensemble_ML';
end
V2DIR  = fullfile(ROOT,'data_gen');
OLDDIR = fullfile(ROOT,'Supporting_info','MATCONT_ensemble','Output_Datasets');
MLDIR  = fullfile(ROOT,'ML_label_comparison','labels_Exact_Class_v2');
% Single-start labels to compare with the exact ones:
%   'v3'       : the 40,000 parameter sets of the ten-class ensemble, labelled by
%                SingleStart_v3.m (folder A) -> same parameter sets as Figs 2-5
%   'original' : the 39,646 parameter sets of the original single-start pipeline
SINGLE_START = 'v3';
V3DIR  = fullfile(ROOT,'data_gen','MatCont7p6','Ensemble_v3');
OUT    = fullfile(ROOT,'Comparison_figures');
% Inside the supporting package (folder B): write the figures to 3_results/a_comparison_figures
% and read the ML outputs from 3_results/c_ML_label_comparison.
BDIR   = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if exist(fullfile(BDIR,'3_results','a_comparison_figures'),'dir')
    OUT   = fullfile(BDIR,'3_results','a_comparison_figures');
    MLDIR = fullfile(BDIR,'3_results','c_ML_label_comparison','labels_Exact_Class_v2');
end
W      = 1e6;            % ATOT window of the v2 labels (Settings sheet: [0 1e6])
NCONV  = 125;            % parameter sets per iteration for the ode15s convergence test (Fig SS.1c)
USE_CACHE = true;        % reuse cached exact analysis / normal forms
if ~exist(OUT,'dir'), mkdir(OUT); end
COLS = {'KTOT','PTOT','k1','k11','k2','k3','k33','k4','k5','k55','k6','k7','k77','k8'};
CLS  = {'LP only, S-shaped','LP only, Z-shaped','BP only','LP + BP, x_2^{LP} > x_2^{BP}', ...
        'LP + BP, x_2^{LP} < x_2^{BP}','no singularity, biphasic','LP on detached branch','no singularity, monotone'};

% colours (same palette as the Python figures)
C.ink=[11 11 11]/255; C.ink2=[82 81 78]/255; C.muted=[137 135 129]/255; C.grid=[225 224 217]/255;
C.blue=hex2rgb('#2a78d6'); C.lblue=hex2rgb('#86b6ef'); C.orange=hex2rgb('#eb6834'); C.aqua=hex2rgb('#1baf7a');
C.yellow=hex2rgb('#eda100'); C.red=hex2rgb('#e34948'); C.grey=hex2rgb('#b5b3ab'); C.lgrey=hex2rgb('#cfcdc6');
C.violet=hex2rgb('#4a3aa7'); C.pink=hex2rgb('#e87ba4'); C.green=hex2rgb('#008300'); C.dblue=hex2rgb('#104281');
CCOL=[C.blue;C.orange;C.aqua;C.yellow;C.pink;C.green;C.red;C.grey];        % class colours 1..8
set(groot,'defaultAxesFontName','Helvetica','defaultTextFontName','Helvetica', ...
    'defaultAxesFontSize',8,'defaultTextFontSize',8,'defaultAxesTickDir','out', ...
    'defaultAxesBox','off','defaultAxesLineWidth',0.6,'defaultFigureColor','w');
LOG=fopen(fullfile(OUT,'figure_numbers.txt'),'w');
logf=@(varargin) fprintf_both(LOG,varargin{:});
logf('Figure numbers (MATLAB version), %s\nROOT = %s\nATOT window W = %g\n\n', datestr(now), ROOT, W);

%% ------------------------------- DATA -------------------------------------------------
fprintf('Loading data ...\n');
% v2 ensemble = the ML input (ML_Ready_Dataset_v2_iterN.xlsx; Iteration 3 has 9,969 rows)
V=[]; 
for it=1:4
    f=fullfile(V2DIR,sprintf('Ensemble_iter%d_v2_1e6',it),'Excel_Files',sprintf('ML_Ready_Dataset_v2_iter%d.xlsx',it));
    opts=detectImportOptions(f,'Sheet','ML_Ready_Data');
    opts=setvartype(opts,intersect([COLS,{'Exact_Class_v2','Has_Detached_Branch'}],opts.VariableNames),'double');
    T=readtable(f,opts);
    T.iter=it*ones(height(T),1);
    T=T(:,[COLS,{'Exact_Class_v2','Has_Detached_Branch','iter'}]);
    V=[V;T]; %#ok<AGROW>
end
Pv=V{:,COLS}; itv=V.iter; labv=V.Exact_Class_v2; detv=logical(V.Has_Detached_Branch);
logf('v2 ensemble: %d parameter sets (%s per iteration)\n', size(Pv,1), mat2str(accumarray(itv,1)'));

% original ensemble (labels of MATCONT without the inverse method)
O=[];
for it=1:4
    if strcmp(SINGLE_START,'v3')
        f=fullfile(V3DIR,sprintf('Ensemble_iter%d_v3_1e6',it),'SingleStart','ML_Ready_Dataset_singlestart.csv');
    else
        f=fullfile(OLDDIR,sprintf('Ensemble_iter%d',it),'Excel_Files','ML_Ready_Dataset.csv');
    end
    opts=detectImportOptions(f);
    sv=intersect({'LP_Concentrations','BP_Concentrations','ATOT_LP','ATOT_BP','Multiclass_Label_7'},opts.VariableNames);
    opts=setvartype(opts,sv,'string');
    % force numeric columns to double (a column that is empty in one iteration file is
    % otherwise read as text, and the four tables then cannot be concatenated)
    nv=intersect([COLS,{'Class_Number_7','ATOT_NO_SING','ATOT_LPBP_LP','ATOT_LPBP_BP','ATOT_MAX_BIPHASIC','ATOT_NO_BIPHASIC','CONC_MAX_BIPHASIC','ATOT_reached','Run_seconds','CPU_seconds'}],opts.VariableNames);
    opts=setvartype(opts,nv,'double');
    T=readtable(f,opts);
    T.iter=it*ones(height(T),1);
    T.oldBP=str2double(T.ATOT_BP);
    if ~ismember('ATOT_reached',T.Properties.VariableNames), T.ATOT_reached=nan(height(T),1); end
    T=T(:,[COLS,{'Class_Number_7','ATOT_NO_SING','ATOT_LPBP_LP','oldBP','ATOT_reached','iter'}]);
    O=[O;T]; %#ok<AGROW>
end
O=O(O.Class_Number_7~=9,:);                     % 9 = failed or negative single-start runs
Po=O{:,COLS}; ito=O.iter; oldc=O.Class_Number_7; reach=O.ATOT_NO_SING;
logf('single-start ensemble (%s): %d parameter sets (%s per iteration)\n\n', SINGLE_START, size(Po,1), mat2str(accumarray(ito,1)'));

%% ------------------------------- EXACT ANALYSIS (cached) ------------------------------
cf=fullfile(OUT,sprintf('cache_exact_%s.mat',SINGLE_START));
if USE_CACHE && exist(cf,'file')
    load(cf,'Sv','So');
else
    fprintf('Exact analysis of the v2 ensemble ...\n');      Sv=ex_summary(Pv,'v2');
    fprintf('Exact analysis of the original ensemble ...\n'); So=ex_summary(Po,'original');
    save(cf,'Sv','So');
end
wv=ex_window_all(Sv,W); wo=ex_window_all(So,W);
logf('Agreement of this implementation with Exact_Class_v2: %d of %d sets (%.3f%%)\n\n', ...
    sum(wv==labv), numel(labv), 100*mean(wv==labv));

cf=fullfile(OUT,'cache_normal_forms.mat');
if USE_CACHE && exist(cf,'file')
    load(cf,'NF');
else
    fprintf('Normal forms of the v2 ensemble ...\n'); NF=nf_all(Pv,itv,'normal forms');
    save(cf,'NF');
end

%% =====================================================================================
%% Fig ML.1  ML_fig1_class_by_iteration  -- class composition of each iteration (v2)
%% =====================================================================================
cbi=zeros(8,4); for it=1:4, cbi(:,it)=accumarray(labv(itv==it),1,[8 1]); end
N=sum(cbi,1); Wp=100*cbi./N;
fig=newfig(9,3.9); ax=axes(fig,'Position',[0.17 0.16 0.55 0.78]); hold(ax,'on');
hb=barh(ax,1:4,Wp','stacked','BarWidth',0.66,'EdgeColor','w','LineWidth',1.2);
for k=1:8, hb(k).FaceColor=CCOL(k,:); end
left=zeros(1,4);
for k=1:8
    for i=1:4
        if Wp(k,i)>=6
            tc='w'; if any(k==[4 8]), tc=C.ink; end
            text(ax,left(i)+Wp(k,i)/2,i,sprintf('%d\n%s',k,commas(cbi(k,i))),'HorizontalAlignment','center', ...
                'VerticalAlignment','middle','FontSize',7,'Color',tc);
        end
    end
    left=left+Wp(k,:);
end
lab={'Iteration 1, S-dominant','Iteration 2, Z-dominant','Iteration 3, full symmetry','Iteration 4, binding asymmetry'};
set(ax,'YTick',1:4,'YTickLabel',strcat(lab,{'\newline(n = '},arrayfun(@commas,N,'uni',0),{')'}),'YDir','reverse', ...
    'XLim',[0 100],'TickLength',[0 0]);
xlabel(ax,'share of parameter sets in the iteration (%)');
legend(ax,hb,strcat(arrayfun(@(k) sprintf('%d  ',k),1:8,'uni',0),CLS),'Location','eastoutside','Box','off','FontSize',7);
save_fig(fig,OUT,'ML_fig1_class_by_iteration');
logf('Fig ML.1: class-by-iteration counts (rows class 1-8, columns iteration 1-4):\n'); logf('%s\n\n',mat2str(cbi));

%% =====================================================================================
%% Fig ML.2  ML_fig2_performance  -- LightGBM accuracy by feature set (ML outputs)
%% =====================================================================================
Tp=readtable(fullfile(MLDIR,'ablation_pooled.xlsx')); Ti=readtable(fullfile(MLDIR,'ablation_by_iteration.xlsx'));
Tp=Tp(strcmp(Tp.Model,'LightGBM'),:); Ti=Ti(strcmp(Ti.Model,'LightGBM'),:);
sets='ABCDEF'; G=5; acc=nan(G,6); sd=acc; nfe=acc; base=nan(G,1);
for j=1:6
    r=startsWith(Tp.FeatureSet,[sets(j) ':']); acc(1,j)=Tp.Acc_mean(r); sd(1,j)=Tp.Acc_sd(r); nfe(1,j)=Tp.nFeatures(r); base(1)=Tp.Baseline_TestMajority(r);
    for it=1:4
        r=startsWith(Ti.FeatureSet,[sets(j) ':']) & Ti.Iteration==it;
        acc(it+1,j)=Ti.Acc_mean(r); sd(it+1,j)=Ti.Acc_sd(r); nfe(it+1,j)=Ti.nFeatures(r); base(it+1)=Ti.Baseline_TestMajority(r);
    end
end
setlab={'A  14 params + 91 log-ratios','B  A minus 7 constrained ratios','C  14 log-parameters', ...
        'D  JSD + broken stick (original)','E  CLR of 14 log-parameters','F  CLR of rates + log(K_{TOT}/P_{TOT})'};
scol=[C.blue;C.lblue;C.orange;C.grey;C.aqua;C.green];
fig=newfig(11,5.2); ax=axes(fig,'Position',[0.08 0.38 0.9 0.56]); hold(ax,'on');
hb=bar(ax,1:G,acc,'grouped','BarWidth',0.85,'EdgeColor','none');
for j=1:6
    hb(j).FaceColor=scol(j,:); xe=hb(j).XEndPoints;
    errorbar(ax,xe(:),acc(:,j),sd(:,j),'LineStyle','none','Color',C.ink2,'LineWidth',0.7,'CapSize',0);
    for g=1:G
        tc='w'; if any(j==[2 4]), tc=C.ink; end
        text(ax,xe(g),22.5,fmtn(nfe(g,j)),'Rotation',90,'FontSize',6.3,'Color',tc,'VerticalAlignment','middle');
    end
end
for g=1:G, plot(ax,g+[-0.43 0.43],base(g)*[1 1],'--','Color',C.ink,'LineWidth',1.4); end
plot(ax,[1.5 1.5],[20 103],'-','Color',[0.76 0.76 0.72],'LineWidth',0.8);
grp={'Pooled','Iteration 1','Iteration 2','Iteration 3','Iteration 4'}; sub={'8 classes','classes 1, 7, 8','classes 2, 6, 7','classes 3, 8','classes 4, 5, 8'};
set(ax,'XTick',1:G,'XTickLabel',cellfun(@(a,b,c) sprintf('%s\\newline(%s)\\newlinebaseline %.1f%%',a,b,c),grp,sub,num2cell(base'),'uni',0), ...
    'YLim',[20 103],'XLim',[0.5 G+0.5],'YGrid','on','GridColor',C.grid,'GridAlpha',1);
ylabel(ax,'held-out accuracy (%)');
hl=plot(ax,nan,nan,'--','Color',C.ink,'LineWidth',1.4);
lg=legend(ax,[hl hb],[{'majority-class baseline'},setlab],'NumColumns',3,'Box','off','FontSize',7.5, ...
    'Position',[0.04 0.01 0.92 0.14]); title(lg,'LightGBM, feature set (numbers at the base of each bar: features used)','FontWeight','normal');
save_fig(fig,OUT,'ML_fig2_performance');
logf('Fig ML.2: LightGBM accuracy (rows pooled, iter 1-4; columns sets A-F):\n%s\nfeatures:\n%s\nbaselines %s\n\n',mat2str(acc,4),mat2str(nfe,3),mat2str(base',4));

%% =====================================================================================
%% Fig ML.4  ML_fig4_analytic_conditions  -- top-ranked features encode exact conditions (v2)
%% =====================================================================================
lr=@(a,b) log10(a./b);
X1=lr(Pv(:,14).*Pv(:,2),Pv(:,5).*Pv(:,1)); X2=lr(Pv(:,8).*Pv(:,1),Pv(:,11).*Pv(:,2)); X3=lr(Pv(:,8),Pv(:,5));
% colour = existence of a detached branch at ANY ATOT (exact analysis). The flag Has_Detached_Branch of the
% v2 files counts only detached branches whose fold lies inside the window (ATOT <= 1e6).
pan={1,X1,Sv.has_det,'log_{10}(k_8P_{TOT}/k_2K_{TOT})','detached branch','(a) Iteration 1: detached branch iff k_8P_{TOT} > k_2K_{TOT}',C.red;
     2,X2,Sv.has_det,'log_{10}(k_4K_{TOT}/k_6P_{TOT})','detached branch','(b) Iteration 2: detached branch iff k_4K_{TOT} > k_6P_{TOT}',C.red;
     3,X3,Sv.has_BP,'log_{10}(k_4/k_2)','branch point','(c) Iteration 3: branch point iff k_4 > k_2',C.aqua;
     4,X3,Sv.has_BP,'log_{10}(k_4/k_2)','branch point','(d) Iteration 4: branch point iff k_4 > k_2',C.aqua};
fig=newfig(12,4.0);
for q=1:4
    ax=axes(fig,'Position',[0.06+0.245*(q-1) 0.3 0.19 0.52]); hold(ax,'on');
    m=itv==pan{q,1}; xv=pan{q,2}(m); fl=pan{q,3}(m);
    bw=(max(xv)-min(xv))/44; e=(floor(min(xv)/bw):ceil(max(xv)/bw))*bw;     % bin edges aligned to zero
    n0=histcounts(xv(~fl),e); n1=histcounts(xv(fl),e); ctr=(e(1:end-1)+e(2:end))/2;
    hb=bar(ax,ctr,[n0;n1]','stacked','BarWidth',1,'EdgeColor','w','LineWidth',0.3);
    hb(1).FaceColor=C.grey; hb(2).FaceColor=pan{q,7};
    xline(ax,0,'Color',C.ink,'LineWidth',1.2);
    xlabel(ax,pan{q,4}); ptitle(ax,pan{q,6}); ax.YGrid='on'; ax.GridColor=C.grid; ax.GridAlpha=1;
    if q==1, ylabel(ax,'parameter sets'); end
    legend(ax,hb,{sprintf('no %s (%s)',pan{q,5},commas(sum(~fl))),sprintf('%s (%s)',pan{q,5},commas(sum(fl)))}, ...
        'Box','off','FontSize',6.8,'Location','southoutside');
    wrong=sum(fl & xv<=0)+sum(~fl & xv>0);
    logf('Fig ML.4(%c): iteration %d, %s %d / none %d; sets on the wrong side of zero: %d\n','a'+q-1,pan{q,1},pan{q,5},sum(fl),sum(~fl),wrong);
end
logf('\n');
save_fig(fig,OUT,'ML_fig4_analytic_conditions');

%% =====================================================================================
%% Fig ML.6  ML_fig5_confusion  -- pooled confusion matrix (first split of the ML run)
%% =====================================================================================
CM=[850 0 0 0 0 0 123 37; 0 549 0 0 0 10 156 0; 0 0 1491 0 0 0 0 0; 0 0 0 760 9 0 0 0;
    0 0 0 9 723 0 0 0; 0 25 0 0 0 321 19 0; 136 177 0 0 0 26 2727 8; 81 0 0 0 0 0 12 3742];
try   % prefer the matrix stored in the ML figure, if it can be read
    hf=openfig(fullfile(MLDIR,'figures','confusion_pooled.fig'),'invisible');
    cc=findall(hf,'Type','ConfusionMatrixChart');
    if ~isempty(cc) && isequal(size(cc(1).NormalizedValues),[8 8]), CM=double(cc(1).NormalizedValues); end
    close(hf);
catch
end
Rn=CM./sum(CM,2);
short={'1 LP only (S)','2 LP only (Z)','3 BP only','4 LP+BP (>)','5 LP+BP (<)','6 biphasic','7 LP detached','8 monotone'};
fig=newfig(7.4,6.4); ax=axes(fig,'Position',[0.3 0.24 0.66 0.72]);
imagesc(ax,Rn,[0 1]); colormap(ax,lincmap({'#fcfcfb','#cde2fb','#86b6ef','#2a78d6','#104281'}));
axis(ax,'image'); hold(ax,'on');
for i=1:8, for j=1:8
    if CM(i,j)>0
        tc=C.ink; if Rn(i,j)>0.55, tc='w'; end
        text(ax,j,i,sprintf('%.1f%%\n(%s)',100*Rn(i,j),commas(CM(i,j))),'HorizontalAlignment','center','FontSize',5.8,'Color',tc);
    end
end, end
for k=0.5:1:8.5, plot(ax,[0.5 8.5],[k k],'w-','LineWidth',1.5); plot(ax,[k k],[0.5 8.5],'w-','LineWidth',1.5); end
set(ax,'XTick',1:8,'XTickLabel',short,'XTickLabelRotation',40,'YTick',1:8, ...
    'YTickLabel',strcat(short,{'  (n='},arrayfun(@commas,sum(CM,2)','uni',0),{')'}),'TickLength',[0 0]);
xlabel(ax,'predicted class'); ylabel(ax,'true class');
save_fig(fig,OUT,'ML_fig5_confusion');
logf('Fig ML.6: n = %d test sets, accuracy %.2f%%; class 1->2: %d, class 2->1: %d\n\n',sum(CM(:)),100*trace(CM)/sum(CM(:)),CM(1,2),CM(2,1));

%% ML figures drawn by the ML script itself: copied unchanged
copyfile(fullfile(MLDIR,'figures','importance_by_iteration.png'),fullfile(OUT,'ML_importance_by_iteration.png'));
copyfile(fullfile(MLDIR,'figures','mechanistic_scatter_Iter4_Class4_vs_5.png'),fullfile(OUT,'ML_mechanistic_scatter_Iter4_Class4_vs_5.png'));
copyfile(fullfile(MLDIR,'figures','importance_by_iteration.fig'),fullfile(OUT,'ML_importance_by_iteration.fig'));
copyfile(fullfile(MLDIR,'figures','mechanistic_scatter_Iter4_Class4_vs_5.fig'),fullfile(OUT,'ML_mechanistic_scatter_Iter4_Class4_vs_5.fig'));

%% =====================================================================================
%% Fig MAP.1  MAP_fig1_class_map  -- classes in mechanistic coordinates (v2)
%% =====================================================================================
xm=log10(Pv(:,14).*Pv(:,2)./(Pv(:,5).*Pv(:,1))); ym=log10(Pv(:,8).*Pv(:,1)./(Pv(:,11).*Pv(:,2)));
rng(0,'twister'); bgi=randperm(numel(xm),min(12000,numel(xm)));
fig=newfig(12.4,6.4);
for c=1:8
    r=floor((c-1)/4); q=mod(c-1,4);
    ax=axes(fig,'Position',[0.06+0.235*q 0.56-0.47*r 0.2 0.33]); hold(ax,'on');
    scatter(ax,xm(bgi),ym(bgi),1.2,C.lgrey,'filled');
    m=labv==c; scatter(ax,xm(m),ym(m),2.2,C.blue,'filled','MarkerFaceAlpha',0.55);
    xline(ax,0,'Color',C.ink,'LineWidth',0.7); yline(ax,0,'Color',C.ink,'LineWidth',0.7);
    its=unique(itv(m))'; set(ax,'XLim',[-4.2 5.8],'YLim',[-4.6 5.8]);
    ptitle(ax,sprintf('%d  %s\nn = %s (Iteration %s)',c,CLS{c},commas(sum(m)),strjoin(arrayfun(@num2str,its,'uni',0),', ')));
    if r==1, xlabel(ax,{'log_{10}(k_8P_{TOT}/k_2K_{TOT})','site 1: phosphatase vs kinase'}); end
    if q==0, ylabel(ax,{'log_{10}(k_4K_{TOT}/k_6P_{TOT})','site 2: kinase vs phosphatase'}); end
    logf('Fig MAP.1: class %d: n = %d, iterations %s\n',c,sum(m),mat2str(its));
end
logf('\n'); save_fig(fig,OUT,'MAP_fig1_class_map');

%% =====================================================================================
%% Fig MAP.2  MAP_fig2_embedding_audit  -- embedding vs supervised accuracy (ML outputs)
%% =====================================================================================
Te=readtable(fullfile(MLDIR,'embedding_audit.xlsx'));
g=strcmp(Te.Scope,'global') & strcmp(Te.Representation,'raw log-params');
supA=Tp.Acc_mean(startsWith(Tp.FeatureSet,'A:'));
baseG=[Te.Baseline_Iteration(g) Te.Baseline_Class(g)]; embG=[Te.Embedding_kNN_Iteration(g) Te.Embedding_kNN_Class(g)];
supG=[100 supA];   % LightGBM predicts the sampling iteration with 100% accuracy (ML summary)
fig=newfig(11,3.6);
ax=axes(fig,'Position',[0.07 0.3 0.28 0.6]); hold(ax,'on');
hb=bar(ax,1:2,[baseG;embG;supG]','grouped','EdgeColor','none'); hb(1).FaceColor=C.lgrey; hb(2).FaceColor=C.lblue; hb(3).FaceColor=C.blue;
for j=1:3, xe=hb(j).XEndPoints; v=hb(j).YData; for i=1:2, text(ax,xe(i),v(i)+1.5,sprintf('%.0f',v(i)),'HorizontalAlignment','center','FontSize',7,'Color',C.ink2,'VerticalAlignment','bottom'); end, end
set(ax,'XTick',1:2,'XTickLabel',{'sampling\newlineiteration','bifurcation\newlineclass'},'XTickLabelRotation',0,'YLim',[0 118],'YGrid','on','GridColor',C.grid,'GridAlpha',1);
ylabel(ax,'accuracy (%)'); ptitle(ax,'(a) what a global embedding recovers');
legend(ax,{'majority baseline','kNN in 2-D embedding','supervised (LightGBM, set A)'},'Box','off','FontSize',7,'Location','southoutside');
reg={'KTOT>PTOT','KTOT<PTOT','KTOT=PTOT'}; regl={'K_{TOT}>P_{TOT}\newline(Iteration 1)','K_{TOT}<P_{TOT}\newline(Iteration 2)','K_{TOT}=P_{TOT}\newline(Iterations 3-4)'};
B=zeros(1,3); E=zeros(3,3); Sp=zeros(3,3);
for i=1:3
    r=find(strcmp(Te.Scope,reg{i})); B(i)=Te.Baseline_Class(r(1)); E(i,:)=Te.Embedding_kNN_Class(r)'; Sp(i,:)=Te.Supervised_Class(r)';
end
ax=axes(fig,'Position',[0.45 0.3 0.52 0.6]); hold(ax,'on');
hb=bar(ax,1:3,[B' mean(E,2) mean(Sp,2)],'grouped','EdgeColor','none'); hb(1).FaceColor=C.lgrey; hb(2).FaceColor=C.lblue; hb(3).FaceColor=C.blue;
xe=hb(2).XEndPoints; errorbar(ax,xe(:),mean(E,2),mean(E,2)-min(E,[],2),max(E,[],2)-mean(E,2),'LineStyle','none','Color',C.ink2,'CapSize',3);
xe3=hb(3).XEndPoints; errorbar(ax,xe3(:),mean(Sp,2),mean(Sp,2)-min(Sp,[],2),max(Sp,[],2)-mean(Sp,2),'LineStyle','none','Color',C.ink2,'CapSize',3);
vals=[B' mean(E,2) mean(Sp,2)];
for j=1:3, xe=hb(j).XEndPoints; for i=1:3, text(ax,xe(i),vals(i,j)+3,sprintf('%.0f',vals(i,j)),'HorizontalAlignment','center','FontSize',7,'Color',C.ink2,'VerticalAlignment','bottom'); end, end
set(ax,'XTick',1:3,'XTickLabel',regl,'YLim',[0 112],'YGrid','on','GridColor',C.grid,'GridAlpha',1);
ylabel(ax,'class accuracy (%)'); ptitle(ax,'(b) re-embedding each regime separately (error bars: range over 3 representations)');
legend(ax,{'majority baseline','kNN in 2-D embedding','supervised (LightGBM)'},'Box','off','FontSize',7,'Location','southoutside','NumColumns',3);
save_fig(fig,OUT,'MAP_fig2_embedding_audit');
logf('Fig MAP.2: global embedding: iteration %.1f%% (baseline %.1f%%), class %.1f%% (baseline %.1f%%), supervised class %.2f%%\n\n',embG(1),baseG(1),embG(2),baseG(2),supA);

%% =====================================================================================
%% Fig BP.1  BP_fig1_reproducibility  -- MATCONT without the inverse method vs exact (original ensemble)
%% =====================================================================================
mb=~isnan(O.oldBP) & ~isnan(So.ATOT_BP); bx=So.ATOT_BP(mb); by=O.oldBP(mb); brel=abs(by./bx-1);
ml=ismember(oldc,[4 5]) & ito==4 & ~isnan(O.ATOT_LPBP_LP) & ~isnan(So.LP_curve_AT);
lx=So.LP_curve_AT(ml); ly=O.ATOT_LPBP_LP(ml); lrel=abs(ly./lx-1);
fig=newfig(9.6,4.2);
ax=axes(fig,'Position',[0.1 0.15 0.36 0.72]); hold(ax,'on');
scatter(ax,lx,ly,7,C.orange,'filled','MarkerFaceAlpha',0.55);
scatter(ax,bx,by,7,C.aqua,'filled','MarkerFaceAlpha',0.55);
lim=[0.8*min([bx;lx]) 1.25*max([bx;lx])]; plot(ax,lim,lim,'--','Color',C.ink,'LineWidth',0.7);
tk=10.^(ceil(log10(lim(1))):floor(log10(lim(2))));
set(ax,'XScale','log','YScale','log','XLim',lim,'YLim',lim,'XTick',tk,'YTick',tk);
xlabel(ax,'A_{TOT}, exact (closed form)'); ylabel(ax,{'A_{TOT}, MATCONT','without inverse method'}); ptitle(ax,'(a) location of each singularity');
legend(ax,{sprintf('limit points, Classes 7-8 (n=%s)',commas(numel(lx))),sprintf('branch points (n=%s)',commas(numel(bx))),'identity (MATCONT = exact)'},'Box','off','Location','northwest','FontSize',7);
ax=axes(fig,'Position',[0.57 0.15 0.4 0.72]); hold(ax,'on');
e=logspace(-12,-1,45);
histogram(ax,max(brel,1e-12),e,'FaceColor',C.aqua,'EdgeColor','none','FaceAlpha',0.85);
histogram(ax,max(lrel,1e-12),e,'FaceColor',C.orange,'EdgeColor','none','FaceAlpha',0.75);
set(ax,'XScale','log','YGrid','on','GridColor',C.grid,'GridAlpha',1);
xlabel(ax,'relative difference, MATCONT vs exact'); ylabel(ax,'singularities'); ptitle(ax,'(b) relative difference');
legend(ax,{'branch points','limit points'},'Box','off','Location','northwest','FontSize',7);
text(ax,0.03,0.70,sprintf('median: branch points %s\nmedian: limit points %s\nbranch points < 10^{-5}: %.1f%%\nlimit points < 10^{-4}: %.1f%%', ...
    sci(median(brel)),sci(median(lrel)),100*mean(brel<1e-5),100*mean(lrel<1e-4)),'Units','normalized','FontSize',6.5,'Color',C.ink2,'VerticalAlignment','top');
save_fig(fig,OUT,'BP_fig1_reproducibility');
logf('Fig BP.1: %d branch points (median rel. diff %.2e, %.2f%% < 1e-5); %d limit points (median %.2e, %.2f%% < 1e-4)\n\n', ...
    numel(bx),median(brel),100*mean(brel<1e-5),numel(lx),median(lrel),100*mean(lrel<1e-4));

%% =====================================================================================
%% Fig BP.2  BP_fig2_bp_mechanism  -- pitchfork, Class 4 and Class 5 examples (v2)
%% =====================================================================================
sbp=Sv.s_BP;
pick={find(itv==3 & labv==3 & Sv.ATOT_BP>=40 & Sv.ATOT_BP<=300,1), ...
      find(itv==4 & labv==4 & Sv.ATOT_BP>=40 & Sv.ATOT_BP<=300 & sbp<-0.25 & sbp>-0.6,1), ...
      find(itv==4 & labv==5 & Sv.ATOT_BP>=40 & Sv.ATOT_BP<=300 & sbp>0.25 & sbp<0.6,1)};
fb={find(itv==3 & labv==3,1), find(itv==4 & labv==4,1), find(itv==4 & labv==5,1)};   % fall-backs
for q=1:3, if isempty(pick{q}), pick{q}=fb{q}; end, end
ttl={{'\bf(A)','\rmIteration 3, full symmetry:','pitchfork (s = 0)'}, ...
     {'\bf(B)','\rmIteration 4, Class 7 (s < 0):','fold above the branch point'}, ...
     {'\bf(C)','\rmIteration 4, Class 8 (s > 0):','fold below the branch point'}};
fig=newfig(13,4.8);
for q=1:3
    k=pick{q}; p=Pv(k,:); o=ex_analyse(p); c=ex_consts(p); BP=o.ATOT_BP; xmax=2*BP;
    ax=axes(fig,'Position',[0.06+0.325*(q-1) 0.22 0.26 0.6]); hold(ax,'on');
    L=o.br.line; m=L.AT<xmax; draw_stab(ax,p,L.AT(m),L.x2(m),L.x(:,m),C.blue,2);
    Cb=o.br.curve; m=Cb.AT<xmax;
    for arm={Cb.u>1, Cb.u<1}
        mm=m & arm{1}; draw_stab(ax,p,Cb.AT(mm),Cb.x2(mm),Cb.x(:,mm),C.orange,2);
    end
    plot(ax,BP,o.x2_BP,'s','MarkerFaceColor',C.ink,'MarkerEdgeColor','w','MarkerSize',7);
    text(ax,BP*1.03,o.x2_BP*0.62,'BP','FontSize',8,'VerticalAlignment','top');
    f=Cb.fidx(Cb.AT(Cb.fidx)<xmax);
    if ~isempty(f)
        plot(ax,Cb.AT(f),Cb.x2(f),'o','MarkerFaceColor',C.red,'MarkerEdgeColor','w','MarkerSize',7);
        text(ax,Cb.AT(f(1))*0.9,Cb.x2(f(1)),'LP ','Color',C.red,'HorizontalAlignment','right');
    end
    yy=Cb.x2(m); set(ax,'YScale','log','XLim',[0 xmax],'YLim',[min(1e-1,0.3*min(yy)) 20*max(yy)],'YGrid','on','XGrid','on','GridColor',C.grid,'GridAlpha',1);
    up=m&Cb.u>1; [~,ia]=max(Cb.AT(up)); xu=Cb.AT(up); yu=Cb.x2(up);
    text(ax,xu(ia)*0.98,yu(ia)*1.5,'u > 1 (kinase-rich)','HorizontalAlignment','right','Color',C.orange,'FontSize',7);
    lo=m&Cb.u<1; [~,ib]=max(Cb.AT(lo)); xl=Cb.AT(lo); yl=Cb.x2(lo);
    text(ax,xl(ib)*0.98,yl(ib)*0.45,'u < 1 (phosphatase-rich)','HorizontalAlignment','right','VerticalAlignment','top','Color',C.orange,'FontSize',7);
    xlabel(ax,'A_{TOT}'); ylabel(ax,{'x_2','(doubly phosphorylated)'});
    title(ax,ttl{q},'FontWeight','normal','FontSize',8.5,'Color',[11 11 11]/255); % centred, bold (A)-(C)
    logf('Fig BP.2(%c): v2 iteration %d row %d: ATOT_BP = %.1f, s = %.3f, folds at ATOT = %s\n','a'+q-1,itv(k),k-find(itv==itv(k),1)+1,BP,o.s_BP,mat2str(Cb.AT(f),4));
end
ax=axes(fig,'Position',[0 0 1 0.06],'Visible','off'); hold(ax,'on');
h=[plot(ax,nan,nan,'-','Color',C.blue,'LineWidth',2), plot(ax,nan,nan,'-','Color',C.orange,'LineWidth',2), ...
   plot(ax,nan,nan,'-','Color',C.ink2,'LineWidth',2), plot(ax,nan,nan,'--','Color',C.ink2,'LineWidth',2), ...
   plot(ax,nan,nan,'s','MarkerFaceColor',C.ink,'MarkerEdgeColor',C.ink), plot(ax,nan,nan,'o','MarkerFaceColor',C.red,'MarkerEdgeColor',C.red)];
legend(ax,h,{'u = 1 branch (K = P)','crossing branch A = 1/(\alpha u)','stable','unstable','branch point (closed form)','fold (LP)'}, ...
    'Orientation','horizontal','Box','off','Location','north','FontSize',7.5);
logf('\n'); save_fig(fig,OUT,'BP_fig2_bp_mechanism');

%% =====================================================================================
%% Fig BP.3  BP_fig4_normal_forms  -- normal-form coefficients (v2)
%% =====================================================================================
bpT=NF(strcmp(NF.kind,'BP'),:); lpT=NF(strcmp(NF.kind,'LP'),:);
fig=newfig(13,4.6);
ax=axes(fig,'Position',[0.06 0.24 0.26 0.6]); hold(ax,'on'); e=logspace(-17,0,52);
b3=bpT(bpT.iter==3,:); b4=bpT(bpT.iter==4,:);
histogram(ax,max(abs(b3.a_rel),1e-17),e,'FaceColor',C.blue,'EdgeColor','none','FaceAlpha',0.9);
histogram(ax,max(abs(b4.a_rel),1e-17),e,'FaceColor',C.orange,'EdgeColor','none','FaceAlpha',0.9);
set(ax,'XScale','log','YGrid','on','GridColor',C.grid,'GridAlpha',1); ax.YLim=[0 1.45*ax.YLim(2)];
xlabel(ax,{'|a|, quadratic coefficient of the branching equation','(relative; values < 10^{-17} shown at 10^{-17})'}); ylabel(ax,'branch points');
legend(ax,{sprintf('Iteration 3 (n=%s): pitchfork',commas(height(b3))),sprintf('Iteration 4 (n=%s): transcritical',commas(height(b4)))},'Box','off','Location','northeast','FontSize',7);
ptitle(ax,'(a) branch points: a = 0 only under full symmetry');
ax=axes(fig,'Position',[0.40 0.24 0.25 0.6]); hold(ax,'on');
sl=@(x) sign(x).*log10(1+abs(x)/1e-5);
lb4=labv(b4.row);
for cc={4,C.yellow,'Class 4';5,C.blue,'Class 5'}'
    m=lb4==cc{1}; scatter(ax,sl(b4.s(m)),sl(b4.s_pred(m)),7,cc{2},'filled','MarkerFaceAlpha',0.7);
end
m=~ismember(lb4,[4 5]); scatter(ax,sl(b4.s(m)),sl(b4.s_pred(m)),7,C.grey,'filled','MarkerFaceAlpha',0.7);
lim=sl([-3 3]); plot(ax,lim,lim,'--','Color',C.ink,'LineWidth',0.6); xline(ax,0,'Color',[0.77 0.77 0.72]); yline(ax,0,'Color',[0.77 0.77 0.72]);
tk=[-1 -1e-2 -1e-4 0 1e-4 1e-2 1]; tl={'-10^0','-10^{-2}','-10^{-4}','0','10^{-4}','10^{-2}','10^0'};
set(ax,'XLim',lim,'YLim',lim,'XTick',sl(tk),'XTickLabel',tl,'YTick',sl(tk),'YTickLabel',tl);
xlabel(ax,'s from the exact crossing branch (symmetric log scale)'); ylabel(ax,{'s predicted by','the branching equation'});
legend(ax,{sprintf('Class 4 (n=%s)',commas(sum(lb4==4))),sprintf('Class 5 (n=%s)',commas(sum(lb4==5))),sprintf('other classes (n=%s)',commas(sum(m)))},'Box','off','Location','northwest','FontSize',7);
sgnok=mean(sign(b4.s_pred)==sign(b4.s));
text(ax,0.97,0.04,sprintf('sign agrees in\n%.1f%% of %s sets',100*sgnok,commas(height(b4))),'Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',7,'Color',C.ink2);
ptitle(ax,'(b) Iteration 4: normal form predicts s and its sign');
ax=axes(fig,'Position',[0.74 0.24 0.24 0.6]); hold(ax,'on'); e=logspace(-12,0,49);
brs={'main','main branch, Iter. 1-2',C.blue;'det','detached branch, Iter. 1-2',C.violet;'curve','crossing branch, Iter. 4',C.orange};
for i=1:3
    t=lpT(strcmp(lpT.branch,brs{i,1}),:);
    histogram(ax,abs(t.a_rel),e,'DisplayStyle','stairs','EdgeColor',brs{i,3},'LineWidth',1.6,'DisplayName',sprintf('%s (%s)',brs{i,2},commas(height(t))));
end
set(ax,'XScale','log','YScale','log','YLim',[0.8 3e5],'YGrid','on','GridColor',C.grid,'GridAlpha',1);
xlabel(ax,'|a_{LP}|, fold normal-form coefficient (relative)'); ylabel(ax,'limit points');
lg=legend(ax); set(lg,'Box','off','Location','northwest','FontSize',7); ptitle(ax,sprintf('(c) all %s limit points are non-degenerate',commas(height(lpT))));
save_fig(fig,OUT,'BP_fig4_normal_forms');
relS=abs(b4.s_pred./b4.s-1);
logf(['Fig BP.3 / Table BP.1: branch points: iteration 3 n = %d, max|a| = %.2e (median %.2e); iteration 4 n = %d, |a| %.2e-%.2e (median %.2e)\n' ...
      '  min (b^2-ac)/b^2 = %.3f; max sigma_min/sigma_max of extended Jacobian %.2e, min next %.2e\n' ...
      '  s predicted vs exact: median rel. diff %.1e, sign agreement %.2f%%; Class 4 with s<0: %d of %d; Class 5 with s>0: %d of %d\n'], ...
    height(b3),max(abs(b3.a_rel)),median(abs(b3.a_rel)),height(b4),min(abs(b4.a_rel)),max(abs(b4.a_rel)),median(abs(b4.a_rel)), ...
    min(bpT.disc),max(bpT.sig_rank),min(bpT.sig_next),median(relS),100*sgnok,sum(lb4==4 & b4.s<0),sum(lb4==4),sum(lb4==5 & b4.s>0),sum(lb4==5));
for i=1:3
    t=lpT(strcmp(lpT.branch,brs{i,1}),:);
    logf('  folds, %s: n = %d, |a_LP| %.1e-%.2f (median %.1e), min|beta| %.1e, max sigma_min/sigma_max %.1e\n', ...
        brs{i,2},height(t),min(abs(t.a_rel)),max(abs(t.a_rel)),median(abs(t.a_rel)),min(abs(t.trans)),max(t.sing));
end
logf('  folds whose turning direction matches the computed branch: %d of %d\n\n',sum(lpT.opens==lpT.geom),height(lpT));

%% =====================================================================================
%% Fig EX.2  EX_fig3_detached_criterion  -- closed-form detached-branch criterion (v2, Iter 1-2)
%% =====================================================================================
m12=itv<=2; X=xm(m12); Y=ym(m12); hd=Sv.has_det(m12);   % detached branch at any ATOT (exact analysis)
fig=newfig(4.6,4.2); ax=axes(fig,'Position',[0.16 0.13 0.8 0.78]); hold(ax,'on');
scatter(ax,X(~hd),Y(~hd),2,C.grey,'filled','MarkerFaceAlpha',0.5);
scatter(ax,X(hd),Y(hd),2,C.red,'filled','MarkerFaceAlpha',0.5);
xline(ax,0,'Color',C.ink); yline(ax,0,'Color',C.ink);
text(ax,0.97,0.97,{'k_8P_{TOT} > k_2K_{TOT}','k_4K_{TOT} > k_6P_{TOT}'},'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top');
xlabel(ax,'log_{10}(k_8P_{TOT}/k_2K_{TOT})'); ylabel(ax,'log_{10}(k_4K_{TOT}/k_6P_{TOT})');
legend(ax,{'no detached branch','detached branch'},'Box','off','Location','southwest','FontSize',7.5);
ptitle(ax,sprintf('Iterations 1-2 (n=%s)',commas(sum(m12))));
save_fig(fig,OUT,'EX_fig3_detached_criterion');
crit=X>0 & Y>0;
logf(['Fig EX.2: %d sets with a detached branch (any ATOT), %d without; criterion agrees for %d of %d.\n' ...
      '  (%d sets have a detached branch whose fold lies beyond the window W; the v2 flag Has_Detached_Branch counts only folds inside W.)\n\n'], ...
    sum(hd),sum(~hd),sum(crit==hd),numel(hd),sum(hd & ~detv(m12)));

%% =====================================================================================
%% Fig EX.1  EX_fig1_inverse_method  -- schematic of the inverse method
%% =====================================================================================
draw_cartoon(OUT,C);

%% =====================================================================================
%% Fig EX.3  EX_fig6_missed_example  -- a bistable set labelled "no bifurcation" (original ensemble)
%% =====================================================================================
k=find(ito==2 & oldc==7 & So.has_det & So.det_fold_AT<100,1);
p=Po(k,:); o=ex_analyse(p); Wx=4*So.det_fold_AT(k);
fig=newfig(8.4,3.8);
ax1=axes(fig,'Position',[0.09 0.16 0.39 0.7]); hold(ax1,'on');
b=o.br.main; m=b.AT<Wx; h1=draw_stab(ax1,p,b.AT(m),b.x2(m),b.x(:,m),C.blue,2);
set(ax1,'YScale','log','XLim',[0 Wx]); xlabel(ax1,'A_{TOT}'); ylabel(ax1,'x_2'); ptitle(ax1,'(a) MATCONT without the inverse method');
legend(ax1,h1,{'main branch (from A_{TOT} = 0.1)'},'Box','off','Location','southeast','FontSize',7);
ax2=axes(fig,'Position',[0.58 0.16 0.39 0.7]); hold(ax2,'on');
h1=draw_stab(ax2,p,b.AT(m),b.x2(m),b.x(:,m),C.blue,2);
d=o.br.det; m=d.AT<Wx; h2=draw_stab(ax2,p,d.AT(m),d.x2(m),d.x(:,m),C.red,2);
f=d.fidx(1); plot(ax2,d.AT(f),d.x2(f),'o','MarkerFaceColor',C.red,'MarkerEdgeColor','w','MarkerSize',7);
set(ax2,'YScale','log','XLim',[0 Wx]); xlabel(ax2,'A_{TOT}'); ptitle(ax2,'(b) exact structure: bistable for A_{TOT} > A^{LP}_{TOT}');
legend(ax2,[h1 h2],{'main branch','detached branch'},'Box','off','Location','southeast','FontSize',7);
yl=[min([ax1.YLim(1) ax2.YLim(1)])*0.5 max([ax1.YLim(2) ax2.YLim(2)])]; set([ax1 ax2],'YLim',yl);
save_fig(fig,OUT,'EX_fig6_missed_example');
logf('Fig EX.3: original Iteration-2 set (row %d of the iteration file), fold of the detached branch at ATOT = %.1f\n\n',k-find(ito==2,1)+1,d.AT(f));

%% =====================================================================================
%% Fig EX.4  EX_fig4_class_comparison  -- MATCONT without vs with the inverse method
%%   All 10,000 parameter sets of every iteration: grey = single-start MATCONT (failed or
%%   negative runs as a separate bar "F"), blue = exact ten-class labels (as in Figs 2-5).
%%   Main-branch hysteresis is merged per orientation (1+2, 3+4), because single-start
%%   continuation cannot tell whether a detached branch is also present.
%% =====================================================================================
% columns: tick label | single-start Class_Number_7 values | exact ten-class labels
CAT = {'1+2', 1,     [1 2];
       '3+4', 2,     [3 4];
       '5',   [],    5;
       '6',   3,     6;
       '7',   4,     7;
       '8',   5,     8;
       '9',   6,     9;
       '10',  7,     10;
       'F',   [8 9], []};
KEY = {'1 - LP only, S-shaped', '2 - S-shaped LP + detached LP', '3 - LP only, Z-shaped', ...
       '4 - Z-shaped LP + detached LP', '5 - LP on detached branch only', '6 - BP only', ...
       '7 - LP + BP, x_2^{LP} > x_2^{BP}', '8 - LP + BP, x_2^{LP} < x_2^{BP}', ...
       '9 - no singularity, biphasic', '10 - no singularity, monotone', ...
       'F - single start failed or negative'};
TT  = {'Iteration 1, S-dominant','Iteration 2, Z-dominant', ...
       'Iteration 3, full symmetry','Iteration 4, binding asymmetry'};
FW = 7.2; FH = 7.6;                                   % figure size, inches (~ print size)
in2x = @(v) v/FW; in2y = @(v) v/FH;
fig = newfig(FW, FH);
axW = 2.95; axH = 2.35; axX = [0.62 4.12]; axY = [4.55 1.25];
ssAll = cell(1,4); exAll = cell(1,4);
for it = 1:4
    % single-start labels of all 10,000 sets (Class_Number_7 = 9: failed or negative)
    if strcmp(SINGLE_START,'v3')
        f = fullfile(V3DIR,sprintf('Ensemble_iter%d_v3_1e6',it),'SingleStart','ML_Ready_Dataset_singlestart.csv');
    else
        f = fullfile(OLDDIR,sprintf('Ensemble_iter%d',it),'Excel_Files','ML_Ready_Dataset.csv');
    end
    o = detectImportOptions(f); o.SelectedVariableNames = {'Class_Number_7'};
    o = setvartype(o,'Class_Number_7','double');
    ssAll{it} = readtable(f,o).Class_Number_7;
    % exact ten-class labels of all 10,000 sets (the labels of Figs 2-5)
    f = fullfile(V3DIR,sprintf('Ensemble_iter%d_v3_1e6',it),'MAT_Files','Complete_Bifurcation_Results_v3.mat');
    Sx = load(f,'exact_class_v3'); exAll{it} = Sx.exact_class_v3(:);
end
for it = 1:4
    ss = ssAll{it}; ex = exAll{it}; N = numel(ex);
    vals = zeros(size(CAT,1),2);
    for i = 1:size(CAT,1)
        vals(i,1) = sum(ismember(ss, CAT{i,2}));
        vals(i,2) = sum(ismember(ex, CAT{i,3}));
    end
    if sum(vals(:,1)) ~= N || sum(vals(:,2)) ~= N || numel(ss) ~= N
        warning('Fig EX.4, iteration %d: counts do not add up to N = %d (%d / %d)', it, N, sum(vals(:,1)), sum(vals(:,2)));
    end
    r = floor((it-1)/2)+1; cI = mod(it-1,2)+1;
    ax = axes(fig,'Units','normalized','Position',[in2x(axX(cI)) in2y(axY(r)) in2x(axW) in2y(axH)]);
    hold(ax,'on'); box(ax,'on');
    hb = bar(ax,1:size(CAT,1),vals,'grouped','BarWidth',0.92,'EdgeColor','k','LineWidth',0.6);
    hb(1).FaceColor = C.grey; hb(2).FaceColor = C.blue;
    top = max(vals(:)); yl = top*1.55;
    set(ax,'XLim',[0.4 size(CAT,1)+0.6],'YLim',[0 yl],'XTick',1:size(CAT,1),'XTickLabel',CAT(:,1), ...
        'FontSize',9,'FontWeight','bold','LineWidth',1.2,'TickDir','in','Layer','bottom', ...
        'YGrid','on','GridColor',C.grid,'GridAlpha',1);
    ax.YAxis.Exponent = 0; ytickformat(ax,'%,d');
    % counts above the bars, rotated so that neighbouring labels cannot overlap
    hT = gobjects(0);
    for j = 1:2
        xe = hb(j).XEndPoints;
        for i = 1:size(CAT,1)
            if ~any(vals(i,:)), continue; end      % no label where both bars are zero
            hT(end+1) = text(ax, xe(i), vals(i,j) + top*0.015, commas(vals(i,j)), 'Rotation', 90, ...
                'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',6.5, ...
                'FontWeight','bold','Color',C.ink); %#ok<AGROW>
        end
    end
    % heading inside the box
    msO = 100*sum(ismember(ss,1:5))/N; msN = 100*sum(ismember(ex,1:8))/N;
    hH = text(ax, 0.02, 0.975, {sprintf('%s  (N = %s)', TT{it}, commas(N)), ...
        sprintf('multistationary: %.1f%% \\rightarrow %.1f%%', msO, msN)}, 'Units','normalized', ...
        'HorizontalAlignment','left','VerticalAlignment','top','FontSize',9,'FontWeight','bold', ...
        'BackgroundColor','w','Margin',1);
    % raise the y-limit until the count labels stay below the heading
    drawnow;
    for k = 1:20
        hH.Units = 'data'; eH = hH.Extent; hH.Units = 'normalized';
        E = vertcat(hT.Extent); topLbl = max(E(:,2) + E(:,4));
        if topLbl < eH(2) - 0.01*yl, break; end
        yl = yl*1.08; ylim(ax,[0 yl]); drawnow;
    end
    xlabel(ax,'Classes','FontSize',10,'FontWeight','bold');
    if cI == 1, ylabel(ax,'Number of parameter sets','FontSize',10,'FontWeight','bold'); end
    % panel label (A)-(D), centred above the panel
    annotation(fig,'textbox',[in2x(axX(cI)) in2y(axY(r)+axH+0.04) in2x(axW) in2y(0.28)], ...
        'String',sprintf('(%c)', char('A'+it-1)),'EdgeColor','none','HorizontalAlignment','center', ...
        'VerticalAlignment','bottom','FontSize',13,'FontWeight','bold','Margin',0);
    if it == 1
        lg = legend(ax, hb, {'MATCONT without the inverse method (single start)', ...
            'with the inverse method (exact labels)'}, 'Orientation','horizontal', ...
            'FontSize',9,'FontWeight','bold','Box','off');
        lg.Units = 'normalized'; lg.Position = [in2x(0.3) in2y(FH-0.30) in2x(FW-0.6) in2y(0.25)];
    end
    logf('Fig EX.4: iteration %d (N=%d): multistationary %.1f%% -> %.1f%%; [single start, exact] per column %s: %s\n', ...
        it, N, msO, msN, strjoin(CAT(:,1)',' '), mat2str(vals'));
end
% class key below the panels (three columns)
kc = {KEY(1:4), KEY(5:8), KEY(9:11)}; kx = [0.30 2.45 4.60];
for j = 1:3
    annotation(fig,'textbox',[in2x(kx(j)) in2y(0.0) in2x(2.55) in2y(0.70)],'String',kc{j}, ...
        'EdgeColor','none','HorizontalAlignment','left','VerticalAlignment','top', ...
        'FontSize',8.5,'FontWeight','bold','Margin',0,'Interpreter','tex');
end
logf('\n'); save_fig(fig,OUT,'EX_fig4_class_comparison',1);

%% =====================================================================================
%% Fig SS.1  SS_fig1_single_start_limits  -- limits of MATCONT without the inverse method
%% =====================================================================================
r7=reach(oldc==7 & ~isnan(reach));
mbp=ito>=3 & oldc==7 & ismember(So.exact,[3 4 5]); bpm=So.ATOT_BP(mbp);
s7=ito<=2 & oldc==7 & So.has_det & ~isnan(reach); rr=reach(s7); df=So.det_fold_AT(s7); ins=df<rr;
fig=newfig(12.5,4.6);
ax=axes(fig,'Position',[0.06 0.32 0.25 0.52]); hold(ax,'on'); e=logspace(1,6.3,60);
histogram(ax,r7,e,'FaceColor',C.grey,'EdgeColor','none'); histogram(ax,bpm,e,'FaceColor',C.aqua,'EdgeColor','none','FaceAlpha',0.9);
set(ax,'XScale','log','YScale','log','YGrid','on','GridColor',C.grid,'GridAlpha',1); xline(ax,W,'--','Color',C.ink);
ax.YLim=[0.8 ax.YLim(2)*20]; text(ax,W/1.4,ax.YLim(2)*0.6,sprintf('window\nW = 10^{%d}',round(log10(W))),'FontSize',7,'Color',C.ink2,'VerticalAlignment','top','HorizontalAlignment','right');
xlabel(ax,'A_{TOT}'); ylabel(ax,'parameter sets'); ptitle(ax,{'(a) where MATCONT without the inverse','     method stopped (no bifurcation found)'});
legend(ax,{sprintf('upper end of continuation (%s sets)',commas(numel(r7))),sprintf('missed branch points, Iter. 3-4 (%s)',commas(numel(bpm)))},'Box','off','Location','southoutside','FontSize',7);
ax=axes(fig,'Position',[0.40 0.32 0.24 0.52]); hold(ax,'on');
scatter(ax,rr(~ins),df(~ins),3,C.grey,'filled','MarkerFaceAlpha',0.5); scatter(ax,rr(ins),df(ins),3,C.red,'filled','MarkerFaceAlpha',0.5);
xr=[floor(min(rr)/100)*100, ceil(max(rr)/100)*100];            % axis spans the reach of the data
xx=linspace(xr(1),xr(2),200); plot(ax,xx,xx,'-','Color',C.ink); text(ax,xr(1)+0.02*diff(xr),1.3*xr(1),'fold = reach','HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',7,'Color',C.ink2);
set(ax,'YScale','log','XLim',xr,'YLim',[10 1e6]);
xlabel(ax,'upper end of continuation (A_{TOT})'); ylabel(ax,{'fold of detached','branch (A_{TOT})'});
ptitle(ax,{'(b) detached folds the continuation passed','     (no bifurcation found, Iter. 1-2)'});
legend(ax,{sprintf('fold beyond reach (%s)',commas(sum(~ins))),sprintf('fold inside reach (%s)',commas(sum(ins)))},'Box','off','Location','southoutside','FontSize',7);
% (c) ode15s from x = 0 at ATOT = 0.1 up to t = 5000 (the original initialisation), compared with the exact stable state
cf=fullfile(OUT,sprintf('cache_convergence_%s.mat',SINGLE_START));
if USE_CACHE && exist(cf,'file') && ismember('fres',who('-file',cf)), load(cf,'rd','fres');
else
    rng(7,'twister'); sel=[]; for it=1:4, ii=find(ito==it); sel=[sel; ii(randperm(numel(ii),NCONV))]; end %#ok<AGROW>
    rd=nan(numel(sel),1); fres=nan(numel(sel),1); opt=odeset('RelTol',1e-8,'AbsTol',1e-12);
    for j=1:numel(sel)
        p=Po(sel(j),:); a=0.1;
        try
            [~,Y]=ode15s(@(t,x) ex_rhs(x,a,p),[0 5000],zeros(6,1),opt);
        catch
            [~,Y]=ode23s(@(t,x) ex_rhs(x,a,p),[0 5000],zeros(6,1),opt);
        end
        xe=Y(end,:)'; fres(j)=norm(ex_rhs(xe,a,p),Inf);
        SS=ex_all_ss(p,a); best=Inf;
        for q=1:numel(SS), if SS(q).stable, best=min(best,norm(xe-SS(q).x)/max(norm(SS(q).x),1e-300)); end, end
        rd(j)=best;
    end
    save(cf,'rd','fres');
end
if ~exist('fres','var'), fres=nan(size(rd)); end
ax=axes(fig,'Position',[0.73 0.32 0.25 0.52]); hold(ax,'on');
histogram(ax,log10(max(rd,1e-16)),linspace(-16,0,49),'FaceColor',C.blue,'EdgeColor','none');
xline(ax,-3,'--','Color',C.ink); nc=sum(rd>1e-3); ax.YLim=[0 ax.YLim(2)*1.15];
text(ax,-3.3,ax.YLim(2)*0.95,sprintf('%d of %d runs\nnot converged ->',nc,numel(rd)),'FontSize',7,'Color',C.ink2,'VerticalAlignment','top','HorizontalAlignment','right');
xlabel(ax,{'log_{10} relative distance to','steady state at t = 5000'}); ylabel(ax,'integrations'); ax.YGrid='on'; ax.GridColor=C.grid; ax.GridAlpha=1;
ptitle(ax,{'(c) ode15s from x = 0 at A_{TOT} = 0.1',sprintf('     (%d random sets)',numel(rd))});
save_fig(fig,OUT,'SS_fig1_single_start_limits');
logf(['Fig SS.1: (a) %d single-start Class-7 sets; median upper end of continuation %.0f; %d missed branch points (median ATOT %.0f)\n' ...
      '  (b) %d Iteration-1/2 Class-7 sets with a detached branch: fold inside reach %d, beyond %d\n' ...
      '  (c) %d of %d ode15s runs not converged (relative distance > 1e-3); %d of these flagged by |dx/dt|_inf > 1e-6\n' ...
      '  reach of all single-start runs (median, 5-95%%): %.0f (%.0f-%.0f); of the no-bifurcation runs: %.0f (%.0f-%.0f)\n\n'], ...
    numel(r7),median(r7),numel(bpm),median(bpm),numel(rr),sum(ins),sum(~ins),nc,numel(rd),sum(rd>1e-3 & fres>1e-6), ...
    median(O.ATOT_reached,'omitnan'),prctile_(O.ATOT_reached,5),prctile_(O.ATOT_reached,95),median(r7),prctile_(r7,5),prctile_(r7,95));

%% =====================================================================================
%% Methods "Verification": checks of the parametrisation on the single-start ensemble
%% =====================================================================================
cf=fullfile(OUT,sprintf('cache_verification_%s.mat',SINGLE_START));
if USE_CACHE && exist(cf,'file'), load(cf,'VER');
else
    rs=RandStream('mt19937ar','Seed',11); VER=struct();
    % (i) scaled residual of the ODEs at every steady state from the parametrisation,
    %     1,500 random parameter sets, random ATOT in [0.1, 1e3] (log-uniform)
    sel=randperm(rs,size(Po,1),1500); res=[];
    for j=sel
        p=Po(j,:); a=10^(-1+4*rand(rs));
        SS=ex_all_ss_br(p,a);
        for q=1:numel(SS)
            [F,v]=ex_rhs_rates(SS(q).x,a,p); res(end+1)=norm(F,Inf)/max(max(abs(v)),realmin); %#ok<AGROW>
        end
    end
    VER.res=res;
    % (ii) 600 random sets with a detached branch: steady states at twice the fold ATOT
    id=find(ito<=2 & So.has_det & isfinite(So.det_fold_AT));
    sel=id(randperm(rs,numel(id),min(600,numel(id)))); nss=zeros(numel(sel),1); nst=nss;
    for j=1:numel(sel)
        p=Po(sel(j),:); a=2*So.det_fold_AT(sel(j)); SS=ex_all_ss_br(p,a);
        nss(j)=numel(SS); nst(j)=sum([SS.stable]);
    end
    VER.nss=nss; VER.nst=nst;
    % (iii) 100 of these: integrate from x = 0 at twice the fold ATOT; which stable state is reached?
    sel=sel(1:min(100,numel(sel))); hitdet=false(numel(sel),1); conv=false(numel(sel),1);
    opt=odeset('RelTol',1e-8,'AbsTol',1e-12);
    for j=1:numel(sel)
        p=Po(sel(j),:); a=2*So.det_fold_AT(sel(j)); SS=ex_all_ss_br(p,a);
        try, [~,Y]=ode15s(@(t,x) ex_rhs(x,a,p),[0 5000],zeros(6,1),opt);
        catch, [~,Y]=ode23s(@(t,x) ex_rhs(x,a,p),[0 5000],zeros(6,1),opt); end
        xe=Y(end,:)'; d=inf(numel(SS),1);
        for q=1:numel(SS), if SS(q).stable, d(q)=norm(xe-SS(q).x)/max(norm(SS(q).x),1e-300); end, end
        [dm,qm]=min(d); conv(j)=dm<1e-3; hitdet(j)=conv(j) && strcmp(SS(qm).br,'det');
    end
    VER.hitdet=hitdet; VER.conv=conv;
    save(cf,'VER');
end
logf(['Methods verification: (i) %d steady states from 1,500 random sets: scaled residual median %.1e, max %.1e\n' ...
      '  (ii) %d sets with a detached branch at twice the fold ATOT: %d with exactly three steady states, %d of them with two stable\n' ...
      '  (iii) integration from x = 0 at twice the fold ATOT: %d of %d reach the stable state of the detached branch (%d converged)\n\n'], ...
    numel(VER.res),median(VER.res),max(VER.res),numel(VER.nss),sum(VER.nss==3),sum(VER.nss==3 & VER.nst==2), ...
    sum(VER.hitdet),numel(VER.hitdet),sum(VER.conv));



%% =====================================================================================
%% Fig SS.2  SS_fig2_biphasic  -- robustness of the biphasic label (v2, Class 6)
%% =====================================================================================
k6=find(labv==6); rel=nan(numel(k6),1); pkA=rel; curves=cell(numel(k6),1);
for j=1:numel(k6)
    o=ex_analyse(Pv(k6(j),:)); b=o.br.main; m=b.AT<=W;
    [ATs,ord]=sort(b.AT(m)); x2=b.x2(m); x2=x2(ord); [pk,i]=max(x2);
    amp=pk-min(x2); prom=pk-max(min(x2(1:i)),min(x2(i:end)));
    rel(j)=prom/amp; pkA(j)=ATs(i); curves{j}={ATs,x2/pk,ATs(i)};
end
fig=newfig(12.5,4.3);
ax=axes(fig,'Position',[0.06 0.18 0.25 0.6]); hold(ax,'on'); rng(3,'twister');
for j=randperm(numel(k6),min(8,numel(k6)))
    cv=curves{j}; plot(ax,cv{1},cv{2},'-','Color',[C.blue 0.85],'LineWidth',1.2); plot(ax,cv{3},1,'o','MarkerFaceColor',C.red,'MarkerEdgeColor',C.red,'MarkerSize',4);
end
set(ax,'XScale','log','XLim',[1e-1 W],'YLim',[0 1.08],'YGrid','on','GridColor',C.grid,'GridAlpha',1);
xlabel(ax,'A_{TOT}'); ylabel(ax,'x_2 / peak value'); ptitle(ax,{'(a) exact main branches of eight biphasic sets','     (red: maximum)'});
ax=axes(fig,'Position',[0.39 0.18 0.25 0.6]); hold(ax,'on');
histogram(ax,max(rel,1e-3)*100,logspace(-1.3,2.01,50),'FaceColor',C.blue,'EdgeColor','none');
set(ax,'XScale','log','YScale','log','XLim',[0.05 130],'YLim',[0.8 2e4]);
taus=[0.1 1 5 10]; flips=arrayfun(@(t) sum(rel<t/100),taus);
for t=taus, xline(ax,t,'--','Color',C.ink); text(ax,t*1.1,1.2,sprintf('\\tau = %g%%',t),'Rotation',90,'VerticalAlignment','top','FontSize',6.5,'Color',C.ink2); end
xlabel(ax,'peak prominence (% of branch amplitude)'); ylabel(ax,'parameter sets');
ptitle(ax,{'(b) prominence of the maximum',sprintf('     (%s sets; flips at \\tau = 0.1/1/5/10%%: %d/%d/%d/%d)',commas(numel(k6)),flips)});
ax=axes(fig,'Position',[0.73 0.18 0.25 0.6]); hold(ax,'on');
histogram(ax,pkA,logspace(-1,log10(W),50),'FaceColor',C.blue,'EdgeColor','none');
set(ax,'XScale','log','YGrid','on','GridColor',C.grid,'GridAlpha',1); xline(ax,W,'--','Color',C.ink);
xlabel(ax,'A_{TOT} at the maximum'); ylabel(ax,'parameter sets'); ptitle(ax,'(c) location of the maximum');
save_fig(fig,OUT,'SS_fig2_biphasic');
logf('Fig SS.2: %d Class-6 sets; label flips at tau = 0.1/1/5/10%%: %s; median prominence %.1f%%; peak ATOT 5/50/95%%: %s\n\n', ...
    numel(k6),mat2str(flips),100*median(rel),mat2str(prctile_(pkA,[5 50 95]),3));

fclose(LOG);
fprintf('\nDone. Figures in %s\n',OUT);
%% =====================================================================================
%%  LOCAL FUNCTIONS  --  exact steady-state structure (u = K/P parametrisation)
%% =====================================================================================
function c = ex_consts(p)
% p = [KTOT PTOT k1 k11 k2 k3 k33 k4 k5 k55 k6 k7 k77 k8]
KT=p(1); PT=p(2); k1=p(3); k11=p(4); k2=p(5); k3=p(6); k33=p(7); k4=p(8);
k5=p(9); k55=p(10); k6=p(11); k7=p(12); k77=p(13); k8=p(14);
c.KT=KT; c.PT=PT;
c.KM1=(k11+k2)/k1; c.KM3=(k33+k4)/k3; c.KM5=(k55+k6)/k5; c.KM7=(k77+k8)/k7;
c.a1=k2*c.KM7/(k8*c.KM1); c.a2=k4*c.KM5/(k6*c.KM3); c.r=KT/PT;
c.al=c.a1/c.KM3*(c.r*k4/k6-1); c.be=1/c.KM1*(c.r*k2/k8-1);
if c.al~=0, c.ustar=-c.be/c.al; else, c.ustar=Inf; end
end

function [AT,x,K,P] = ex_state(c,u,A)
% full steady state for given u and A (row vectors); x is 6 x N
u=u(:)'; A=A(:)';
K=c.KT./(1+A./c.KM1+c.a1.*u.*A./c.KM3); P=K./u;
x1=c.a1.*u.*A; x2=c.a2.*u.*x1;
x3=A.*K./c.KM1; x4=x1.*K./c.KM3; x5=x2.*P./c.KM5; x6=x1.*P./c.KM7;
AT=A+x1+x2+x3+x4+x5+x6;
x=[x1;x2;x3;x4;x5;x6];
end

function J = ex_jac(p,AT,x)
% analytic Jacobian of the six ODEs (A, K, P eliminated by conservation)
KT=p(1); PT=p(2); k1=p(3); k11=p(4); k2=p(5); k3=p(6); k33=p(7); k4=p(8);
k5=p(9); k55=p(10); k6=p(11); k7=p(12); k77=p(13); k8=p(14);
x=x(:); x1=x(1); x2=x(2); x3=x(3); x4=x(4); x5=x(5); x6=x(6);
A=AT-sum(x); K=KT-(x3+x4); P=PT-(x5+x6);
dA=-ones(1,6); dK=[0 0 -1 -1 0 0]; dP=[0 0 0 0 -1 -1]; e=eye(6);
v1=k1*(dA*K+A*dK)-k11*e(3,:); v2=k2*e(3,:);
v3=k3*(e(1,:)*K+x1*dK)-k33*e(4,:); v4=k4*e(4,:);
v5=k5*(e(2,:)*P+x2*dP)-k55*e(5,:); v6=k6*e(5,:);
v7=k7*(e(1,:)*P+x1*dP)-k77*e(6,:); v8=k8*e(6,:);
J=[v2-v3+v6-v7; v4-v5; v1-v2; v3-v4; v5-v6; v7-v8];
end

function F = ex_rhs(x,AT,p)
% right-hand side of the six ODEs at fixed ATOT
KT=p(1); PT=p(2); k1=p(3); k11=p(4); k2=p(5); k3=p(6); k33=p(7); k4=p(8);
k5=p(9); k55=p(10); k6=p(11); k7=p(12); k77=p(13); k8=p(14);
x=x(:); x1=x(1); x2=x(2); x3=x(3); x4=x(4); x5=x(5); x6=x(6);
A=AT-sum(x); K=KT-x3-x4; P=PT-x5-x6;
v1=k1*A*K-k11*x3; v2=k2*x3; v3=k3*x1*K-k33*x4; v4=k4*x4;
v5=k5*x2*P-k55*x5; v6=k6*x5; v7=k7*x1*P-k77*x6; v8=k8*x6;
F=[v2-v3+v6-v7; v4-v5; v1-v2; v3-v4; v5-v6; v7-v8];
end

function idx = ex_folds(AT)
% indices of turning points of ATOT along a branch
d=diff(AT); s=sign(d); s(s==0)=1;
idx=find(s(2:end)~=s(1:end-1))+1;
end

function u = ex_grid(lo,hi,n)
% logarithmic grid on (lo,hi) accumulating at finite end points
if nargin<3, n=3000; end
h=n/2;
if isinf(hi)
    u=lo*(1+logspace(-12,14,n));
elseif lo==0
    u=[hi*logspace(-14,-0.3,h), hi*(1-logspace(-0.3,-12,h))];
else
    t=[logspace(-12,-0.31,h), 1-logspace(-0.31,-12,h)];
    if hi/lo>10, u=exp(log(lo)+(log(hi)-log(lo))*t); else, u=lo+(hi-lo)*t; end
end
end

function o = ex_analyse(p)
% complete steady-state structure of one parameter set (all ATOT)
c=ex_consts(p); r=c.r; al=c.al; be=c.be; us=c.ustar;
o.r=r; o.alpha=al; o.beta=be; o.ustar=us;
o.sym = abs(us-r) <= 1e-9*r;
o.has_BP=false; o.ATOT_BP=NaN; o.x2_BP=NaN; o.s_BP=NaN;
br=struct();
if o.sym
    Aline=logspace(-10,8,4000);
    br.line.u=r*ones(size(Aline)); br.line.A=Aline;
    if al>0
        u=logspace(-12,14,6000); br.curve.u=u; br.curve.A=1./(al*u);
        o.has_BP=true;
        [ATbp,xb]=ex_state(c,r,1/al); o.ATOT_BP=ATbp; o.x2_BP=xb(2);
        g=@(v) ex_state(c,v,1./(al*v)); hh=1e-5;
        o.s_BP=(g(r*exp(hh))-g(r*exp(-hh)))/(2*hh)/ATbp;
    end
else
    det=[];
    if al>0 && be<0
        if r<us, main=[0 r]; det=[us Inf]; else, main=[r Inf]; det=[0 us]; end
    elseif (al>0 && be>0) || (al==0 && be>0)
        main=[r Inf];
    elseif (al<0 && be<0) || (al==0 && be<0)
        main=[0 r];
    else
        main=[min(r,us) max(r,us)];
    end
    u=ex_grid(main(1),main(2)); A=(u-r)./(u.*(al*u+be)); ok=A>0;
    br.main.u=u(ok); br.main.A=A(ok);
    if ~isempty(det)
        u=ex_grid(det(1),det(2)); A=(u-r)./(u.*(al*u+be)); ok=A>0;
        br.det.u=u(ok); br.det.A=A(ok);
    end
end
names=fieldnames(br);
for i=1:numel(names)
    b=br.(names{i});
    [AT,x]=ex_state(c,b.u,b.A);
    if strcmp(names{i},'line'), [~,ord]=sort(b.A); else, [~,ord]=sort(b.u); end
    AT=AT(ord); x=x(:,ord); uu=b.u(ord); AA=b.A(ord);
    f=ex_folds(AT);
    if strcmp(names{i},'curve') && abs(o.s_BP)<1e-8
        % under exact symmetry the turning point AT the branch point is the pitchfork, not a fold
        f=f(~((abs(log(uu(f)/r))<0.05) & (abs(AT(f)/o.ATOT_BP-1)<1e-3)));
    end
    b.u=uu; b.A=AA; b.AT=AT; b.x=x; b.x2=x(2,:); b.fidx=f; b.nf=numel(f); b.fold_AT=AT(f);
    br.(names{i})=b;
end
o.br=br;
o.has_det = isfield(br,'det') && ~isempty(br.det.AT);
o.nLP_main=0; o.nLP_curve=0; o.det_fold_AT=NaN; o.main_fold_AT=[];
if isfield(br,'main'), o.nLP_main=br.main.nf; o.main_fold_AT=br.main.fold_AT; end
if isfield(br,'curve'), o.nLP_curve=br.curve.nf; end
if o.has_det && br.det.nf>0, o.det_fold_AT=min(br.det.fold_AT); end
o.biphasic=false; o.peak_AT=NaN;
if isfield(br,'main') && br.main.nf==0 && numel(br.main.AT)>2
    [ATs,ord]=sort(br.main.AT); x2=br.main.x2(ord); [~,i]=max(x2);
    o.biphasic = i>1 && i<numel(x2) && x2(i)>x2(end)*(1+1e-6);
    o.peak_AT=ATs(i);
end
% exact label over all ATOT
o.LP_curve_AT=NaN;
if o.sym
    if ~o.has_BP, o.exact=8;
    elseif o.nLP_curve==0, o.exact=3;
    else
        f=br.curve.fidx; o.LP_curve_AT=min(br.curve.AT(f));
        if min(br.curve.x2(f))>o.x2_BP, o.exact=4; else, o.exact=5; end
    end
else
    if o.nLP_main>0, if p(1)>p(2), o.exact=1; else, o.exact=2; end
    elseif o.has_det, o.exact=7;
    elseif o.biphasic, o.exact=6;
    else, o.exact=8;
    end
end
end

function w = ex_window_label(R,W)
% class when only ATOT <= W is considered (R = row struct of ex_summary fields)
e=R.exact; w=e;
mfm=R.main_fold_min; if isnan(mfm), mfm=Inf; end
dfa=R.det_fold_AT; if isnan(dfa), dfa=Inf; end
pk=R.peak_AT; if isnan(pk), pk=Inf; end
switch e
    case {1,2}
        if mfm<=W, w=e;
        elseif R.has_det && dfa<=W, w=7;
        elseif R.biphasic && pk<=W, w=6;
        else, w=8; end
    case {4,5}
        if R.ATOT_BP>W, w=8; elseif R.LP_curve_AT<=W, w=e; else, w=3; end
    case 3
        if R.ATOT_BP<=W, w=3; else, w=8; end
    case 7
        if dfa<=W, w=7; elseif R.biphasic && pk<=W, w=6; else, w=8; end
end
end

function S = ex_summary(P, label)
% exact analysis of every row of P (N x 14); returns a struct of column vectors
N=size(P,1);
S.exact=zeros(N,1); S.has_BP=false(N,1); S.has_det=false(N,1); S.biphasic=false(N,1);
S.nLP_main=zeros(N,1); S.nLP_curve=zeros(N,1);
f={'det_fold_AT','ATOT_BP','x2_BP','s_BP','LP_curve_AT','main_fold_min','main_fold_max','peak_AT','alpha','beta','r','ustar'};
for i=1:numel(f), S.(f{i})=nan(N,1); end
t0=tic;
for k=1:N
    o=ex_analyse(P(k,:));
    S.exact(k)=o.exact; S.has_BP(k)=o.has_BP; S.has_det(k)=o.has_det; S.biphasic(k)=o.biphasic;
    S.nLP_main(k)=o.nLP_main; S.nLP_curve(k)=o.nLP_curve;
    S.det_fold_AT(k)=o.det_fold_AT; S.ATOT_BP(k)=o.ATOT_BP; S.x2_BP(k)=o.x2_BP; S.s_BP(k)=o.s_BP;
    S.LP_curve_AT(k)=o.LP_curve_AT; S.peak_AT(k)=o.peak_AT;
    if ~isempty(o.main_fold_AT), S.main_fold_min(k)=min(o.main_fold_AT); S.main_fold_max(k)=max(o.main_fold_AT); end
    S.alpha(k)=o.alpha; S.beta(k)=o.beta; S.r(k)=o.r; S.ustar(k)=o.ustar;
    if mod(k,5000)==0, fprintf('  [%s] %d / %d  (%.0f s)\n',label,k,N,toc(t0)); end
end
end

function w = ex_window_all(S,W)
N=numel(S.exact); w=zeros(N,1);
fn=fieldnames(S);
for k=1:N
    R=struct(); for i=1:numel(fn), v=S.(fn{i}); R.(fn{i})=v(k); end
    w(k)=ex_window_label(R,W);
end
end
%% =====================================================================================
%%  LOCAL FUNCTIONS  --  normal-form coefficients at branch points and folds
%% =====================================================================================
function F = nf_F(y,p)
% steady-state equations F(y), y = (x1..x6, ATOT); F is quadratic in y
F=ex_rhs(y(1:6),y(7),p);
end

function [S,R,Js,Bf] = nf_setup(p,y)
% dimensionless coordinates: variables scaled by their value, equations by the largest flux
y=y(:); S=abs(y); S(S==0)=1;
KT=p(1); PT=p(2); k1=p(3); k11=p(4); k2=p(5); k3=p(6); k33=p(7); k4=p(8);
k5=p(9); k55=p(10); k6=p(11); k7=p(12); k77=p(13); k8=p(14);
x1=y(1); x2=y(2); x3=y(3); x4=y(4); x5=y(5); x6=y(6); AT=y(7);
A=AT-sum(y(1:6)); K=KT-x3-x4; P=PT-x5-x6;
v=abs([k1*A*K,k11*x3,k2*x3,k3*x1*K,k33*x4,k4*x4,k5*x2*P,k55*x5,k6*x5,k7*x1*P,k77*x6,k8*x6]);
R=[max(v([3 4 5 9 10 11])); max(v([6 7 8])); max(v([1 2 3])); max(v([4 5 6])); max(v([7 8 9])); max(v([10 11 12]))];
J=ex_jac(p,AT,y(1:6)); Fa=zeros(6,1); Fa(3)=k1*K;
Jfull=[J Fa];
Js=(Jfull.*repmat(S',6,1))./repmat(R,1,7);
z0=zeros(7,1);
Bf=@(a,b) (nf_F((a+b).*S,p)-nf_F(a.*S,p)-nf_F(b.*S,p)+nf_F(z0,p))./R;   % exact for quadratic F
end

function v = nf_dlnu(K,P)
v=zeros(6,1); v(3)=-1/K; v(4)=-1/K; v(5)=1/P; v(6)=1/P;   % d ln(K/P)
end

function out = nf_bp(p,y,dlnu)
[S,~,Js,B]=nf_setup(p,y); y=y(:);
sv=svd(Js);
Jx=Js(:,1:6); [Ux,sx,Vx]=svd(Jx); sx=diag(sx);
phi=[Vx(:,end);0]; psi=Ux(:,end);
if dot(phi(1:6).*S(1:6),dlnu)<0, phi=-phi; end
w=pinv(Jx)*(-Js(:,7)); w=w-dot(w,phi(1:6))*phi(1:6); wv=[w;1];
a=psi'*B(phi,phi); b=psi'*B(phi,wv); c=psi'*B(wv,wv);
if b<0, a=-a; b=-b; c=-c; end
nrm=norm(B(phi,phi));
if abs(c)>0, taus=real(roots([c 2*b a])); else, taus=-a/(2*b); end
best=-Inf; s_pred=NaN;
for t=taus(:)'
    tv=phi+t*wv; dy=tv.*S; dAT=dy(7); dl=dot(dy(1:6),dlnu);
    if abs(dl)>best, best=abs(dl); s_pred=dAT/y(7)/dl; end
end
out=struct('a',a,'b',b,'c',c,'a_rel',a/nrm,'disc',(b*b-a*c)/(b*b), ...
    'sig_rank',sv(end)/sv(1),'sig_next',sv(end-1)/sv(1),'Jx_sing',sx(end)/sx(1),'s_pred',s_pred);
end

function out = nf_lp(p,y,dlnu)
[S,~,Js,B]=nf_setup(p,y);
Jx=Js(:,1:6); [Ux,sx,Vx]=svd(Jx); sx=diag(sx);
q=Vx(:,end); pv=Ux(:,end);
if dot(q.*S(1:6),dlnu)<0, q=-q; end
pv=pv/dot(pv,q);                         % MATCONT normalisation p.q = 1
qq=[q;0]; a=0.5*pv'*B(qq,qq); trans=pv'*Js(:,7);
nrm=norm(B(qq,qq))*norm(pv);
out=struct('a_LP',a,'a_rel',a/(0.5*nrm),'trans',trans,'sing',sx(end)/sx(1),'opens',sign(-a/trans));
end

function T = nf_all(P, iters, label)
% normal-form coefficients at every branch point (symmetric iterations) and every fold
rows={}; t0=tic;
for k=1:size(P,1)
    p=P(k,:); c=ex_consts(p); o=ex_analyse(p);
    if o.has_BP
        [AT,x,K,Pf]=ex_state(c,c.r,1/c.al); y=[x;AT];
        q=nf_bp(p,y,nf_dlnu(K,Pf));
        rows(end+1,:)={k,iters(k),'BP','crossing',q.a_rel,q.disc,q.sig_rank,q.sig_next,q.s_pred,o.s_BP,NaN,NaN,NaN,NaN,AT}; %#ok<AGROW>
    end
    nm=fieldnames(o.br);
    for i=1:numel(nm)
        if strcmp(nm{i},'line'), continue; end
        b=o.br.(nm{i});
        if strcmp(nm{i},'curve'), Af=@(v) 1./(c.al*v); else, Af=@(v) (v-c.r)./(v.*(c.al*v+c.be)); end
        g=@(v) ex_state(c,v,Af(v));
        dg=@(v) g(v*(1+1e-7))-g(v*(1-1e-7));
        for f=b.fidx(:)'
            u0=b.u(max(f-1,1)); u1=b.u(min(f+1,numel(b.u)));
            try, uf=fzero(dg,[u0 u1]); catch, uf=b.u(f); end
            [AT,x,K,Pf]=ex_state(c,uf,Af(uf)); y=[x;AT];
            ismin = g(uf*1.001)>AT && g(uf/1.001)>AT;
            ismax = g(uf*1.001)<AT && g(uf/1.001)<AT;
            geom = 1*ismin - 1*ismax;
            q=nf_lp(p,y,nf_dlnu(K,Pf));
            rows(end+1,:)={k,iters(k),'LP',nm{i},q.a_rel,NaN,NaN,NaN,NaN,NaN,q.trans,q.sing,q.opens,geom,AT}; %#ok<AGROW>
        end
    end
    if mod(k,5000)==0, fprintf('  [%s] %d / %d  (%.0f s)\n',label,k,size(P,1),toc(t0)); end
end
T=cell2table(rows,'VariableNames',{'row','iter','kind','branch','a_rel','disc','sig_rank','sig_next', ...
    's_pred','s','trans','sing','opens','geom','ATOT'});
end

function SS = ex_all_ss(p,a)
% all steady states at ATOT = a (from the parametrisation), with stability
o=ex_analyse(p); c=ex_consts(p); SS=struct('x',{},'stable',{});
nm=fieldnames(o.br);
for i=1:numel(nm)
    b=o.br.(nm{i}); d=b.AT-a; idx=find(sign(d(1:end-1))~=sign(d(2:end)));
    for j=idx(:)'
        if strcmp(nm{i},'line')
            g=@(A) ex_state(c,c.r,A)-a; A0=fzero(g,[b.A(j) b.A(j+1)]); u0=c.r;
        else
            if strcmp(nm{i},'curve'), Af=@(v) 1./(c.al*v); else, Af=@(v) (v-c.r)./(v.*(c.al*v+c.be)); end
            g=@(v) ex_state(c,v,Af(v))-a; u0=fzero(g,[b.u(j) b.u(j+1)]); A0=Af(u0);
        end
        [~,x]=ex_state(c,u0,A0);
        ev=eig(ex_jac(p,a,x));
        SS(end+1)=struct('x',x,'stable',all(real(ev)<0)); %#ok<AGROW>
    end
end
end

function SS = ex_all_ss_br(p,a)
% as ex_all_ss, with the name of the branch of every steady state ('main','det','line','curve')
o=ex_analyse(p); c=ex_consts(p); SS=struct('x',{},'stable',{},'br',{});
nm=fieldnames(o.br);
for i=1:numel(nm)
    b=o.br.(nm{i}); d=b.AT-a; idx=find(sign(d(1:end-1))~=sign(d(2:end)));
    for j=idx(:)'
        if strcmp(nm{i},'line')
            g=@(A) ex_state(c,c.r,A)-a; A0=fzero(g,[b.A(j) b.A(j+1)]); u0=c.r;
        else
            if strcmp(nm{i},'curve'), Af=@(v) 1./(c.al*v); else, Af=@(v) (v-c.r)./(v.*(c.al*v+c.be)); end
            g=@(v) ex_state(c,v,Af(v))-a; u0=fzero(g,[b.u(j) b.u(j+1)]); A0=Af(u0);
        end
        [~,x]=ex_state(c,u0,A0);
        ev=eig(ex_jac(p,a,x));
        SS(end+1)=struct('x',x,'stable',all(real(ev)<0),'br',nm{i}); %#ok<AGROW>
    end
end
end

function [F,v] = ex_rhs_rates(x,AT,p)
% right-hand side and the eight reaction rates (for a scaled residual)
KT=p(1); PT=p(2); k1=p(3); k11=p(4); k2=p(5); k3=p(6); k33=p(7); k4=p(8);
k5=p(9); k55=p(10); k6=p(11); k7=p(12); k77=p(13); k8=p(14);
x=x(:); x1=x(1); x2=x(2); x3=x(3); x4=x(4); x5=x(5); x6=x(6);
A=AT-sum(x); K=KT-x3-x4; P=PT-x5-x6;
v=[k1*A*K, k11*x3, k2*x3, k3*x1*K, k33*x4, k4*x4, k5*x2*P, k55*x5, k6*x5, k7*x1*P, k77*x6, k8*x6];
v1=v(1)-v(2); v2=v(3); v3=v(4)-v(5); v4=v(6); v5=v(7)-v(8); v6=v(9); v7=v(10)-v(11); v8=v(12);
F=[v2-v3+v6-v7; v4-v5; v1-v2; v3-v4; v5-v6; v7-v8];
end

%% =====================================================================================
%%  LOCAL FUNCTIONS  --  figure helpers
%% =====================================================================================
function rgb = hex2rgb(h)
h=strrep(h,'#',''); rgb=[hex2dec(h(1:2)) hex2dec(h(3:4)) hex2dec(h(5:6))]/255;
end

function fprintf_both(fid,varargin)
fprintf(varargin{:}); fprintf(fid,varargin{:});
end

function fig = newfig(w,h)
fig=figure('Units','inches','Position',[0.5 0.5 w h],'Color','w','PaperPositionMode','auto','InvertHardcopy','off');
end

function save_fig(fig,out,name,fscale,gscale)
% MATLAB .fig (editable), vector PDF and 300-dpi PNG with the same base name as the Python figures
%  fscale : font factor of the legibility pass (1 = no legibility pass, figure saved as drawn)
%  gscale : canvas factor (the layout is kept; fonts grow by fscale/gscale relative to the panels)
if nargin<4 || isempty(fscale), fscale=1.45; end
if nargin<5 || isempty(gscale), gscale=1.2; end
pdf=fullfile(out,[name '.pdf']); png=fullfile(out,[name '.png']);
set(fig,'Visible','on');                          % so that the .fig opens visible
if fscale>1, beautify(fig,fscale,gscale); end     % box on, thicker lines, larger fonts, fit pass
savefig(fig,fullfile(out,[name '.fig']),'compact');
try
    exportgraphics(fig,pdf,'ContentType','vector','BackgroundColor','white');
    exportgraphics(fig,png,'Resolution',300,'BackgroundColor','white');
catch
    set(fig,'PaperUnits','inches'); pos=get(fig,'Position'); set(fig,'PaperSize',pos(3:4),'PaperPosition',[0 0 pos(3:4)]);
    print(fig,pdf,'-dpdf','-painters'); print(fig,png,'-dpng','-r300');
end
close(fig);
end

function beautify(fig,f,g)
% Legibility pass applied to every figure before it is saved:
%  (1) the canvas is enlarged by g, keeping the layout (all positions normalized);
%  (2) axes box on (line width 2), fonts x f (axes >= 10 pt), lines 2-3 pt, larger markers,
%      thin grid lines and no minor grid;
%  (3) fit pass: titles and axis labels wider than their panel are wrapped onto two lines
%      (then shrunk if needed); text that sat inside a panel is kept inside it (y-range of
%      bar charts extended, otherwise moved or shrunk); axes, labels and legends are kept
%      inside the figure and off each other.
drawnow;
tf=1+(f-1)*0.75;                                   % text inside the plots grows a little less
axs=findobj(fig,'Type','axes'); lgs=findobj(fig,'Type','legend');
set(axs,'Units','normalized'); if ~isempty(lgs), set(lgs,'Units','normalized'); end
set(fig,'Units','inches'); pos=fig.Position; fig.Position=[pos(1) pos(2) pos(3)*g pos(4)*g]; drawnow;
vis=axs(arrayfun(@(a) strcmp(a.Visible,'on'),axs));
% ---- remember which free texts lie inside their panel ------------------------------------
txt=findobj(fig,'Type','text'); inside=false(size(txt));
for i=1:numel(txt)
    pa=txt(i).Parent;
    if isgraphics(pa,'axes') && strcmp(pa.Visible,'on')
        e=text_ext(txt(i)); inside(i)=e(1)>=-0.01 && e(2)>=-0.01 && e(1)+e(3)<=1.01 && e(2)+e(4)<=1.01;
    end
end
% ---- axes: box, line width, fonts, grid ----------------------------------------------------
for ax=vis(:)'
    ti=ax.TightInset; p=ax.Position;
    op=[p(1)-ti(1) p(2)-ti(2) p(3)+ti(1)+ti(3) p(4)+ti(2)+ti(4)];
    lg=ax.Legend;
    if ~isempty(lg) && isvalid(lg) && contains(lg.Location,'outside')
        lp=lg.Position; x0=min(op(1),lp(1)); y0=min(op(2),lp(2));
        x1=max(op(1)+op(3),lp(1)+lp(3)); y1=max(op(2)+op(4),lp(2)+lp(4)); op=[x0 y0 x1-x0 y1-y0];
    end
    fx=ax.XLabel.FontSize; fy=ax.YLabel.FontSize; ft=ax.Title.FontSize;
    ax.FontSize=max(10,ax.FontSize*f);
    ax.XLabel.FontSize=max(ax.XLabel.FontSize,fx*f); ax.YLabel.FontSize=max(ax.YLabel.FontSize,fy*f);
    ax.Title.FontSize=max(ax.Title.FontSize,max(10.5,ft*f));
    ax.LineWidth=2; ax.TickLength=[0.015 0.015];
    if ~(ischar(ax.XColor) || ischar(ax.YColor)), ax.Box='on'; end   % not for axes with a hidden ruler
    ax.XMinorGrid='off'; ax.YMinorGrid='off';                         % minor grids get heavy with thick axes
    if isprop(ax,'GridLineWidth'), ax.GridLineWidth=0.8; end
    try
        ax.OuterPosition=op; ax.PositionConstraint='outerposition';
    catch
        ax.ActivePositionProperty='outerposition'; ax.OuterPosition=op;
    end
end
% ---- free text, legends ------------------------------------------------------------------
for t=txt(:)', t.FontSize=max(8,t.FontSize*tf); end
for lg=lgs(:)', lg.FontSize=min(12,max(9,lg.FontSize*f)); lg.ItemTokenSize=[20 12]; end
% ---- lines, markers, histograms, error bars, reference lines -----------------------------
for h=findobj(fig,'Type','line')'
    if ~strcmp(h.LineStyle,'none'), h.LineWidth=min(3,max(2,h.LineWidth*1.6)); end
    h.MarkerSize=h.MarkerSize*1.3;
end
for h=findobj(fig,'Type','constantline')', h.LineWidth=min(3,max(2,h.LineWidth*1.6)); end
for h=findobj(fig,'Type','errorbar')', h.LineWidth=2; end
for h=findobj(fig,'Type','histogram')'
    if strcmp(h.DisplayStyle,'stairs'), h.LineWidth=2.5; end
end
for h=findobj(fig,'Type','scatter')', h.SizeData=h.SizeData*1.8; end
for h=findobj(fig,'Type','rectangle')', if h.LineWidth<1.5 && ~isequal(h.EdgeColor,'none'), h.LineWidth=1.5; end, end
for h=findall(fig,'Type','arrowshape')', h.LineWidth=2; h.HeadLength=8; h.HeadWidth=8; end
drawnow;
% ---- fit (a): titles and axis labels no wider than their panel ----------------------------
for ax=vis(:)'
    fit_label(ax.Title,3,1.02); fit_label(ax.XLabel,3,1.05); fit_label(ax.YLabel,4,1.0);
end
% ---- fit (b): text that was inside a panel stays inside -----------------------------------
for i=find(inside(:))'
    t=txt(i); ax=t.Parent; e=text_ext(t); k=0;
    while e(2)+e(4)>0.99 && strcmp(t.Units,'data') && strcmp(ax.YScale,'linear') && strcmp(ax.YDir,'normal') && k<6
        yl=ax.YLim; ax.YLim=[yl(1) yl(1)+(yl(2)-yl(1))*(e(2)+e(4)+0.03)];   % head-room for bar labels
        drawnow; e=text_ext(t); k=k+1;
    end
    while (e(3)>0.97 || e(4)>0.97) && t.FontSize>7, t.FontSize=t.FontSize-0.5; drawnow; e=text_ext(t); end
    d=[0 0];
    if e(1)<0.01, d(1)=0.01-e(1); elseif e(1)+e(3)>0.99, d(1)=0.99-e(1)-e(3); end
    if e(2)<0.01, d(2)=0.01-e(2); elseif e(2)+e(4)>0.99, d(2)=0.99-e(2)-e(4); end
    if any(d)
        u=t.Units; t.Units='normalized'; t.Position(1:2)=t.Position(1:2)+d; t.Units=u;
    end
end
% ---- fit (c): panels inside the figure and off each other ---------------------------------
fit_layout(fig,vis,lgs);
end

function fit_layout(fig,vis,lgs)
m=0.008; n=numel(vis);
hid=findobj(fig,'Type','axes','Visible','off');   % e.g. a hidden axes that only carries a legend
for ax=vis(:)'
    try, ax.PositionConstraint='innerposition'; catch, ax.ActivePositionProperty='position'; end
end
for pass=1:5
    drawnow; B=zeros(n,4); P=B;
    for i=1:n
        ax=vis(i); p=ax.Position; ti=ax.TightInset;
        b=[p(1)-ti(1) p(2)-ti(2) p(1)+p(3)+ti(3) p(2)+p(4)+ti(4)];
        lg=ax.Legend;
        if ~isempty(lg) && isvalid(lg) && contains(lg.Location,'outside')
            l=lg.Position; b=[min(b(1),l(1)) min(b(2),l(2)) max(b(3),l(1)+l(3)) max(b(4),l(2)+l(4))];
        end
        B(i,:)=b; P(i,:)=p;
    end
    P0=P;
    % (i) inside the figure
    for i=1:n
        b=B(i,:);
        if b(1)<m,   d=m-b(1);     P(i,1)=P(i,1)+d; P(i,3)=P(i,3)-d; end
        if b(3)>1-m, P(i,3)=P(i,3)-(b(3)-(1-m)); end
        if b(2)<m,   d=m-b(2);     P(i,2)=P(i,2)+d; P(i,4)=P(i,4)-d; end
        if b(4)>1-m, P(i,4)=P(i,4)-(b(4)-(1-m)); end
    end
    % (ii) panels off each other (side by side -> narrower; stacked -> lower)
    for i=1:n
        for j=1:n
            if i==j, continue; end
            bi=B(i,:); bj=B(j,:);
            vo=min(bi(4),bj(4))-max(bi(2),bj(2)); ho=min(bi(3),bj(3))-max(bi(1),bj(1));
            if vo<=0 || ho<=0, continue; end
            if ho<=vo && bi(1)+bi(3)<bj(1)+bj(3)
                d=(ho+m)/2; P(i,3)=P(i,3)-d; P(j,1)=P(j,1)+d; P(j,3)=P(j,3)-d;
            elseif ho>vo && bi(2)+bi(4)<bj(2)+bj(4)
                d=(vo+m)/2; P(i,4)=P(i,4)-d; P(j,2)=P(j,2)+d; P(j,4)=P(j,4)-d;
            end
        end
    end
    % (iii) free-standing legends (placed by hand, or in a hidden axes) inside the figure and off the panels
    for lg=lgs(:)'
        if ~isvalid(lg), continue; end
        free=strcmp(lg.Location,'none') || any(arrayfun(@(a) ~isempty(a.Legend) && isequal(a.Legend,lg),hid));
        if ~free, continue; end
        l=lg.Position; l0=l;
        l(1)=min(max(l(1),m),1-m-l(3)); l(2)=min(max(l(2),m),1-m-l(4));
        if any(abs(l-l0)>1e-6), lg.Position=l; end
        for i=1:n
            b=B(i,:);
            vo=min(b(4),l(2)+l(4))-max(b(2),l(2)); ho=min(b(3),l(1)+l(3))-max(b(1),l(1));
            if vo<=0 || ho<=0, continue; end
            if l(2)+l(4)/2 < (b(2)+b(4))/2, P(i,2)=P(i,2)+vo+m; P(i,4)=P(i,4)-vo-m;   % legend below
            else, P(i,4)=P(i,4)-vo-m; end                                              % legend above
        end
    end
    P(:,3)=max(P(:,3),0.04); P(:,4)=max(P(:,4),0.04);
    if max(abs(P(:)-P0(:)))<1e-4, break; end
    for i=1:n, vis(i).Position=P(i,:); end
end
end

function fit_label(h,dim,lim)
% wrap a title / axis label onto two lines when it is wider (taller) than its panel, then shrink
if isempty(h.String), return; end
u=h.Units; h.Units='normalized'; drawnow;
if h.Extent(dim)>lim
    s=h.String; s2=wrap2(s);
    if ~isequal(s2,s), h.String=s2; drawnow; end
end
while h.Extent(dim)>lim && h.FontSize>9, h.FontSize=h.FontSize-0.5; drawnow; end
h.Units=u;
end

function s2 = wrap2(s)
% split a one-line (TeX) string at the space nearest its middle, never inside {...}
s2=s;
if iscell(s), if numel(s)~=1, return; end, s=s{1}; end
if ~ischar(s) || size(s,1)>1 || any(s==newline), return; end
depth=0; cand=[];
for i=2:numel(s)-1
    if s(i)=='{', depth=depth+1; elseif s(i)=='}', depth=depth-1;
    elseif s(i)==' ' && depth==0 && s(i-1)~=' ', cand(end+1)=i; end %#ok<AGROW>
end
if isempty(cand), return; end
[~,k]=min(abs(cand-numel(s)/2)); i=cand(k);
s2={strtrim(s(1:i-1)), strtrim(s(i+1:end))};
end

function e = text_ext(t)
% extent of a text object in normalized units of its axes
u=t.Units;
if ~strcmp(u,'normalized'), t.Units='normalized'; e=t.Extent; t.Units=u; else, e=t.Extent; end
end

function ptitle(ax,str)
title(ax,str,'FontWeight','normal','FontSize',8.5,'Color',[11 11 11]/255);
try, ax.TitleHorizontalAlignment='left'; catch, end
end

function s = commas(n)
s=sprintf('%d',round(n)); k=numel(s)-3;
while k>0 && ~(k==1 && s(1)=='-'), s=[s(1:k) ',' s(k+1:end)]; k=k-3; end
end

function s = fmtn(x)
if abs(x-round(x))<1e-9, s=sprintf('%d',round(x)); else, s=sprintf('%.1f',x); end
end

function s = sci(x)
e=floor(log10(abs(x))); s=sprintf('%.1f\\times10^{%d}',x/10^e,e);
end

function q = prctile_(x,p)
x=sort(x(~isnan(x))); n=numel(x); q=zeros(size(p));
for i=1:numel(p), t=1+(n-1)*p(i)/100; lo=floor(t); hi=min(lo+1,n); q(i)=x(lo)+(t-lo)*(x(hi)-x(lo)); end
end

function cm = lincmap(hexes,n)
if nargin<2, n=256; end
K=numel(hexes); P=zeros(K,3); for i=1:K, P(i,:)=hex2rgb(hexes{i}); end
cm=interp1(linspace(0,1,K),P,linspace(0,1,n));
end

function h = draw_stab(ax,p,AT,x2,x,col,lw)
% branch with stability: solid = stable, dashed = unstable (eigenvalues of the analytic Jacobian)
n=numel(AT); st=false(1,n);
for i=1:n, st(i)=all(real(eig(ex_jac(p,AT(i),x(:,i))))<0); end
y1=x2; y1(~st)=NaN; y2=x2; y2(st)=NaN;
h=plot(ax,AT,y1,'-','Color',col,'LineWidth',lw);
plot(ax,AT,y2,'--','Color',col,'LineWidth',lw);
end

function col = ucol(u)
% colour that encodes u = K/P (log scale, 1e-2 .. 1e5) in the schematic
cm=lincmap({'#1c5cab','#2a78d6','#1baf7a','#eda100','#eb6834','#e34948'});
t=(log10(u)+2)/7; t=min(max(t,0),1); col=cm(1+round(t*255),:);
end

function draw_cartoon(OUT,C)
% Fig EX.1: forward pipeline (a) and the inverse method (b); example S-dominant set
p=[114.78,17.251,6.328e-3,0.2131,1.543e-3,1.992e-2,8.245e-2,0.2087,2.945,1.688e-3,6.060e-3,1.080,7.265e-3,0.2224];
c=ex_consts(p); r=c.r; us=c.ustar; o=ex_analyse(p);
FW=7.5; FH=4.6; asp=FH/FW;
fig=newfig(FW,FH);
bg=axes(fig,'Position',[0 0 1 1],'XLim',[0 1],'YLim',[0 1],'Visible','off'); hold(bg,'on');
rbox=@(x,y,w,h,fc,ec) rectangle(bg,'Position',[x y w h],'Curvature',[0.06 0.1],'FaceColor',fc,'EdgeColor',ec,'LineWidth',0.7);
arr=@(x0,y0,x1,y1,col) annotation(fig,'arrow',[x0 x1],[y0 y1],'Color',col,'LineWidth',1,'HeadLength',5,'HeadWidth',5);
circ=@(x,y,rad,fc) rectangle(bg,'Position',[x-rad*asp y-rad 2*rad*asp 2*rad],'Curvature',[1 1],'FaceColor',fc,'EdgeColor','w','LineWidth',0.8);
T=@(x,y,s,varargin) text(bg,x,y,s,varargin{:});
% (a) forward pipeline
T(0,0.975,'(a)','FontSize',9,'FontWeight','bold','VerticalAlignment','top');
T(0.04,0.975,'MATCONT without the inverse method (single-start labels)','FontSize',8.5,'FontWeight','bold','VerticalAlignment','top');
steps={'Fix A_{TOT} = 0.1',{'Integrate ODEs','from X_0 = 0'},{'Continue in A_{TOT}','(single run)'},{'Label from the','branch followed'}};
bx=[0.01 0.205 0.40 0.595]; bw=0.165; by=0.80; bh=0.105;
for i=1:4
    rbox(bx(i),by,bw,bh,[0.957 0.953 0.937],C.grid);
    T(bx(i)+bw/2,by+bh/2,steps{i},'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',7.2);
    if i<4, arr(bx(i)+bw+0.004,by+bh/2,bx(i+1)-0.004,by+bh/2,C.muted); end
end
arr(bx(4)+bw+0.004,by+bh/2,0.79,by+bh/2,C.muted);
ax=axes(fig,'Position',[0.815 0.795 0.175 0.135]); hold(ax,'on');
b=o.br.main; k=b.AT<150; plot(ax,b.AT(k),b.x2(k),'-','Color',C.blue,'LineWidth',1.4);
b=o.br.det; k=b.AT<150; plot(ax,b.AT(k),b.x2(k),':','Color',C.muted,'LineWidth',0.9);
set(ax,'YScale','log','XTick',[],'YTick',[],'XColor',[0.6 0.6 0.56],'YColor',[0.6 0.6 0.56]);
xlabel(ax,'A_{TOT}','FontSize',6.5); ylabel(ax,'x_2','FontSize',6.5);
text(ax,0.98,0.52,'branch followed','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',6,'Color',C.blue);
text(ax,0.98,1.0,'detached branch (not reached)','Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',6,'Color',C.ink2);
% (b) inverse method
T(0,0.735,'(b)','FontSize',9,'FontWeight','bold','VerticalAlignment','top');
T(0.04,0.735,'MATCONT with the inverse method (exact labels)','FontSize',8.5,'FontWeight','bold','VerticalAlignment','top');
X=[0 0.2525 0.505 0.7575]; Wc=0.2375; Y=0.0; H=0.665;
for i=1:4, rbox(X(i)+0.002,Y+0.002,Wc-0.004,H,'w',hex2rgb('#b9d4f5')); end
for i=1:3, arr(X(i)+Wc+0.002,Y+H*0.5,X(i+1)-0.002,Y+H*0.5,C.blue); end
heads={'Choose u = K/P','Solve for A exactly','Rebuild the state','Place on the diagram'};
for i=1:4
    circ(X(i)+0.022,Y+H-0.045,0.022,C.blue);
    T(X(i)+0.022,Y+H-0.045,num2str(i),'HorizontalAlignment','center','Color','w','FontSize',7.5,'FontWeight','bold');
    T(X(i)+0.044,Y+H-0.045,heads{i},'FontSize',7.8,'FontWeight','bold');
end
foot={{'Each grid point is a candidate','steady state.'},{'Grid points with A < 0 (\times) are not','steady states; no root finding.'}, ...
      {'Conservation of substrate used in','reverse: A_{TOT} is an output.'},{'All branches appear, including the','detached one; folds are turning','points of A_{TOT}(u).'}};
for i=1:4, T(X(i)+0.012,Y+0.022,foot{i},'FontSize',6.6,'Color',C.ink2,'VerticalAlignment','bottom'); end
% step 1
T(X(1)+0.012,Y+H-0.10,{'ratio of free kinase to free','phosphatase; no X_0, no A_{TOT}'},'FontSize',6.6,'Color',C.ink2,'VerticalAlignment','top');
cy=Y+0.40;
circ(X(1)+0.055,cy,0.042,C.orange); T(X(1)+0.055,cy,'K','HorizontalAlignment','center','Color','w','FontSize',9,'FontWeight','bold');
T(X(1)+0.098,cy,':','HorizontalAlignment','center','FontSize',11,'Color',C.ink2);
circ(X(1)+0.135,cy,0.03,C.aqua); T(X(1)+0.135,cy,'P','HorizontalAlignment','center','Color','w','FontSize',8,'FontWeight','bold');
T(X(1)+0.175,cy,'= u','FontSize',9);
ax=axes(fig,'Position',[X(1)+0.02 Y+0.21 Wc-0.04 0.05]); hold(ax,'on');
ug=logspace(-2,5,22); for j=1:numel(ug), plot(ax,ug(j),0,'o','MarkerFaceColor',ucol(ug(j)),'MarkerEdgeColor','w','MarkerSize',4); end
set(ax,'XScale','log','XLim',[5e-3 2e5],'YLim',[-1 1],'YTick',[],'YColor','none','XTick',[1e-2 1e0 1e2 1e4],'FontSize',6.5);
xlabel(ax,'u (grid; colour = u throughout)','FontSize',6.5);
% step 2
T(X(2)+Wc/2,Y+H-0.125,'$A(u)=\frac{u-r}{u\,(\alpha u+\beta)}$','Interpreter','latex','HorizontalAlignment','center','FontSize',10);
ax=axes(fig,'Position',[X(2)+0.065 Y+0.155 Wc-0.08 0.27]); hold(ax,'on');
u=logspace(-2,5,3000); A=(u-r)./(u.*(c.al*u+c.be));
patch(ax,[r us us r],[1e-4 1e-4 1e5 1e5],[0.925 0.92 0.9],'EdgeColor','none');
m=A>0; plot(ax,u(m&u<r),A(m&u<r),'-','Color',C.grid,'LineWidth',1.2); plot(ax,u(m&u>us),A(m&u>us),'-','Color',C.grid,'LineWidth',1.2);
Ag=(ug-r)./(ug.*(c.al*ug+c.be)); ok=Ag>0;
for j=find(ok), plot(ax,ug(j),Ag(j),'o','MarkerFaceColor',ucol(ug(j)),'MarkerEdgeColor','w','MarkerSize',3.5); end
plot(ax,ug(~ok),3e-4*ones(1,sum(~ok)),'x','Color',C.muted,'MarkerSize',4);
text(ax,sqrt(r*us),1e-2,'A < 0','Rotation',90,'HorizontalAlignment','center','FontSize',6.3,'Color',C.ink2);
text(ax,r,1.4e5,'r','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',6.3,'Color',C.ink2);
text(ax,us,1.4e5,'u^*','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',6.3,'Color',C.ink2);
set(ax,'XScale','log','YScale','log','XLim',[1e-2 1e5],'YLim',[1e-4 1e5],'XTick',[1e-2 1e1 1e4],'YTick',[1e-3 1e1 1e5],'FontSize',6.5);
xlabel(ax,'u','FontSize',6.8); ylabel(ax,'free substrate A','FontSize',6.8);
% step 3
T(X(3)+0.012,Y+H-0.10,{'K, P and all substrate forms','follow from u and A'},'FontSize',6.6,'Color',C.ink2,'VerticalAlignment','top');
uu=1e3; AA=(uu-r)/(uu*(c.al*uu+c.be)); [~,x]=ex_state(c,uu,AA);
vals=[AA;x]; labs={'A','x_1','x_2','x_3','x_4','x_5','x_6'};
cols={'#cde2fb','#9ec5f4','#2a78d6','#ffd6c2','#f5a785','#b8ecd6','#6fd3a8'};
ax=axes(fig,'Position',[X(3)+0.03 Y+0.13 0.05 0.33],'Visible','off'); hold(ax,'on'); bot=0; tot=sum(vals);
for j=1:7
    hgt=max(vals(j)/tot,0.09); rectangle(ax,'Position',[0 bot 1 hgt-0.012],'FaceColor',hex2rgb(cols{j}),'EdgeColor','none');
    text(ax,1.15,bot+hgt/2,labs{j},'FontSize',6.5); bot=bot+hgt;
end
set(ax,'XLim',[0 1.7],'YLim',[0 bot]);
T(X(3)+0.145,Y+0.315,'sum','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',6.6,'Color',C.ink2);
arr(X(3)+0.115,Y+0.30,X(3)+0.175,Y+0.30,C.ink);
T(X(3)+0.205,Y+0.30,'A_{TOT}','HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
% step 4
ax=axes(fig,'Position',[X(4)+0.065 Y+0.17 Wc-0.08 0.36]); hold(ax,'on');
for nm={'main','det'}, b=o.br.(nm{1}); k=b.AT<160; plot(ax,b.AT(k),b.x2(k),'-','Color',C.grid,'LineWidth',1); end
ug4=logspace(-2,5,90); ug4=ug4(((ug4-r)./(ug4.*(c.al*ug4+c.be)))>0);
for uval=ug4
    Au=(uval-r)/(uval*(c.al*uval+c.be)); [ATu,xu]=ex_state(c,uval,Au);
    if ATu<160
        st=all(real(eig(ex_jac(p,ATu,xu)))<0); col=ucol(uval);
        if st, plot(ax,ATu,xu(2),'o','MarkerFaceColor',col,'MarkerEdgeColor','w','MarkerSize',3.5);
        else,  plot(ax,ATu,xu(2),'o','MarkerFaceColor','w','MarkerEdgeColor',col,'MarkerSize',3.5); end
    end
end
b=o.br.det; f=b.fidx(1); plot(ax,b.AT(f),b.x2(f),'p','MarkerFaceColor',C.ink,'MarkerEdgeColor',C.ink,'MarkerSize',6);
text(ax,4,b.x2(f)/60,'fold','FontSize',6.3);
bm=o.br.main; k=bm.AT<160; text(ax,155,max(bm.x2(k))*4,'main branch','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',6.3,'Color',C.ink2);
bd=o.br.det; k=bd.AT<160; text(ax,155,max(bd.x2(k))*1.6,'detached branch','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',6.3,'Color',C.ink2);
set(ax,'YScale','log','XLim',[0 160],'YLim',[1e-7 3e4],'YTick',[1e-6 1e-3 1e0 1e3],'XTick',[0 50 100 150],'FontSize',6.3);
xlabel(ax,'A_{TOT}','FontSize',6.8); ylabel(ax,'output x_2','FontSize',6.8);
h=[plot(ax,nan,nan,'o','MarkerFaceColor',C.muted,'MarkerEdgeColor',C.muted,'MarkerSize',3.5), ...
   plot(ax,nan,nan,'o','MarkerFaceColor','w','MarkerEdgeColor',C.muted,'MarkerSize',3.5), ...
   plot(ax,nan,nan,'p','MarkerFaceColor',C.ink,'MarkerEdgeColor',C.ink,'MarkerSize',5)];
legend(ax,h,{'stable','unstable','fold'},'Orientation','horizontal','Box','off','FontSize',6,'Location','south');
save_fig(fig,OUT,'EX_fig1_inverse_method',1);   % schematic drawn at print size: saved as drawn
end
