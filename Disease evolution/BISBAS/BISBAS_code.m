clearvars;
set(0,'DefaultFigureWindowStyle','docked');
T = readtable('DatiDaAnalizzare_Restricted_2024-06-10.xlsx');
T.BMI = T.WEIG./((T.HEIG/100).^2);
T.Duration = T.AGE-T.AGED;
T(T.BIS==0,:)=[];
T.GEND = categorical(T.GEND);
T.COUR = categorical(T.COUR);
eqCOUR=categorical(T.COUR, ...
	["Benigna" "Ricadute - remissioni" "Secondaria progressiva" "Progressiva con ricadute" "Primaria progressiva"], ...
	["RR" "RR" "SP" "SP" "PP"],'ordinal',true);
T=addvars(T,eqCOUR,'after','COUR');
T(isundefined(T.eqCOUR),:)=[];
%% nan sistemazione
w = sum(ismissing(T),2)>10; %1
T(w,:)=[];
T(isundefined(T.COUR),:)=[];
w2 = sum(ismissing(T))>15;T(:,w2)=[];
figure(1);clf;
imagesc(ismissing(T))
w3 = sum(ismissing(T))>0;w3f = find(w3);
for i=1:numel(w3f)
    T(isnan(T.(w3f(i))),w3f(i))=median(T(:,w3f(i)),'omitmissing');
end
imagesc(ismissing(T))

%% rescaling 
t = removevars(T,["CodicePaz_" "Data" "Paziente" "AGEO" "AGED" "AGEV" "WEIG" "HEIG"]);
for i=1:size(t,2)
    if string(class(t.(i)))=="double"
        t.(i)=(t.(i)-min(t.(i)))./(max(t.(i))-min(t.(i)));
    else
        t.(i)=t.(i);
    end
end
cut_off=(38-min(T.MFISTOT))./(max(T.MFISTOT)-min(T.MFISTOT));
%% definizine gruppi fatigue e fatigue free
fg = t.MFISTOT>cut_off;
FG = T(fg,:);
FF = T(~fg,:);
FG_norm = t(fg,:);FF_norm = t(~fg,:);
T = addvars(T,fg);

%% demografica di inizio 
demo_cat = ["GEND"   "COUR"];
demo = ["BIS" "BAS_Drive" "BAS_FunSeeking" "BAS_Global" "AGE"  "BMI"    "RELAPS"    "ISTRU"    "ABILHTOT"    "EDINBTOT"    "EDSSTOT2"    "FIM_TOT"    "HADSSUB1"    "HADSSUB2"    "LIFETOT"    "MOCATOT"    "OAB_QTOT"  "SDMTTOT"];

for i=1:numel(demo)
    [meanFG(i,2),meanFG(i,1)]= std(FG.(demo(i)),'omitmissing');
    [meanFF(i,2),meanFF(i,1)]= std(FF.(demo(i)),'omitmissing');
end
meanFG = array2table(meanFG);meanFG.Properties.VariableNames=["MEAN fg" "SD fg"];meanFG = addvars(meanFG,demo','Before','MEAN fg');
meanFF = array2table(meanFF);meanFF.Properties.VariableNames=["MEAN ff" "SD ff"];meanFF = addvars(meanFF,demo','Before','MEAN ff');
MEAN = [meanFG,meanFF(:,2:3)];
clear meanFG meanFF;

[FG_democat(1,:),~] = histcounts(FG.GEND);
[FG_democat(2,:),democat] = histcounts(FF.GEND);
FG_democat = array2table(FG_democat);FG_democat.Properties.VariableNames=democat;FG_democat=addvars(FG_democat,["FG";"FF"],'Before','F');

[FG_cour(1,:),~]=histcounts(FG.COUR);
[FG_cour(2,:),groupcat]=histcounts(FF.COUR);
FG_cour = array2table(FG_cour);FG_cour.Properties.VariableNames=groupcat;FG_cour=addvars(FG_cour,["FG";"FF"],'Before','Benigna');
DEMOCAT = [FG_democat,FG_cour(:,2:end)];

clear FG_cour FG_democat;
% [C,order] = confusionmat(T.GEND,T.fg);
[hTOT, chiTOT, pTOT]=crosstab(T.GEND,T.fg);
[mat, chi, p]=crosstab(T.COUR,T.fg);
% writetable(MEAN,'demo1.xlsx')
% writetable(DEMOCAT,'democat.xlsx')
%% correlazione semplice 

tcorr = removevars(t,demo_cat);
[rho, p_rho]=corr(table2array(tcorr));
f=figure(2);clf;
% tiledlayout('flow')
% ax(1)=nexttile;
imagesc(rho)
ax(1)=gca;
colormap(ax(1),'parula')
colorbar
axis square
ax(1).XTick=1:size(tcorr,2);ax(1).XTickLabel=tcorr.Properties.VariableNames;ax(1).XTickLabelRotation=45;ax(1).TickLabelInterpreter='none';
ax(1).YTick=1:size(tcorr,2);ax(1).YTickLabel=tcorr.Properties.VariableNames;
title(ax(1),'Pearson correlation \rho ')
% ax(2)=nexttile;
% imagesc(p_rho<0.05/numel(rho))
% colormap(ax(2),'bone')
% title(ax(2),'P value - Bonferroni correction')

% picturewidth=45;
% hw_ratio=0.5;
% set(findall(f,'-property','FontSize'),'FontSize',9)
% set(findall(f,'-property','Box'),'Box', 'off')
% % set(findall(ht,'-property','Interpreter'),'Interpreter','none')
% % set(findall(ht,'-property','TickLabelInterpreter'),'TickLabelInterpreter','none')
% set(f,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos=get(f,'Position');
% set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[pos(3) pos(4)])
% print(f, 'pearson','-dpdf','-painters','-fillpage');
% print(f, 'pearson','-dpng','-painters');
%% confronto variabili tra fg e ff
for i=1:size(tcorr,2)
    [h(i),p(i),ci(:,i),stats]=ttest2(table2array(tcorr(~fg,i)),table2array(tcorr(fg,i))); %unpaired samples 
end
pvalue = array2table(p');pvalue = addvars(pvalue, tcorr.Properties.VariableNames','Before','Var1');
pvalue = addvars(pvalue,repmat("not significant",size(pvalue,1),1),'After','Var1');ID = pvalue.Var1<0.05;
pvalue.Var3(ID)="SIG";
pvalue.Properties.VariableNames=["variable" "p value" "significance (.05)"];
% writetable(pvalue,'pvalue_fgff.xlsx')
%% linear regression su all sample 

ychosen = "MFISTOT";
xvariable = ["BIS" "BAS_Drive" "BAS_FunSeeking" "BAS_Global" "AGE"  "BMI"  "GEND"   "COUR"  "RELAPS"    "ISTRU"    "ABILHTOT"    "EDINBTOT"    "EDSSTOT2"    "FIM_TOT"    "HADSSUB1"    "HADSSUB2"    "LIFETOT"    "MOCATOT"    "OAB_QTOT"   "SDMTTOT"];
Y = t.(ychosen);
X = t(:,xvariable);

for i=1:size(X,2)
    unimdl{i}=fitlm((X(:,i)), Y);
    P_UNIVARIATA(i,3)=unimdl{i}.Coefficients.pValue(2);
    P_UNIVARIATA(i,1)=unimdl{i}.Coefficients.Estimate(2);
    P_UNIVARIATA(i,2)=unimdl{i}.Coefficients.SE(2);
end

P_UNIVARIATA=array2table(P_UNIVARIATA);P_UNIVARIATA.Properties.VariableNames=["beta" "SE" "p value"];
P_UNIVARIATA=addvars(P_UNIVARIATA,xvariable','Before','beta');

%% multilple linear regression 
xmulti = P_UNIVARIATA.Var1(P_UNIVARIATA.("p value")<0.05)';
multiMDL=fitlm(X(:,xmulti),Y);

P_MULTIVARIATA(:,1)=multiMDL.Coefficients.Estimate;
P_MULTIVARIATA(:,2)=multiMDL.Coefficients.SE;
P_MULTIVARIATA(:,3)=multiMDL.Coefficients.pValue;

P_MULTIVARIATA=array2table(P_MULTIVARIATA);P_MULTIVARIATA.Properties.VariableNames=["beta" "SE" "p value"];
P_MULTIVARIATA=addvars(P_MULTIVARIATA,["Intercept" xmulti]','Before','beta');

% Number of observations: 523, Error degrees of freedom: 507
% Root Mean Squared Error: 0.148
% R-squared: 0.53,  Adjusted R-Squared: 0.516
% F-statistic vs. constant model: 38.1, p-value = 2.89e-73


%% clustering k means ? 
xclust = ["BIS" "BAS_Drive" "BAS_FunSeeking" "BAS_Global" "HADSSUB1" "HADSSUB2"];
eval_pre = evalclusters(table2array(t(:,xclust)),'linkage',"CalinskiHarabasz","KList",1:10);
f=figure(3);clf;
plot(eval_pre)

%danno tutti 2...
 % vediamo se bisbas clusterizza in fg e ff
%% clusyering 
x = table2array(t(:,xclust));

linktype='ward';
distance='euclidean';
Z=linkage(x,linktype,distance); % matrice di pca reducted
distcoff=0.6*max(Z(:,3));
% vettore con label dei cluster per ogni paziente
C=clusterdata(x,'Criterion','distance','Cutoff',distcoff,'linkage',linktype,'distance',distance);
f=figure(4);clf;
ax=nexttile;
[hd,~,torder]=dendrogram(Z,0,'ColorThreshold',distcoff,'orientation','left');
xline(distcoff,'r');
ax.YTickLabel="";
ax.YDir='reverse';
xlabel('euclidean distance');
%%
for i=1:4
    for j=1:4
        % crtab = [sum(T.GEND(C==j)=='F'),sum(T.GEND(C==j)=='M');sum(T.GEND(C==i)=='F'),sum(T.GEND(C==i)=='M')];
        crtab = [sum(T.eqCOUR(C==j)=='RR'),sum(T.eqCOUR(C==j)=='SP'),sum(T.eqCOUR(C==j)=='PP');sum(T.eqCOUR(C==i)=='RR'),sum(T.eqCOUR(C==i)=='SP'),sum(T.eqCOUR(C==i)=='PP')];
        [h,p2gof(i,j),stats] = chi2cont(crtab);
      
    end
end
%% descrivo i cluster 

for i=1:max(C)
    d{i} = sum(C==i); %scopro la dimensione della cell i --> dipenda da quanta gente è nel cluster
    indixes_change = find(C == i); % trovo gli indici dei soggetti che appartengono all'i-esimo cluster
    phenotype{i}=t(indixes_change,:); % demografica dei cluster
end

means = zeros(size(t,2),2*max(C)); 
mediana = zeros(size(t,2),2*max(C)); 


for j = 1:max(C) % righe per cluster     
        for i = 1:size(t,2) % colonne per variabili
            if string(class(t.(i)))=="double"
                [means(i,2*j),means(i,2*j-1)]=std(table2array(phenotype{j}(:,i)));
                mediana(i,2*j-1)=median(table2array(phenotype{j}(:,i)));
                mediana(i,2*j)=iqr(table2array(phenotype{j}(:,i)));
           
            end
        end

end
means = array2table(means);means.Properties.VariableNames=["mean cluster 1" "SD cluster 1" "mean cluster 2" "SD cluster 2" "mean cluster 3" "SD cluster 3" "mean cluster 4" "SD cluster 4"];
means = addvars(means, t.Properties.VariableNames','Before','mean cluster 1');

mediana = array2table(mediana);mediana.Properties.VariableNames=["median cluster 1" "IQR cluster 1" "median cluster 2" "IQR cluster 2" "median cluster 3" "IQR cluster 3" "median cluster 4" "IQR cluster 4"];
mediana = addvars(mediana, t.Properties.VariableNames','Before','median cluster 1');

% writetable(means,'cluster_means.xlsx')
%% DIFFERENZE SIGNIFICATIVE FRA I CLUSTER ?

for i=1:size(t,2)
    if string(class(t.(i)))=="double"
       panova(i,2) = anova1(table2array(t(:,i)),C,'off');
    else
       panova(i,2)=0;
    end
end
panova = array2table(panova);panova.Properties.VariableNames=["variable" "p value (anova)"];
panova.variable=t.Properties.VariableNames';
xttest = panova.variable(panova.("p value (anova)")<0.05);

p_cluster = cell(numel(xttest),1);
for i=1:numel(xttest)
    if string(class(t.(xttest{i})))=="double"
    for j=1:max(C)
        for k=1:max(C)
            [~,p_cluster{i}(j,k)]=ttest2(table2array(phenotype{j}(:,xttest(i))),table2array(phenotype{k}(:,xttest(i))));
        end
    end
    end
    
end
% writecell(p_cluster,'pttest_cluster.xlsx')
%% vediamo nei vari cluster quanti affaticati ci sono 
f=figure(5);clf;
% ht = tiledlayout('flow');
cmap4=cbrewer2('Set3',4);
% for i = 1:max(C)
       % nexttile;
       bar(table2array(mean(phenotype{i}(:,xclust))),'FaceColor',cmap4(i,:))       
       ax=gca;
       ax.XTick=1:numel(xclust);ax.XTickLabel=xclust;ax.XTickLabelRotation=45;ax.TickLabelInterpreter='none';
       percent = round(d{i}/numel(C)*100);
       fatigue = round(sum(phenotype{i}.MFISTOT>cut_off)/d{i}*100);
       % title(['Group ', num2str(i),' (', num2str(d{i}),' pwMS and fatigue ',  num2str(fatigue) ' %)']);
       title(['Group ', num2str(2),' (', num2str(d{i}),' pwMS and fatigue ',  num2str(fatigue) ' %)']);

       ylabel('Mean value of each predictor');
ylim([0 0.8])
       % end
% title('Cluster')
%linkaxes

picturewidth=45;
hw_ratio=0.5;
set(findall(f,'-property','FontSize'),'FontSize',12)
set(findall(f,'-property','Box'),'Box', 'off')
% set(findall(ht,'-property','Interpreter'),'Interpreter','none')
% set(findall(ht,'-property','TickLabelInterpreter'),'TickLabelInterpreter','none')
set(f,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos=get(f,'Position');
set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[pos(3) pos(4)])
print(f, '2','-dpdf','-painters','-fillpage');
print(f, '2','-dpng','-painters');
%%
% category(:,3)=C;
% category =array2table(category);
% category.category3=str2double(category.category3);
% category.Properties.VariableNames=["Subject" "Real group" "cluster"];
% 
% % cluster 1 fatigue free
% acc_fg = sum(table2array(category(fg,"cluster"))==2)/sum(fg); % 60%
% acc_ff =  sum(table2array(category(~fg,"cluster"))==1)/sum(~fg); % 50
% acc = (sum(table2array(category(~fg,"cluster"))==1)+ sum(table2array(category(fg,"cluster"))==2))/numel(fg); % 55%
% 
% rocObj = rocmetrics(category.("Real group"),scores(:,1),names{1});
% plot(rocObj,ShowConfidenceIntervals=true) 
% 
% 
% writetable(means_dtc,'medie.xlsx')
% writematrix(sd_dtc,'SD.xlsx')





































%% fatigue VS non fatigue group 
% category = string(zeros(numel(fg),3));
% category(fg,2)="Fatigue Group";
% category(~fg,2)="Fatigue-free group";
% category(:,1)=T.CodicePaz_;
%% sistemazione table
% id_BISBAS = find(contains(T.Properties.VariableNames,'BIS') | contains(T.Properties.VariableNames,'BAS'));
% BISBAS = T(:,id_BISBAS);
% BISBASfg = T(fg,id_BISBAS);BISBASff = T(~fg,id_BISBAS);
% 
% descrittiva = zeros(numel(id_BISBAS),3);
% for i=1:numel(id_BISBAS)
%     [descrittiva(i,2),descrittiva(i,1)]=std(table2array(BISBAS(:,i)),'omitmissing');
%     descrittiva(i,3)=median(table2array(BISBAS(:,i)),'omitmissing');
% end
% descrittiva = array2table(descrittiva);descrittiva.Properties.VariableNames=["Mean" "SD" "Median"];
% descrittiva = addvars(descrittiva,BISBAS.Properties.VariableNames','Before','Mean');
% descrittiva.Var1=string(descrittiva.Var1);
% writetable(descrittiva,'descrittiva_demo.xlsx')




%% linear regression univariate nel fatigue group
% ychosen = "MFISTOT";
% xvariable = ["BIS" "BAS_Drive" "BAS_FunSeeking" "BAS_Global" "AGE"  "BMI"  "GEND"   "COUR"  "RELAPS"    "ISTRU"    "ABILHTOT"    "EDINBTOT"    "EDSSTOT2"    "FIM_TOT"    "HADSSUB1"    "HADSSUB2"    "LIFETOT"    "MOCATOT"    "OAB_QTOT"   "SDMTTOT"];
% Y = FG(:,contains(T.Properties.VariableNames,"MFIS"));
% X = FG(:,xvariable);
% 
% for i=1:size(Y,2)
%     Y.(i)=(Y.(i)-min(Y.(i)))./(max(Y.(i))-min(Y.(i)));
% end
% 
% for i=1:numel(xvariable)
%     if string(class(X.(xvariable(i))))=="double"
%         X.(i)=(X.(i)-min(X.(i)))./(max(X.(i))-min(X.(i)));
%     else
%         X.(i)=X.(i);
%     end
% end
% Xtot = [X,Y(:,ychosen)];
% Xtot.GEND=categorical(Xtot.GEND);Xtot.COUR=categorical(Xtot.COUR);
% % linear regression univariate 
% for i=1:size(Xtot,2)
%     FGmdl{i}=fitlm(table2array(Xtot(:,i)), Y.(ychosen));
%     P_UNIVARIATAfg(i)=FGmdl{i}.Coefficients.pValue(2);
%     beta_univariatafg(i)=FGmdl{i}.Coefficients.Estimate(2);
% end
% P_UNIVARIATAfg=array2table(P_UNIVARIATAfg');P_UNIVARIATAfg(end,:)=[];
% P_UNIVARIATAfg=addvars(P_UNIVARIATAfg,xvariable','Before','Var1');
% 
% xmulti = P_UNIVARIATAfg.Var1_1(P_UNIVARIATAfg.Var1<0.05);
% FGmultimdl=fitlm(table2array(Xtot(:,xmulti')), Y.(ychosen));
% resmultiFG=[FGmultimdl.Coefficients.Estimate,FGmultimdl.Coefficients.SE,FGmultimdl.Coefficients.pValue];
% resmultiFG=array2table(resmultiFG);resmultiFG=addvars(resmultiFG,["intercept"; xmulti],'Before','resmultiFG1')
% % writetable(resmultiFG,'multivar_fG.xlsx')
% %% linear regression univariate nei fatigue free
% Yff = FF(:,contains(T.Properties.VariableNames,"MFIS"));
% Xff = FF(:,xvariable);
% for i=1:size(Yff,2)
%     % if string(class(Y.(xvariable(i))))=="double"
%         Yff.(i)=(Yff.(i)-min(Yff.(i)))./(max(Yff.(i))-min(Yff.(i)));
%     % else
%     %     X.(i)=X.(i);
%     % end
% end
% 
% for i=1:numel(xvariable)
%     if string(class(Xff.(xvariable(i))))=="double"
%         Xff.(i)=(Xff.(i)-min(Xff.(i)))./(max(Xff.(i))-min(Xff.(i)));
%     else
%         Xff.(i)=Xff.(i);
%     end
% end
% 
% Xtot = [Xff,Yff(:,ychosen)];
% Xtot.GEND=categorical(Xtot.GEND);Xtot.COUR=categorical(Xtot.COUR);
% % FFmdl=fitlm(Xtot, ychosen);
% for i=1:size(Xtot,2)
%     FFmdl{i}=fitlm(table2array(Xtot(:,i)), Yff.(ychosen));
%     P_UNIVARIATAff(i)=FGmdl{i}.Coefficients.pValue(2);
%     beta_univariataff(i)=FGmdl{i}.Coefficients.Estimate(2);
% end
% P_UNIVARIATAff=array2table(P_UNIVARIATAff');P_UNIVARIATAff(end,:)=[];
% P_UNIVARIATAff=addvars(P_UNIVARIATAff,xvariable','Before','Var1');
% 
% xmulti = P_UNIVARIATAff.Var1_1(P_UNIVARIATAff.Var1<0.05);
% FFmultimdl=fitlm(table2array(Xtot(:,xmulti')), Yff.(ychosen));
% resmultiFF=[FFmultimdl.Coefficients.Estimate,FFmultimdl.Coefficients.SE,FFmultimdl.Coefficients.pValue];
% resmultiFF=array2table(resmultiFF);resmultiFF=addvars(resmultiFF,["intercept"; xmulti],'Before','resmultiFF1')
% 
% % writetable(resmultiFF,'multivar_fF.xlsx')
% 
% %% multivariata 
% % possibilità di usare come variabili esplicative la BIS e BAS per
% % spiegare la fatica percepita 
% xnames = ["Intercept" "BIS" "BAS_Drive" "BAS_FunSeeking" "BAS_Global"];
% Y = T(:,contains(T.Properties.VariableNames,"MFIS"));
% X = BISBAS(:,["BIS" "BAS_Drive" "BAS_FunSeeking" "BAS_Global"]);
% for i=1:size(Y,2)
%     % if string(class(Y.(xvariable(i))))=="double"
%         Y.(i)=(Y.(i)-min(Y.(i)))./(max(Y.(i))-min(Y.(i)));
%         X.(i)=(X.(i)-min(X.(i)))./(max(X.(i))-min(X.(i)));
%     % else
%     %     X.(i)=X.(i);
%     % end
% end
% 
% for i=1:size(Y,2)
%     [b(:,i),bint(:,2*i-1:2*i),r,rint,stats_multi(i,:)] = regress(table2array(Y(:,i)),[ones(size(X,1),1),table2array(X)]);
%     multimdl = fitlm(table2array(X),table2array(Y(:,i)))
% end
% b = array2table(b);b.Properties.VariableNames=Y.Properties.VariableNames;b = addvars(b,xnames','Before','MFISSUB1');
% stats_multi = array2table(stats_multi);stats_multi.Properties.VariableNames=["R^2" "F-score" "p value" "error variance estimate"];
% %% logistic? 
% % vedere se bis e bas classificano affatiacato vs non affaticato
% x = table2array(X);y=table2array(Y);
% [rho, p_rho]=corr(x, y);
% figure(1);clf;
% tiledlayout('flow')
% ax(1)=nexttile;
% imagesc(rho)
% colormap(ax(1),"parula")
% colorbar
% ax(1).XTick=1:4;ax(1).XTickLabel=X.Properties.VariableNames;
% ax(1).YTick=1:4;ax(1).YTickLabel=Y.Properties.VariableNames;
% ax(2)=nexttile;
% imagesc(p_rho<0.05/numel(rho))
% colormap(ax(2),"bone")
% ax(2).XTick=1:4;ax(2).XTickLabel=X.Properties.VariableNames;
% ax(2).YTick=1:4;ax(2).YTickLabel=Y.Properties.VariableNames;
% 
% %% clustering ?
% figure(1);clf;
% tiledlayout('flow');
% nexttile
% scatter(BISBAS.BIS, BISBAS.BAS_Global)
% [rho, p]=corrcoef(BISBAS.BIS, BISBAS.BAS_Global);
% title(sprintf('Linear correlation \\rho=%1.2f (p < .001)',rho(1,2)))
% ylabel('BAS')
% xlabel('BIS')
% nexttile
% scatter(BISBAS.BIS, BISBAS.BAS_Drive)
% ylabel('BAS DRIVE')
% xlabel('BIS')
% nexttile
% scatter(BISBAS.BIS, BISBAS.BAS_FunSeeking)
% ylabel('BAS FUNSEEKING')
% xlabel('BIS')
% nexttile
% scatter(BISBAS.BIS, BISBAS.BAS_RewardResp)
% ylabel('BAS REWARD RESPONSE')
% xlabel('BIS')
% 
% % not so prone to be clustered
% 
% figure(2);clf;
% tiledlayout('flow');
% nexttile
% scatter(BISBAS.BIS, T.MFISSUB1)
% % [rho, p]=corrcoef(BISBAS.BIS, BISBAS.BAS_Global);
% % title(sprintf('Linear correlation \\rho=%1.2f (p < .001)',rho(1,2)))
% ylabel('MFISSUB1')
% xlabel('BIS')
% nexttile
% scatter(BISBAS.BIS, T.MFISSUB2)
% ylabel('MFISSUB2')
% xlabel('BIS')
% nexttile
% scatter(BISBAS.BIS, T.MFISSUB3)
% ylabel('MFISSUB3')
% xlabel('BIS')
% nexttile
% scatter(BISBAS.BIS, T.MFISTOT)
% ylabel('MFIS')
% xlabel('BIS')
% 
% 
% figure(3);clf;
% tiledlayout('flow');
% nexttile
% scatter(BISBAS.BAS_Global, T.MFISSUB1)
% % [rho, p]=corrcoef(BISBAS.BIS, BISBAS.BAS_Global);
% % title(sprintf('Linear correlation \\rho=%1.2f (p < .001)',rho(1,2)))
% ylabel('MFISSUB1')
% xlabel('BaS')
% nexttile
% scatter(BISBAS.BAS_Global, T.MFISSUB2)
% ylabel('MFISSUB2')
% xlabel('BaS')
% nexttile
% scatter(BISBAS.BAS_Global, T.MFISSUB3)
% ylabel('MFISSUB3')
% xlabel('BaS')
% nexttile
% scatter(BISBAS.BAS_Global, T.MFISTOT)
% ylabel('MFIS')
% xlabel('BaS')
% 
% %% dividiamo in 4?
% rng(1);
% idC=kmeans(table2array(BISBAS(:,25:28)),4);
% 
% figure(4);clf;
% tiledlayout('flow');
% nexttile
% gscatter(BISBAS.BIS, BISBAS.BAS_Global,idC)
% [rho, p]=corrcoef(BISBAS.BIS, BISBAS.BAS_Global);
% title(sprintf('Linear correlation \\rho=%1.2f (p < .001)',rho(1,2)))
% ylabel('BAS')
% xlabel('BIS')
% nexttile
% gscatter(BISBAS.BIS, BISBAS.BAS_Drive,idC)
% ylabel('BAS DRIVE')
% xlabel('BIS')
% nexttile
% gscatter(BISBAS.BIS, BISBAS.BAS_FunSeeking,idC)
% ylabel('BAS FUNSEEKING')
% xlabel('BIS')
% nexttile
% gscatter(BISBAS.BIS, BISBAS.BAS_RewardResp,idC)
% ylabel('BAS REWARD RESPONSE')
% xlabel('BIS')
% 
% % not so prone to be clustered
% 
% figure(5);clf;
% tiledlayout('flow');
% nexttile
% gscatter(BISBAS.BIS, T.MFISSUB1,idC)
% % [rho, p]=corrcoef(BISBAS.BIS, BISBAS.BAS_Global);
% % title(sprintf('Linear correlation \\rho=%1.2f (p < .001)',rho(1,2)))
% ylabel('MFISSUB1')
% xlabel('BIS')
% nexttile
% gscatter(BISBAS.BIS, T.MFISSUB2,idC)
% ylabel('MFISSUB2')
% xlabel('BIS')
% nexttile
% gscatter(BISBAS.BIS, T.MFISSUB3,idC)
% ylabel('MFISSUB3')
% xlabel('BIS')
% nexttile
% gscatter(BISBAS.BIS, T.MFISTOT,idC)
% ylabel('MFIS')
% xlabel('BIS')
% 
% 
% figure(6);clf;
% tiledlayout('flow');
% nexttile
% gscatter(BISBAS.BAS_Global, T.MFISSUB1,idC)
% % [rho, p]=corrcoef(BISBAS.BIS, BISBAS.BAS_Global);
% % title(sprintf('Linear correlation \\rho=%1.2f (p < .001)',rho(1,2)))
% ylabel('MFISSUB1')
% xlabel('BaS')
% nexttile
% gscatter(BISBAS.BAS_Global, T.MFISSUB2,idC)
% ylabel('MFISSUB2')
% xlabel('BaS')
% nexttile
% gscatter(BISBAS.BAS_Global, T.MFISSUB3,idC)
% ylabel('MFISSUB3')
% xlabel('BaS')
% nexttile
% gscatter(BISBAS.BAS_Global, T.MFISTOT,idC)
% ylabel('MFIS')
% xlabel('BaS')
% 


