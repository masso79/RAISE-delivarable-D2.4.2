clearvars;
T=readtable('PROMOPROMS_EQUI.xlsx');
set(0,'DefaultFigureWindowStyle','docked')
%% Descrizione del campione 

% costruzione tabella demografica
DEMO = {'AGE' 'DISDUR_Diagnosis' 'ISTRU' 'EDSSTOT2'};
CAT = {'GEND' 'COUR'};
EQUI = {'SOM' 'VIS' 'VEST' 'PREF' 'COMPOSITE'};
Tdemo = T(:,[DEMO CAT]);
% mean e SD delle var demografiche numeriche con relativo intervallo di min
% e max
demo = zeros(size(DEMO,2),2);MIN = zeros(size(DEMO,2),2);
for i=1:size(DEMO,2)
    [demo(i,2),demo(i,1)]=std(Tdemo.(i),'omitnan');
    MIN(i,:) = minmax(Tdemo.(i)');
end
demo=[demo,MIN];clear MIN;
nameDEMO = Tdemo.Properties.VariableNames';
demo = array2table(demo); demo = addvars(demo,nameDEMO(1:4),'Before','demo1');
newNames = ["Vars","Mean","SD","Minimum", "Maximum"];
demo = renamevars(demo,1:5,newNames);
% categorical variables
cat{1} = histcounts(categorical((Tdemo.GEND)));
[cat{2},cour] = histcounts(categorical((Tdemo.COUR)));

% patient reported outcomes
PRO = {'ABILHTOT' 'FIM_SUB1' 'FIM_SUB2' 'FIM_SUB3' 'FIM_SUB4' 'FIM_SUB5'    'FIM_SUB6'    'FIM_TOT'    'HADSSUB1'   'HADSSUB2'    'HADSTOT' 'LIFETOT'   'MFISSUB1'    'MFISSUB2'    'MFISSUB3'    'MFISTOT'    'MOCA001'    'MOCA002'   'MOCA003' 'MOCA004'    'MOCA005'    'MOCA006' 'MOCA007'    'MOCA008'    'MOCA009'    'MOCA010'   'MOCA011'    'MOCATOT'    'OAB_QTOT'    'PASATTOT'    'SDMTTOT'};
Tpro=T(:,PRO);
demoPRO=zeros(size(PRO,2),2);MIN = zeros(size(PRO,2),2);
for i=1:size(PRO,2)
    [demoPRO(i,2),demoPRO(i,1)]=std(Tpro.(i),'omitnan');
    MIN(i,:) = minmax(Tpro.(i)');
end
demoPRO=[demoPRO,MIN];clear MIN;
namePRO = Tpro.Properties.VariableNames';
demoPRO = array2table(demoPRO); demoPRO = addvars(demoPRO,namePRO(1:end),'Before','demoPRO1');
demoPRO = renamevars(demoPRO,1:5,newNames);
T.PASATTOT(isnan(T.PASATTOT))=median(T.PASATTOT,'omitmissing');
%% Kolmogorov Smirnov test for normality

% non paramteric test to see if data follows a given distribution (nor only
% normality)
tF = T(categorical(T.GEND)=='F',:); % matrice con tutte le info delle sole donne
tM = T(categorical(T.GEND)=='M',:); % matrice con tutte info dei soli uomini
tRR = T(categorical(T.COUR)=='Ricadute - remissioni',:); % matrice con tutte info dei RR
tPP = T(categorical(T.COUR)=='Primaria progressiva ',:); % matrice con tutte info dei PP
tSP = T(categorical(T.COUR)=='Secondaria progressiva ',:); % matrice con tutte info SP
tP = [tPP;tSP];
tcell = {tF;tM;tRR;tPP;tSP;tP};
% Valutazione della normalità dei risultati di equitest delle catgoriche
for i=1:size(EQUI,2)
    [h(i,1),p(i,3)] = kstest(tF.(EQUI{i}));
    [h(i,2),p(i,4)] = kstest(tM.(EQUI{i}));
    [~,p(i,5)] = kstest(tRR.(EQUI{i}));
    [~,p(i,6)] = kstest(tPP.(EQUI{i}));
    [~,p(i,7)] = kstest(tSP.(EQUI{i}));
    [~,p(i,8)] = kstest(tP.(EQUI{i}));
end
% Valutazione della normalità dei risultati di equitest di tutti
for i=1:size(EQUI,2)    
    [~,pall(i)] = kstest(T.(EQUI{i}));
    
end
% Valutazione della normalità dei PROs
for i=1:size(PRO,2)  
    [~,pallPRO(i)] = kstest(T.(PRO{i}));
end
% tutti normali 
% tutti normalmente distribuiti
categ = ["F","M","RR","PP","SP","PROGRESSIVO"];
meanRES = zeros(size(EQUI,2),numel(categ)*2);
% medie per gruppi
for i=1:size(EQUI,2) % sui 5 test dell'equitest
    for j=1:numel(categ)
        [meanRES(i,j*2),meanRES(i,j*2-1)]=std(tcell{j}.(EQUI{i}),'omitmissing');
    end
end

meanRES = array2table(meanRES);
meanRES = addvars(meanRES,EQUI','Before','meanRES1');clear categ;

for i=1:5
    meanRES.Properties.VariableNames(2*i+1)=EQUI(i);
end
meanRES.Properties.VariableNames=["Categories","Mean F","SD F","Mean M","SD M","Mean RR","SD RR","Mean PP","SD PP","Mean SP","SD SP","Mean PROGRESSIVO","SD PROGRESSIVO"];

% t test
pt = zeros(size(EQUI,2),5);
for i=1:size(EQUI,2)    
    [~,pt(i,1)] = ttest2(tcell{1}.(EQUI{i}),tcell{2}.(EQUI{i})); 
    [~,pt(i,2)] = ttest2(tRR.(EQUI{i}),tPP.(EQUI{i})); 
    [~,pt(i,3)] = ttest2(tRR.(EQUI{i}),tSP.(EQUI{i})); 
    [~,pt(i,4)] = ttest2(tSP.(EQUI{i}),tPP.(EQUI{i})); 
    [~,pt(i,5)] = ttest2(tRR.(EQUI{i}),tP.(EQUI{i}));
end
pt = array2table(pt);pt = addvars(pt,EQUI','Before','pt1');
pt = renamevars(pt,1:6,["Equitest res","M/F","RR/PP","RR/SP","PP/SP","RR/PROG"]);

%% Pearson correlation

rho = zeros(size(PRO,2),size(EQUI,2));pval = zeros(size(PRO,2),size(EQUI,2));
T.PREF(isnan(T.PREF))=median(T.PREF,'omitmissing');
for i=1:size(EQUI,2)
    for j=1:size(PRO,2)
        [rho(j,i),pval(j,i)] = corr(T.(EQUI{i}),T.(PRO{j}));
    end
end
p_thresh=pval<0.05;
p_threshB=pval < 0.05/numel(rho);
f=figure(1);clf;
subplot(1,3,1);imagesc(rho);axis image;xlabel('EQUITEST');ylabel('PROs');
ax=gca;
ax.XTick=1:5;ax.XTickLabel=string(EQUI);ax.XTickLabelRotation=45;
ax.YTick=1:31;ax.YTickLabel=string(PRO);
ax.TickLabelInterpreter='none';
colorbar
title('Pearson correlation - \rho')

subplot(1,3,2);imagesc(p_thresh);axis image;xlabel('EQUITEST');ylabel('PROs');
ax=gca;
ax.XTick=1:5;ax.XTickLabel=string(EQUI);ax.XTickLabelRotation=45;
ax.YTick=1:31;ax.YTickLabel=string(PRO);
ax.TickLabelInterpreter='none';

title('Pearson correlation - P value')


subplot(1,3,3);imagesc(p_threshB);axis image;xlabel('EQUITEST');ylabel('PROs');
ax=gca;
ax.XTick=1:5;ax.XTickLabel=string(EQUI);ax.XTickLabelRotation=45;
ax.YTick=1:31;ax.YTickLabel=string(PRO);
ax.TickLabelInterpreter='none';

title('Pearson correlation - P value (Bonferroni correction)')

%% Fisher exact test 
% [hTOT, chiTOT, pTOT]=crosstab(Realclasses,predictions)
Fpv=zeros(size(PRO,2),size(EQUI,2));
for i=1:size(EQUI,2)
   equi_coff(i)=median(T.(EQUI{i}));
   equibin=discretize(T.(EQUI{i}),[-inf, equi_coff(i),Inf],'categorical'); % binarizzare l'informazione (categorizzare)
    for j=1:size(PRO,2)
        pro_coff(j)=median(T.(PRO{j}));
        probin=discretize(T.(PRO{j}),[-inf, pro_coff(j),Inf],'categorical');
        % conttlb=crosstab(equibin,probin);
        conttlb=crosstab(equibin,probin)
        [~,Fpv(j,i)]=fishertest(conttlb);
    end
end

figure(7);clf;
imagesc(-log10(Fpv));axis image;xlabel('AMY');ylabel('TAU'); % valori del tipo 0.0003 che è 3 per 10 alla -4, quindi rimane solo l'esposnente, sopra il 2 è già bene
ax=gca;
ax.YTick=1:numel(PRO);ax.YTickLabel=PRO;
ax.XTick=1:numel(EQUI);ax.XTickLabel=EQUI;ax.XTickLabelRotation=45;
ax.TickLabelInterpreter='none';
title('Fisher test p-value [-log_{10}]','FontWeight','normal');
colormap("hot");
colorbar;
ax.CLim=[2,3.5];

p_threshB=Fpv < 0.01;
%% scatter plot delle coppie di coordinate significative per Bonferroni
f=figure(4);clf;
tiledlayout('flow')
for i=1:size(EQUI,2)
    for j=1:size(PRO,2)
        if p_threshB(j,i)==1
            nexttile;
            scatter(T.(EQUI{i}),T.(PRO{j}),'filled')
            xlabel(EQUI{i},'Interpreter','none')
            ylabel(PRO{j},'Interpreter','none')
            title(sprintf('Linear correlation \\rho = %1.2f  [%1.1e]',rho(j,i),pval(j,i)),'fontweight','normal','fontsize',11)
        end
    end
end
T.PREF(isnan(T.PREF))=median(T.PREF,'omitmissing');
%% CLUSTER SULLA BASE DEI RISULTATI EQUI TEST
% valuto dai dati i numero ideale di cluster (massimo valore dell'indice)
eva = evalclusters(table2array(T(:,[EQUI(1:4)])),'kmeans','CalinskiHarabasz','KList',1:6);
figure(2);clf;
plot(eva)

% faccio i cluster con il valore valutato dal calinsky harabasz
k=eva.OptimalK; % viene un cluster da 1 persona
rng default; % For reproducibility
Tclust = T(:,[EQUI(1:4)]);
idC=kmeans(table2array(Tclust),k,'OnlinePhase','on','Replicates',25);
for i =1:max(idC)
    sz_cluster(i,1) = sum(idC==i);
end
Tclust = addvars(Tclust,idC);
%% valuto quale feature ha pesato di più 
figure(3);clf;
tiledlayout('flow')
mod = TreeBagger(150,Tclust,'idC',Surrogate="on",OOBPredictorImportance="on");
nexttile;
plot(oobError(mod))
title('Error in function of the # of trees')
[imp,ord] = sort(mod.OOBPermutedPredictorDeltaError); 
ax=nexttile;
bar(imp)
ax.XTick=1:4;ax.XTickLabel=Tclust.Properties.VariableNames(ord);ax.XTickLabelRotation=45;
title('Feature importance')
grid on;

%% valuto la composizione descrittiva dei cluster
jbar = find(sz_cluster==1); % il clusyer j bar non va
T(idC==jbar,:)=[];idC(idC==jbar)=[];

Tc = T(:,[CAT,DEMO,PRO,EQUI]);
for i=1:numel(unique(idC))
    for j=1:size(idC)
        if idC(j)==4
            idC(j)=3;
        else
            idC(j)=idC(j);
        end
    end
end
equi = cell(max(idC),1);impaired = cell(max(idC),1);
for j=1:1:max(idC)
    equi{j}=Tc(idC==j,:);
    impaired{j}=T(idC==j,["SOM_Impaired" "VIS_Impaired" "VEST_Impaired" "TOT_Impaired"]);
end
%% FOCUS SU SESSO E COURSE
X = crosstab(idC,T.GEND);
clear h;
for i=1:max(idC)
    for j=1:max(idC)
        if j>i
            [h(i,j),pX(i,j),~] = fishertest([X(i,:);X(j,:)]);
        end
    end
end

T.COUR = string(T.COUR);
for i=1:size(T.COUR,1)
    if T.COUR(i)=="Secondaria progressiva"
        T.COUR(i)="Primaria progressiva";
    else 
        T.COUR(i)=T.COUR(i);
    end
end
x2 = crosstab(idC,T.COUR);
% T2COUR = T.COUR()
% x2 = [X2(:,1),X2(:,3)];
for i=1:max(idC)
    for j=1:max(idC)
        if j>i
            [h2(i,j),pX2(i,j),~] = fishertest([x2(i,:);x2(j,:)]);
        end
    end
end

%%


% medie e dispersione per i tre cluster
clustSD = table(max(idC),size(Tc,2));clustM = table(max(idC),size(Tc,2));
for j=1:numel(equi)
    for i=3:size(Tc,2)
        [clustSD(j,i),clustM(j,i)]=std(equi{j}(:,i));        
    end
end

profile = Tc.Properties.VariableNames;
clustSD.Properties.VariableNames=profile;
clustM.Properties.VariableNames=profile;
panova= zeros(size(Tc,2),1);
for i=3:size(Tc,2)
    panova(i) = anova1(table2array(Tc(:,i)),idC,'off');
end
panova_tresh=panova<0.05/3; % Dunn-Bonferroni correction
panova = array2table(panova);
panova=addvars(panova,string(Tc.Properties.VariableNames'),'before','panova');

panova=addvars(panova,panova_tresh,'After','panova');
name_sig=panova.Var1(panova_tresh(1:end-5));
PVAL = panova.panova(panova_tresh);
clear mat

for i=1:size(clustM,2)
    for j=1:size(name_sig,1)
        if clustM.Properties.VariableNames{i}==name_sig(j)
            mat{j,1}=clustM.(i)';
            matSD{j,1}=clustSD.(i)';
         
            % resume = cat(1,resume,clustM.(i)');
           % resume(j,1)=string(clustM.Properties.VariableNames{i});
           % resume(j,2:4)=clustM.(i)';
        end
    end
end

mat = cell2mat(mat);matSD = cell2mat(matSD);
Mat = array2table(mat);Mat = addvars(Mat,name_sig,'Before','mat1');
MATsd = array2table(matSD);
Mat = addvars(Mat,MATsd.matSD1,'after','mat1');
Mat = addvars(Mat,MATsd.matSD2,'after','mat2');
Mat = addvars(Mat,MATsd.matSD3,'after','mat3');

% for i=1:numel(CAT)
%     cat_clust=T.(CAT{i})
%     
% end

%% POST HOC ANALYSIS
pt = cell(numel(name_sig),2);
for i=1:max(idC)
    for j=1:max(idC)
        for k=3:numel(name_sig)
            [~,pt{k,1}(i,j)] = ttest2(equi{i}.(name_sig{k}),equi{j}.(name_sig{k}));
        end
    end
end
for k=3:numel(name_sig)
    pt{k,2}=name_sig(k);
end

%% VISUALIZZIAMO QUANTI IMPAIRED NEI 3 CLUSTER

f=figure(6);clf;
ht = tiledlayout(3,numel(equi));
for j=1:numel(equi) % cluster
    ax(j)=nexttile;
    for i=1:3 % numer feat clust
        prob_imp = sum(table2array(impaired{j}(:,i)))/size(impaired{j},1);
        bar(i,prob_imp);axis equal square;
        hold on
       
    end   
    h = gca;
    h.XTick=1:3;
    h.XTickLabel = T(:,EQUI).Properties.VariableNames(1:3);
    h.XTickLabelRotation = 45;
    h.TickLabelInterpreter = "none";
    ylabel('Probabilità di essere impaired')
    title(ax(j),['cluster ',num2str(j)])
    hold off
end
linkaxes([ax(1) ax(2) ax(3)],'y')

for j=1:numel(equi) % cluster
    ax(j+3)=nexttile;  
    % bar(1,mean(equi{j}.PREF));axis equal square;
    hold on;
    bar(1,mean(equi{j}.COMPOSITE));axis equal square;
    ylabel('mean COMPOSITE score')
    h = gca;
    h.XTick=1;
    h.XTickLabel = T(:,EQUI).Properties.VariableNames(5);
    h.XTickLabelRotation = 45;
    h.TickLabelInterpreter = "none";
end
linkaxes([ax(4) ax(5) ax(6)],'y')

for j=1:numel(equi) % cluster
    ax(j+6)=nexttile;  
    % bar(1,mean(equi{j}.PREF));axis equal square;
    hold on;
    bar(1,mean(equi{j}.PREF,'omitnan'));axis equal square;
    ylabel('mean score PREF')
    h = gca;
    h.XTick=1;
    h.XTickLabel = T(:,EQUI).Properties.VariableNames(4);
    h.XTickLabelRotation = 45;
    h.TickLabelInterpreter = "none";
end
linkaxes([ax(7) ax(8) ax(9)],'y')
% picturewidth=35;
% hw_ratio=0.4;
% set(findall(f,'-property','FontSize'),'FontSize',11)
% set(findall(f,'-property','Box'),'Box', 'off')
% % set(findall(ht,'-property','Interpreter'),'Interpreter','none')
% % set(findall(ht,'-property','TickLabelInterpreter'),'TickLabelInterpreter','none')
% set(f,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos=get(f,'Position');
% set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[pos(3) pos(4)])
% print(f, 'clustprefALL','-dpdf','-painters','-fillpage');
% print(f, 'clustprefALL','-dpng','-painters');

%% other comorbidities

cat_morbius=categorize_comorbidities(T.OTHERDIS);
impairment=sum(cat_morbius,2);
T.impairment=categorical(impairment.sum,0:4,{'none' 'low' 'medium' 'high' 'severe'},'ordinal',true);
figure(5);clf;
tiledlayout(2,numel(unique(idC)));
for j=1:numel(equi) % cluster
    ax(j)=nexttile;
    for i=1:3 % numer feat clust
        prob_imp = sum(table2array(impaired{j}(:,i)))/size(impaired{j},1);
        bar(i,prob_imp);axis equal square;
        hold on       
    end 
    grid on;
    h = gca;
    h.XTick=1:3;
    h.XTickLabel = T(:,EQUI).Properties.VariableNames(1:3);
    h.XTickLabelRotation = 45;
    h.TickLabelInterpreter = "none";
    ylabel('Probabilità di essere impaired')
    title(ax(j),['cluster ',num2str(j)])
    hold off
end
linkaxes([ax(1) ax(2) ax(3)],'y')

for j=1:numel(unique(idC))
    Tw = T(idC == j,:);
    ax(j+3) = nexttile;
    histogram(Tw.impairment,'Normalization','probability');
    grid on;
    ylabel('Normalised frequency of comorbidities')
    title(ax(j),sprintf('cluster %i',j))
end
linkaxes([ax(4) ax(5) ax(6)],'y')
