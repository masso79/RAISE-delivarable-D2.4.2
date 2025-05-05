%% table loads and defaults

set(0,'DefaultFigureWindowStyle','docked');
Nvars = 10;
datarange='A2:J284';
varrange='A1:J1';

% import data for time points 1 and 2 with standardised values + non
% standardised
opts = spreadsheetImportOptions('Sheet','T1_sd','NumVariables', Nvars,'datarange',datarange,'VariableNamesRange',varrange,'VariableTypes',repmat({'double'},[1,Nvars]));
T1_sd = readtable('PHENOTYPE.xlsx',opts);

opts = spreadsheetImportOptions('Sheet','T2_sd','NumVariables', Nvars,'datarange',datarange,'VariableNamesRange',varrange,'VariableTypes',repmat({'double'},[1,Nvars]));
T2_sd = readtable('PHENOTYPE.xlsx',opts);

opts = spreadsheetImportOptions('Sheet','T1','NumVariables', Nvars,'datarange',datarange,'VariableNamesRange',varrange,'VariableTypes',repmat({'double'},[1,Nvars]));
T1 = readtable('PHENOTYPE.xlsx',opts);

opts = spreadsheetImportOptions('Sheet','T2','NumVariables', Nvars,'datarange',datarange,'VariableNamesRange',varrange,'VariableTypes',repmat({'double'},[1,Nvars]));
T2 = readtable('PHENOTYPE.xlsx',opts);
%%
opts = spreadsheetImportOptions('Sheet','matrice','NumVariables',130,'datarange','A2:DZ284','VariableNamesRange','A1:DZ1');
tot = readtable('PHENOTYPE.xlsx',opts);


%% rename vars standardised

newNames = ["VISUO EXEC.","NAMING", "ATTENTION" ,"LANGUAGE","ABSTRACTION","DELAYED RECALL","ORIENTATION","IPS","ANXIETY","DEPRESSION"];

T1_sd = renamevars(T1_sd,1:width(T1_sd),newNames);
T2_sd = renamevars(T2_sd,1:width(T2_sd),newNames);
T1 = renamevars(T1,1:width(T1),newNames);
T2 = renamevars(T2,1:width(T2),newNames);

variations = T1.Properties.VariableNames;

%% valori mancanti standardized
% individuazione di valori nan nella table T1 ed eliminazione dei pazienti
% con valori assenti
a = isnan(table2array(T1_sd));pat_thresh = sum(a,2)>0;

% matrici con 272 pazienti (da 283)
t1_sd = table2array(T1_sd(~pat_thresh,:)); t2_sd = table2array(T2_sd(~pat_thresh,:)); t1=table2array(T1(~pat_thresh,:));t2=table2array(T2(~pat_thresh,:));

%tabelle con 272 pazienti
T1final_sd = T1_sd(~pat_thresh,:); T2final_sd = T2_sd(~pat_thresh,:); T1final = T1(~pat_thresh,:);T2final = T2(~pat_thresh,:);
tot = tot(~pat_thresh,:);

%%

tot.GEND=categorical(tot.GEND);
histogram(tot.GEND)
%% silhouette
%valutazione better number of clusters from data 
figure(1);clf;
evaluation = evalclusters(t1_sd,"kmeans","silhouette","KList",2:6);
plot(evaluation)
title("Silhouette indix for k=2,3,4,5,6 clusters")
% clusters = evaluation.OptimalY;

%% pca + clustering 

[pc,sco,~,~,expl]=pca(zscore(t1_sd)); % pca richiede zscore, e fare z score su dati ordinali non è corretto
maxeig=find(cumsum(expl)<85,1,'last');

% [iso,mapping]=compute_mapping(zscore(t1_sd), 'Isomap',maxeig);


%clustering su pca coordinates (solo a scopo rappresentativo!!)
figure(2);clf;
kT1_pca=kmeans(sco(:,1:maxeig),4,'OnlinePhase','on','Replicates',10);
% kiso=kmeans(iso(:,1:maxeig),4,'OnlinePhase','on','Replicates',10);

hg1=gscatter(sco(:,1),sco(:,2),kT1_pca);axis equal
title("Clustering con k da silhouette evaluation")
xlabel("PC1");ylabel("PC2");

%clustering pazienti al tempo t1 standard
kT1_sd = kmeans(t1_sd,4,'OnlinePhase','on','Replicates',10); % ok fare cluster qua? meglio sui non?
kT1_tab = array2table(kT1_sd);
%% second control about silhouette an number of cluster choice

figure(3);clf;
silhouette(t1,kT1_sd)

%% clustering phenotypes t1

cmap4=cbrewer2('Set3',10);

matrix = zeros(size(t1,1),Nvars+1); % matrice di 272 x 11, con aggiunta colonna indici cluster
matrix(:,1)= kT1_sd; % cluster choice for the subject
matrix(:,2:Nvars+1)=t1; % subjects with not standardised values

phenotype = cell(1,4);
CI = cell(1,4); % ogni cella sarà un vettore di 10 elementi 
MD = cell(1,4);
d = cell(1,4);

for i=1:4
d{i} = sum(ismember(matrix(:,1),i)); %scopro la dimensione della cell i --> dipenda da quanta gente è nel cluster
indixes_change = find(ismember(matrix(:,1),i)); % trovo gli indici dei soggetti che appartengono all'i-esimo cluster
phenotype{i}=zeros(d{i},10); % creo matrice con 10 colonne e d righe, quindi una matrice per ogni cluster
phenotype{i}=matrix(indixes_change,2:11); % riempio la matrice con i soggetti giusti
impairments_threshold_under = [4 3 6 3 1 4 6 35]; impairments_threshold_over = [7 7]; % in ordine di new names
for j=1:8
    CI{i}(:,j) = phenotype{i}(:,j)<impairments_threshold_under(1,j);  % cognitive impairment
    for l=1:2
    MD{i}(:,l) = phenotype{i}(:,l+8)>impairments_threshold_over(1,l); % mood disorder
    end
end
end

% merge 
total_imp = cell(1,4); %matrici di 0 e 1 per ogni soggetto del cluster impaired nella feature
SOMME = cell(1,4); %totale di soggetti impaired i quela feature indiepndente da chi è 
p = cell(1,4); %probabilità per ogni colonna di avere soggetto impaired nel fenotipo 

for i=1:4
total_imp{i}=double([CI{i} MD{i}]);
% total_imp{i}=array2table(total_imp{i});
% total_imp{i}.Properties.VariableNames=newNames;
SOMME{i} = sum(total_imp{i}); % di default somma sulle colonne
p{i}=SOMME{i}/d{i}; % probabilità di avere un soggetto impaired nell feature considerata
end

% plot phenotypes trend
pheno_ht=figure(4);clf;
pheno_ht1 = tiledlayout(2,2);

% estetica grafici
for i = 1:4
   nexttile;
   bar(p{i},'FaceColor',cmap4(i,:));
   set(gca, 'XTickLabel', newNames);
   ylim([0, 1]); yline(0.5, '--', 'Color', 'k', 'LineWidth', .5); % Linea tratteggiata a +0.5
   title(['Phenotype ', num2str(i)]);
   percent = round(d{i}/size(t1,1)*100);
   ylabel('Probability');
end
title(pheno_ht1,'Time point 1')

% picturewidth=35;
% hw_ratio=0.4;
% set(findall(pheno_ht,'-property','FontSize'),'FontSize',10)
% set(findall(pheno_ht,'-property','Box'),'Box', 'off')
% % set(findall(pheno_ht,'-property','Interpreter'),'Interpreter','none')
% % set(findall(pheno_ht,'-property','TickLabelInterpreter'),'TickLabelInterpreter','none')
% set(pheno_ht,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos=get(pheno_ht,'Position');
% set(pheno_ht,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[pos(3) pos(4)])
% print(pheno_ht, 'Phenot1','-dpdf','-painters','-fillpage');
% print(pheno_ht, 'Phenot1','-dpng','-painters');

%% phenotypes t2

matrix2 = zeros(size(t2,1),Nvars+1); % matrice di 272 x 11, con aggiunta colonna indici cluster
matrix2(:,1)= kT1_sd;  % applio la clusterizzazione trovata al t1
matrix2(:,2:Nvars+1)=t2; % inserisco dati t2

phenotype2 = cell(1,4);
CI2 = cell(1,4); % ogni cella sarà un vettore di 10 elementi 
MD2 = cell(1,4);
d2 = cell(1,4);

for i=1:4
    d2{i} = sum(ismember(matrix2(:,1),i)); %scopro la dimensione della cell i --> dipenda da quanta gente è nel cluster
    indixes_change2 = find(ismember(matrix2(:,1),i)); % trovo gli indici delle righe del i esimo cluster
    phenotype2{i}=zeros(d2{i},10);
    phenotype2{i}=matrix2(indixes_change2,2:11);
    impairments_threshold_under = [4 3 6 3 1 4 6 35]; impairments_threshold_over = [7 7]; % in ordine di new names
for j=1:8
    CI2{i}(:,j) = phenotype2{i}(:,j)<impairments_threshold_under(1,j);
    for l=1:2
        MD2{i}(:,l) = phenotype2{i}(:,l+8)>impairments_threshold_over(1,l);
    end
end
end

% merge 
total_imp2 = cell(1,4); %matrici di 0 e 1 per ogni soggetto del cluser impaired nella feature
SOMME2 = cell(1,4); %totale di soggetti impaired i quela feature indiepndente da chi è 
p2 = cell(1,4); %probailità per ogni colonna di avere soggetto impaired nel fenotipo 

for i=1:4
total_imp2{i}=double([CI2{i} MD2{i}]);
% total_imp{i}=array2table(total_imp{i});
% total_imp{i}.Properties.VariableNames=newNames;
SOMME2{i} = sum(total_imp2{i});
p2{i}=SOMME2{i}/d2{i};
end

% plot phenotypes trend
fig5=figure(5);clf;
pheno_ht2 = tiledlayout(2,2);


for i = 1:4
   nexttile;
   bar(p2{i},'FaceColor',cmap4(i,:));
   set(gca, 'XTickLabel', newNames);
   ylim([0, 1]); yline(0.5, '--', 'Color', 'k', 'LineWidth', .5); % Linea tratteggiata a +0.5
   title(['Phenotype ', num2str(i)]);
   percent = round(d2{i}/size(t2,1)*100);
   ylabel('Probability');
end
title(pheno_ht2,'Time point 2')

% picturewidth=35;
% hw_ratio=0.4;
% set(findall(fig5,'-property','FontSize'),'FontSize',10)
% set(findall(fig5,'-property','Box'),'Box', 'off')
% % set(findall(pheno_ht,'-property','Interpreter'),'Interpreter','none')
% % set(findall(pheno_ht,'-property','TickLabelInterpreter'),'TickLabelInterpreter','none')
% set(fig5,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos=get(fig5,'Position');
% set(fig5,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[pos(3) pos(4)])
% print(fig5, 'Phenot2','-dpdf','-painters','-fillpage');
% print(fig5, 'Phenot2','-dpng','-painters');
   %% costruzione cell phenotype sd

    matrix_sd = zeros(size(t1_sd,1),Nvars+1); matrix2_sd = zeros(size(t2_sd,1),Nvars+1);% t1_sd e t2_sd
    matrix_sd(:,1)= kT1_sd;matrix2_sd(:,1)= kT1_sd;
    matrix_sd(:,2:Nvars+1)=t1_sd;matrix2_sd(:,2:Nvars+1)=t2_sd; % qua seleziono il vettore di clustering che uso (su sd)

    phenotype_sd = cell(1,4); phenotype2_sd = cell(1,4);
    % d2 = cell(1,4);

    for i=1:4
    d{i} = sum(ismember(matrix_sd(:,1),i)); %scopro la dimensione della cell i --> dipenda da quanta gente è nel cluster
    indixes_change = find(ismember(matrix_sd(:,1),i)); % trovo gli indici delle righe del i esimo cluster
    phenotype_sd{i}=zeros(d{i},10);phenotype2_sd{i}=zeros(d{i},10);
    phenotype_sd{i}=matrix_sd(indixes_change,2:11);phenotype2_sd{i}=matrix2_sd(indixes_change,2:11);
    end

%% Line plots with T-scores

figure(6);clf;
htot = tiledlayout(4,3);
t=[0 365]; y=zeros(1,2);
for i=1:10 % fissa la feature , e per ogni feature ho un plot a sè 
    nexttile;
    hold on;
for j=1:4 %in ogni feature plot rappresento i 4 fenotipi e le loro variazioni
    y(1,1) = mean(phenotype_sd{j}(:,i),'omitnan');
    y(1,2) = mean(phenotype2_sd{j}(:,i),'omitnan');         
    hold on;
    plot(t,y,'-','LineWidth',1)
    xlim([0,364]);
    title(variations(i));
     
end
xlabel('Time');ylabel('Normalised score')
end
legend(['Phenotype', num2str(1)],['Phenotype', num2str(2)],['Phenotype', num2str(3)],['Phenotype', num2str(4)]);

%% change score distribution
cmap2=cbrewer2('Set1',10);

figure(7);clf;
tiledlayout(1,2);
ax(1)=nexttile;
hold on;
for i=1:numel(variations)
    histogram(t1_sd(:,i),'FaceColor',cmap2(i,:));
end
legend(variations,'Location','northwest');
xlabel('Baseline standardized score T1');ylabel('Proportion');

ax(2)=nexttile;
hold on;
for i=1:numel(variations)
    histogram(t2_sd(:,i),'FaceColor',cmap2(i,:));
end
legend(variations,'Location','northwest')
xlabel('Baseline standardized score T2');ylabel('Proportion');

%% costruire i box plot

% se seguiamo reference, loro hanno raggruppato i pazienti con k means
% e già ne userei un altro magari, visto che a priori non sappiamo i gruppi
% vuol dire trovare di nuovo i fenotipi
%si possono anche fare due prove, una con k=4 e l'altra libera
% e poi di questi fare i box plot per igni fenotipo

%kpca2 vettore che codifica il cluster di ogni paziente
cmap3=cbrewer2('Set3',10);

figure(8);clf;
ht1=tiledlayout(2,2);
for i=1:4
% phenotype = find(kT1_sd==i);
ax(i)=nexttile;
boxplot(phenotype_sd{i},variations,'BoxStyle','filled','ColorGroup',cmap3,'Symbol','.');
% binEdges = 25:5:75;
title(['Phenotype ', num2str(i)] );
end
title(ht1, 'Box plots - T1 sd distribution PwMS');

figure(9);clf;
ht2=tiledlayout(2,2);
for i=1:4
% phenotype = find(kT1_sd==i);
ax(i)=nexttile;
boxplot(phenotype2_sd{i},variations,'BoxStyle','filled','ColorGroup',cmap3,'Symbol','.');
% boxchart(phenotype2_sd{i},variations,'BoxWidth', 0.3)
% binEdges = 25:5:75;
title(['Phenotype ', num2str(i)] );
end
title(ht2, 'Box plots - T2 sd distribution PwMS');

% figure(12);
% ht2=tiledlayout(2,2);
% for i=1:4
% phenotype = find(kT1_sd==i); %uso la clusterizzazione dei pazienti effettuata a T1 e guardo negli stessi gruppi cosa varia di più
% ax(i)=nexttile;
% boxplot(DELTAt(phenotype,:),variations,'BoxStyle','filled','ColorGroup',cmap3,'Symbol','.');
% title(ax(i),i);
% end
% title(ht2, 'Box plots - Delta distribution PwMS');



%% Demographic information 

opts_age = spreadsheetImportOptions('Sheet','AGE','VariableNamesRange','A1','DataRange','A2:A284','VariableTypes',repmat({'double'},1,1));
T1_age = readtable('DEMOGRAPHIC1.xlsx',opts_age);

T1_age.AGE = categorical(T1_age.AGE);
t1_age = T1_age(~pat_thresh,:);

opts_edss = spreadsheetImportOptions('Sheet','EDSS','VariableNamesRange','A1','DataRange','A2:A284','VariableTypes',repmat({'double'},1,1));
T1_edss = readtable('DEMOGRAPHIC1.xlsx',opts_edss);

T1_edss.EDSS = categorical(T1_edss.EDSS);
t1_edss = T1_edss(~pat_thresh,:);

opts_edu = spreadsheetImportOptions('Sheet','EDU','VariableNamesRange','A1','DataRange','A2:A284','VariableTypes',repmat({'double'},1,1));
T1_edu = readtable('DEMOGRAPHIC1.xlsx',opts_edu);

T1_edu.EDU = categorical(T1_edu.EDU);
t1_edu = T1_edu(~pat_thresh,:);

opts_last = spreadsheetImportOptions('Sheet','LAST','VariableNamesRange','A1','DataRange','A2:A284','VariableTypes',repmat({'double'},1,1));
T1_last = readtable('DEMOGRAPHIC1.xlsx',opts_last);

T1_last.LAST = categorical(T1_last.LAST);
t1_last = T1_last(~pat_thresh,:);

opts_cour = spreadsheetImportOptions('Sheet','COURSE','VariableNamesRange','A1','DataRange','A2:A284');
T1_cour = readtable('DEMOGRAPHIC1.xlsx',opts_cour);

T1_cour.COUR = categorical(T1_cour.COUR);
t1_cour = T1_cour(~pat_thresh,:);

opts_gend = spreadsheetImportOptions('Sheet','GEND','VariableNamesRange','A1','DataRange','A2:A284');
T1_gend = readtable('DEMOGRAPHIC1.xlsx',opts_gend);

T1_gend.GEND = categorical(T1_gend.GEND);
t1_gend = T1_gend(~pat_thresh,:);


table_tot_cat=[t1_age.AGE t1_edss.EDSS t1_edu.EDU t1_last.LAST];
table_tot = str2double(cellstr(table_tot_cat));

cat_table = [t1_gend.GEND t1_cour.COUR];

%% histo 
demo2 = cell(1,4); %4 celle per 4 fenotipi

Cat2 = ["GENDER" , "DISEASE COURSE"];
dim_demo2 = cell(1,4);
% means2 = zeros(4,2); %ogni riga sta per le medie di un fenotipo, le colonne stanno per le categorie
% err2 = zeros(4,2);

for i=1:4
    dim_demo2{i}=sum(ismember(kT1_sd(:,1),i));  %scopro la dimensione della cell i --> dipenda da quanta gente è nel cluster
    indixes_change_demo = find(ismember(kT1_sd(:,1),i));
    
    for j=1:2
        demo2{i}(:,j)=cat_table(indixes_change_demo,j);
       
    end
end
%% withney solo per sesso e corso
MF = zeros(2,4);
PP_SP_RR = zeros(3,4);  j=1;

for i=1:4
   
        MF(1,i)=sum(ismember(demo2{i}(:,j),'M'));MF(2,i)=sum(ismember(demo2{i}(:,j),'F'));
        PP_SP_RR(1,i)=sum(ismember(demo2{i}(:,j+1),'PPMS'));PP_SP_RR(2,i)=sum(ismember(demo2{i}(:,j+1),'SPMS'));PP_SP_RR(3,i)=sum(ismember(demo2{i}(:,j+1),'RRMS'));
end

P_mann2 = cell(1,2); 
% P_tab2 = cell(1,2); Pttest = cell(1,2);
P_mann2{1} = zeros(4,4);P_tab2 = cell(1,1);
P_mann2{2} = zeros(4,4);
% Pheno_names = ["Phenotype 1","Phenotype 2","Phenotype 3","Phenotype 4" ];
for i=1:4
    for j=1:4
            % [~, Pttest{k}(i,j)] = ttest2(MF(:,i), MF(:,k));
            P_mann2{1}(i,j) = ranksum(MF(:,i), MF(:,j)); %Mann-Witnhey U test non parametric
            P_tab2 = array2table(P_mann2{1});
            P_mann2{2}(i,j)= ranksum(PP_SP_RR(:,i), PP_SP_RR(:,j));
            % P_tab2{k}= renamevars(Pttest{k},1:4,Pheno_names); 
    end
end

%% histo 
demo = cell(1,4); %4 celle per 4 fenotipi
DEMO = cell(1,4);
Cat = ["AGE","EDSS", "YEARS OF EDUCATION" ,"DISEASE DURATION"];
dim_demo = cell(1,4);
means = zeros(4,4); %ogni riga sta per le medie di un fenotipo, le colonne stanno per le categorie
err = zeros(4,4);

for i=1:4
    dim_demo{i}=sum(ismember(kT1_sd(:,1),i));  %scopro la dimensione della cell i --> dipenda da quanta gente è nel cluster
    indixes_change_demo = find(ismember(kT1_sd(:,1),i));
    demo{i}=zeros(dim_demo{i},4);
    for j=1:4
        demo{i}(:,j)=str2double(cellstr(table_tot_cat(indixes_change_demo,j)));
        [err(i,j),means(i,j)]=std(demo{i}(:,j));
    end
end
Means = array2table(means);Err = array2table(err);
Means.Properties.VariableNames=Cat;Err.Properties.VariableNames=Cat;
        
for i=1:4 % DEMOGRAPHIC 
    figure(9+i);clf;
    ht(i)=tiledlayout(2,2);
        for j=1:4
        ax(j)=nexttile;
        histogram(demo{j}(:,i),'FaceColor',cmap4(j,:))
        ylabel(dim_demo{j})
        title(['Phenotype ', num2str(j)])
        end
    title(ht(i),Cat{i});
end

%% verifica normalità e omoschedasticità
% ipotesi nulla che i dati vengano da una distribuzione normale
% test di Jarque Bera
p_norm = zeros(4,4);h_jb = zeros(4,4);
for j=1:4 %fenotipi
    for i=1:4 %demografiche
        [h_jb(j,i),p_norm(j,i)] = jbtest(demo{j}(:,i));
    end
end
P_norm = array2table(p_norm);
P_norm.Properties.VariableNames = Cat;
%% KRUSKAL-WALLIS (ANOVA NON PARAM VERSION)
    
   % kruskal wallis non parametrico per valutare se i dati demografici 
   % vengano tutti da una stessa distribuzione (H0) opppure no

   p_krusk = zeros(1,4); p_anova = zeros(1,4);
  
   for i = 1:4 %fenotipo
        for j=1:4 %categorie
        [p_krusk(1,i), ~, ~] = anova1('kruskalwallis',str2double(cellstr(table_tot_cat(:,i))), kT1_sd, 'off');
        [p_anova(1,i), ~, ~] = anova1(str2double(cellstr(table_tot_cat(:,i))), kT1_sd, 'off');
        end
   end

%% MANN-WITHNEY U TEST (T TEST NON PARAM VERSION)
        
% Is there an effect in the population? t test per ogni combinazione di gruppi
% there is a significant difference among the mean of the demographic
% features that says that these phenotypes groupe together 
% people with different demographic features?

P_mann = cell(1,4); P_tab = cell(1,4); P_ttest=cell(1,4);
Pheno_names = ["Phenotype 1","Phenotype 2","Phenotype 3","Phenotype 4" ];
for k=1:4 % indice che corre fra age, edss, edu, last
    P_mann{k}=zeros(4,4);
    P_ttest{k}=zeros(4,4);
    for j = 1:4  
        for i =1:4
            [~, P_ttest{k}(i,j)] = ttest2(demo{i}(:,k), demo{j}(:,k));
            P_mann{k}(i,j) = ranksum(demo{i}(:,k), demo{j}(:,k)); %Mann-Witnhey U test non parametric
            P_tab{k} = array2table(P_mann{k});
            P_tab{k}= renamevars(P_tab{k},1:4,Pheno_names); 
            
        end
    end
end

% same issue, but I confer the mean of each demo feature of 1 phenotype
% with all the population, not only another single cluster. 
Ptot_mann = zeros(4,4);Ptot_ttest = zeros(4,4);
   
for k=1:4 % age, edss, edu, last
    for i =1:4 %fenotipo
        Ptot_mann(i,k) = ranksum(demo{i}(:,k), table_tot(:,k));
        Ptot_ttest(i,k) = ttest2(demo{i}(:,k), table_tot(:,k));
    end
end
Ptot_table = array2table(Ptot_mann);
Ptot_table.Properties.VariableNames = ["Phenotype 1 vs. ALL","Phenotype 2 vs. ALL","Phenotype 3 vs. ALL","Phenotype 4 vs. ALL" ];
Cat_table = array2table(Cat');Cat_table.Var1 = categorical(Cat_table.Var1);

Ptot_table = addvars(Ptot_table, Cat_table.Var1);Ptot_table = renamevars(Ptot_table,5,"Demography");
%% P-value image demo infomation

num = [1 2 3 4];

figure(14);clf;
hp1 = tiledlayout(2,2);

for i = 1:4
    ax(i)=nexttile;
    imagesc(P_mann{i});axis equal;
    ax(i)=gca;
    ax(i).XTick=1:4;ax(i).XTickLabel=(num);
    ax(i).YTick=1:4;ax(i).YTickLabel=(num);
    ax(i).TickLabelInterpreter='none';
    title(Cat(i),'FontWeight','normal');
    colormap('parula');
    colorbar;
    % ax(i).CLim=[0,1];
end
title(hp1,'Mann Withney')
% quali differenze sono significative?

figure(15);clf;
hp2 = tiledlayout(2,2);

for i = 1:4
    ax(i)=nexttile;
    imagesc(P_ttest{i});axis equal;
    ax(i)=gca;
    ax(i).XTick=1:4;ax(i).XTickLabel=(num);
    ax(i).YTick=1:4;ax(i).YTickLabel=(num);
    ax(i).TickLabelInterpreter='none';
    title(Cat(i),'FontWeight','normal');
    colormap('parula');
    colorbar;
    % ax(i).CLim=[0,1];
end
title(hp2,'T test')
% quali differenze sono significative?

P_th_mann = cell(1,4);
for k=1:4
    P_th_mann{k} = P_mann{k}<0.05;
end
P_th_ttest = cell(1,4);
for k=1:4
    P_th_ttest{k} = P_ttest{k}<0.05;
end


figure(16);clf;
h_th1 = tiledlayout(2,2);
for i=1:4
    ax(i)=nexttile;
    imagesc(P_th_mann{i});axis image;
    ax(i)=gca;
    ax(i).XTick=1:4;ax(i).XTickLabel=(num);
    ax(i).YTick=1:4;ax(i).YTickLabel=(num);
    ax(i).TickLabelInterpreter='none';
    title(Cat(i),'FontWeight','normal');
    colormap('bone')
end
title(h_th1,'Mann withney test (0.05)')

figure(17);clf;
h_th2 = tiledlayout(2,2);
for i=1:4
    ax(i)=nexttile;
    imagesc(P_th_ttest{i});axis image;
    ax(i)=gca;
    ax(i).XTick=1:4;ax(i).XTickLabel=(num);
    ax(i).YTick=1:4;ax(i).YTickLabel=(num);
    ax(i).TickLabelInterpreter='none';
    title(Cat(i),'FontWeight','normal');
    colormap('bone')
end
title(h_th2,'T test (0.05)')
% 14-17
%% versus tutta la popolazione:
figure(18);clf;
    nexttile;
    imagesc(Ptot_ttest);axis equal;
    ax=gca;
    ax.XTick=1:4;ax.XTickLabel=(Cat);
    ax.YTick=1:4;ax.YTickLabel=(num);
    ax.TickLabelInterpreter='none';
    title('Phenotype i vs. ALL - T test','FontWeight','normal');
    colormap('parula');
    colorbar;

    figure(19);clf;
    nexttile;
    imagesc(Ptot_mann);axis equal;
    ax=gca;
    ax.XTick=1:4;ax.XTickLabel=(Cat);
    ax.YTick=1:4;ax.YTickLabel=(num);
    ax.TickLabelInterpreter='none';
    title('Phenotype i vs. ALL - Mann Withney test','FontWeight','normal');
    colormap('parula');
    colorbar;
    
% quali differenze sono significative?

Ptot_th_ttest = Ptot_ttest<0.05;
figure(20);clf;

    nexttile;
    imagesc(Ptot_th_ttest);axis equal;
    ax=gca;
    ax.XTick=1:4;ax.XTickLabel=(Cat);
    ax.YTick=1:4;ax.YTickLabel=(num);
    ax.TickLabelInterpreter='none';
    title('Phenotype i vs. ALL - T test (0.05)','FontWeight','normal');
    colormap('bone')

    Ptot_th_mann = Ptot_mann<0.05;
    figure(21);clf;

    nexttile;
    imagesc(Ptot_th_mann);axis equal;
    ax=gca;
    ax.XTick=1:4;ax.XTickLabel=(Cat);
    ax.YTick=1:4;ax.YTickLabel=(num);
    ax.TickLabelInterpreter='none';
    title('Phenotype i vs. ALL - Mann Withney (0.05)','FontWeight','normal');
    colormap('bone')
%% T test allo stesso tempo confer medie delle stesse features differenze
Pmean = zeros(4,10);Pmean_wilc = zeros(4,10);
% confer le medie delle feature nello stesso fenotipo t1 e t2 --> campioni
% appaiati (t test appaiato è ok)
% normalità?
for k=1:10 %corre lungo le features
    for i = 1:4  
        [~, Pmean(i,k)] = ttest2(phenotype_sd{i}(:,k), phenotype2_sd{i}(:,k));
        Pmean_wilc(i,k) = signrank(phenotype_sd{i}(:,k), phenotype2_sd{i}(:,k)); % signed rank wilcoxon test for paired samples
    end
end

Pmean_th = Pmean<0.05;Pmean_wilc_th = Pmean_wilc<0.05;

% Wilcoxon valuta se la differenza tra le due distribuzioni possa venire da
% una ditribuzione con MEDIANA 0.
% il t test o quello di mann witney valutano se le due popolazioni possano
% provenire da due distribuzioni uguali o diverse. e valutano le MEDIE.
%% P value imaging

figure(22);clf;
imagesc(Pmean);
% axis equal;
ax=gca;
ax.XTick=1:10;ax.XTickLabel=(variations);
ax.YTick=1:4;ax.YTickLabel=(num);
ax.TickLabelInterpreter='none';
title('T1 vs T2 for each feature-T test','FontWeight','normal');
colormap('parula')
colorbar

figure(23);clf;
imagesc(Pmean_th);
% axis equal;
ax=gca;
ax.XTick=1:10;ax.XTickLabel=(variations);
ax.YTick=1:4;ax.YTickLabel=(num);
ax.TickLabelInterpreter='none';
title('T1 vs T2 for each feature - T test (0.05)','FontWeight','normal');
colormap('bone')

figure(24);clf;
imagesc(Pmean_wilc_th);
% axis equal;
ax=gca;
ax.XTick=1:10;ax.XTickLabel=(variations);
ax.YTick=1:4;ax.YTickLabel=(num);
ax.TickLabelInterpreter='none';
title('T1 vs T2 for each feature - Wilcoxon signed rank (0.05)','FontWeight','normal');
colormap('bone')

        
%% Cohen's D coefficient

Dcoeff = cell(1,4); Dcoeff_correct = cell(1,4); %matrice di coefficienti di Cohen

% it describes the number of standard deviations that separates two groups,
% people in a group are confered at t1 and t2 for each feature evaluated 

for i=1:10 %  per ogni feature
    for j=1:4 % in ogni fenotipo
        Dcoeff{j}(i,:)= meanEffectSize(phenotype_sd{j}(:,i),phenotype2_sd{j}(:,i),"Effect","cohen","Paired",true);
    end
end

for i = 1:4
    if d{i}<50
        Dcoeff_correct{i}(:,1) = Dcoeff{i}.Effect*sqrt((d{i}-2)/d{i})*((d{i}-3)/(d{i}-2.25));
    else
        Dcoeff_correct{i}(:,1) = Dcoeff{i}.Effect;
    end
end


fig=figure(25);clf;
htot2=tiledlayout(2,2);
for i=1:4
    nexttile;
    Dcoeff_group = (Dcoeff_correct{i}(:, 1));
    for j=1:10
        b = bar(j, Dcoeff_group(j), 'FaceColor', cmap4(j, :)); % Usa j per selezionare il colore appropriato dal cmap4
        hold on;
    end
    yline(0.5, '--', 'Color', 'k', 'LineWidth', .5); % Linea tratteggiata a +0.5
    yline(-0.5, '--', 'Color', 'k', 'LineWidth', .5);
    hold off;
    % title(['Phenotype ', num2str(i)]); % Titolo del subplot basato sul numero del gruppo
end
xlabel(htot2,'Features');
ylabel(htot2,'Cohen''s D');
lg=legend(variations);
lg.Location="bestoutside";

% picturewidth=35;
% hw_ratio=0.4;
% set(findall(fig,'-property','FontSize'),'FontSize',10)
% set(findall(fig,'-property','Box'),'Box', 'off')
% % set(findall(pheno_ht,'-property','Interpreter'),'Interpreter','none')
% % set(findall(pheno_ht,'-property','TickLabelInterpreter'),'TickLabelInterpreter','none')
% set(fig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos=get(fig,'Position');
% set(fig,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[pos(3) pos(4)])
% print(fig, 'cohen_leg','-dpdf','-painters','-fillpage');
% print(fig, 'cohen_leg','-dpng','-painters');

%% section disease course

d_course = cell(1,4); %4 celle per 4 fenotipi


for i=1:4 %fenotipo
    % dim_demo{i}=sum(ismember(kT1_sd(:,1),i));  %scopro la dimensione della cell i --> dipenda da quanta gente è nel cluster
    d_course{i}=zeros(dim_demo{i},4);
    indixes_change_demo = find(ismember(kT1_sd(:,1),i));
    d_course{i}=t1_cour.COUR(indixes_change_demo);
    
end


figure(26);clf;
ht_cour=tiledlayout(2,2);
   for j=1:4
       ax(j)=nexttile;
       histogram(d_course{j},'FaceColor',cmap4(j,:))
       ylabel(dim_demo{j})
       title(['Phenotype ', num2str(j)])
   end
title(ht_cour,"DISEASE COURSE");










              
              
  

 







   
  
  









    
