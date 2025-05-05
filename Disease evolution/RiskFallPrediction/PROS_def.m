%
%
%
% analisi preliminare tabella PROMOPROMS
%
%
clearvars;
set(0,'DefaultFigureWindowStyle','docked');
load('data_promo2ID.mat');
%% sistemazione variabili
T = Torig;
clear Torig;
eqCOUR=categorical(T.COUR, ...
	["Benigna" "Ricadute - remissioni" "Secondaria progressiva" "Progressiva con ricadute" "Primaria progressiva"], ...
	["RR" "RR" "SP" "SP" "PP"],'ordinal',true);
T=addvars(T,eqCOUR,'after','COUR');
T(isundefined(T.eqCOUR),:)=[];
T.Paziente=string(T.Paziente);
T.BMI=T.WEIGHT./((T.HEIGHT/100).^2);
MAXFIM=126;
MAXEDSS=10;
T.FIM_TOT(T.FIM_TOT>MAXFIM)=MAXFIM;
missingedss=isnan(T.EDSS006); % i nan nelle singole edss danno 0 di edss totale
T(missingedss,:)=[];

%% 
% for j=1:numel(up2)
% 	w=T.Paziente==up2(j);
%     wf = find(w);
%     Tw = T(w,:);
%     FUP(j,1)= Tw.fup(end);
% end
% [MEAN,SD]=std(FUP);
%% costruzione del dataset per la RF 
tstar=1;                                                                    % definiamo 1y di fup come baseline. con l'ipotesi che ci siano ~ 3 visite / anno per tutti
nvisit=tstar*3;
theorvisittime=years(linspace(0,tstar,nvisit));
Ddependent={'Delta_EDSSTOT' 'Delta_FIM_TOT' 'Delta_LIFETOT' 'Delta_MFISTOT'}; % differenti scale da predire
Bdependent={'EDSSTOT' 'FIM_TOT' 'MFISTOT' 'LIFETOT'};
independent0={'EDSSTOT' 'FIM_TOT' 'BMI' 'ABILHTOT'  'OAB_QTOT' 'HADSTOT' 'LIFETOT' 'MFISTOT' 'SDMTTOT'}; % predittori finali selezionati

clabel=["base" "slope"];
fupslabelfun=@(x,n) sprintf('%s_%s',x,clabel(n));
independent=cell(2,1);
for j=1:2 % baseline e trend
	independent{j}=cellfun(@(x) fupslabelfun(x,j),independent0,'UniformOutput',false);
end

problem={'deltaT'};proof = {'dummy'};
aux={'Paziente'};
external = {'eqCOUR' 'OtherDisease' 'RELAPS'};
comorb = {'metabolic' 'cardiovascular' 'psychic' 'onco-neuro' 'other'};
demo = {'AGE' 'GEND' 'Duration_from_DD'};
%% Ttrain con 1 y fup
Ttrain=table();
up=unique(T.Paziente);
for j=1:numel(up)
	w=T.Paziente==up(j);
	Tpat=T(w,:);
	[~,pos]=sort(Tpat.Date,'ascend');
	Tpat=Tpat(pos,:); % matrice del paziente ordinata temporalmente
	if range(Tpat.Date)>years(tstar) % se il range di fup del paziente è maggiore di 1 anno
		visitbeforetstar=(Tpat.Date-min(Tpat.Date))<=years(tstar); % rendo le viiste in formato anno giusto       
		visittime=duration(Tpat.Date-min(Tpat.Date),'format','y');  % prendo le visite prima dell'anno		
		% var=zeros(numel(independent0),numel(theorvisittime));% costruisco una matrice con 4 colonne (numero di viSiste) e righe tante quante le variabili predittori
		var=zeros(numel(independent0),2);% costruisco una matrice con 2 colonne (baseline e trend) e righe tante quante le variabili predittori
		for v=1:numel(independent0) % loppo sulle variabili predittrici
			tstarvisittime=visittime(visitbeforetstar);% prendo le date delle visite pre anno
			visit=Tpat.(independent0{v})(visitbeforetstar); % di ogni variabile prendo i valori pre anno
			nan=isnan(visit);% vedo i nan della variabile
			switch numel(visit) % per ogni numero di visite effettive che ha il paziente
				case 1 % se ha una visita sola
					% var(v,:)=repmat(visit,[1,nvisit]);% riempio i 3 valori con il primo valore
					var(v,:)=[NaN,NaN];% BASELINE e DERIVATA NULLA
				case 2 % se ha due visite 
					switch nnz(nan)  % numero di missin values
						case 0 % NESSUN MISSING VALUE
							% var(v,:)=interp1(tstarvisittime,visit,theorvisittime,'nearest','extrap');
							lm=polyfit(years(tstarvisittime),visit,1);
							var(v,:)=[visit(1),lm(1)];
						otherwise
							var(v,:)=[NaN,NaN];
                    end % se sono tutti nan lasciamo i nan e poi il soggetto se ne andrà 
				otherwise
					if nnz(~nan)>=2
					    p=polyfit(years(tstarvisittime(~nan)),visit(~nan),1);
					    f=find(visit(~nan),1,'first');
					    var(v,:)=[visit(f),p(1)];
					else
						var(v,:)=[NaN,NaN];
                    end
			 end
		end
		deltaT=visittime(~visitbeforetstar)-years(tstar);  % GIA' COSTRUISCE MATRICE A PARTIRE DA TSTAR IN POI
		Vdelta=zeros(nnz(~visitbeforetstar),numel(Bdependent));
		for v=1:numel(Bdependent)
			%Vdelta(:,v)=table2array(Tpat(~visitbeforetstar,Bdependent(v))-Tpat(find(visitbeforetstar,1,'last'),Bdependent(v)));
            Vdelta(:,v)=table2array( Tpat(~visitbeforetstar,Bdependent(v))-mean(Tpat(visitbeforetstar,Bdependent(v))) );
		end
		Cdelta=cumsum(Vdelta)./years(deltaT);  % il delta edss e delta fim sono medie mobili normalizzate per il tempo (variazione cumulativa)
		varbase=repmat(var(:)',[nnz(~visitbeforetstar),1]);
		Tpatbase=table(deltaT);
		Tpatbase=[Tpatbase,array2table([varbase,Cdelta],'variablenames',[independent{:},Ddependent])];
		Tpatbase=[Tpat(~visitbeforetstar,[aux,demo,external]),Tpatbase];
		Ttrain=[Ttrain;Tpatbase];
	end
end

%% pulizia dei missing & outliers
FOCUS = "Delta_EDSSTOT"; % scelta del task
Ttrain(Ttrain.RELAPS>0,:)=[]; % eliminato chi ha vauto ricadute nei 4 mesi precedenti
missinglines=any(ismissing(Ttrain(:,[demo,independent{:},FOCUS])),2); % missin della table finale
Ttrain(missinglines,:)=[];

% eliminiamo outliers
out_id = Ttrain.(FOCUS) < quantile(Ttrain.(FOCUS),0.05) | Ttrain.(FOCUS)> quantile(Ttrain.(FOCUS),0.99);
Ttrain(out_id,:)=[];

%% dummy + comorbidities

% COMORBIDITY
cat_morbius=categorize_comorbidities(Ttrain.OtherDisease);
comorbidities=sum(cat_morbius,2);
Ttrain.comorbidities=categorical(comorbidities.sum,0:4,{'none' 'low' 'medium' 'high' 'severe'},'ordinal',true);
% DUMMY VARIABLE
Ttrain.dummy= rand(height(Ttrain),1);  
%% Definizione classi 

risk_class=cell(size(Ttrain,1),1);

id_lowbase = Ttrain.EDSSTOT_base==0; % 11
id_medbase = Ttrain.EDSSTOT_base>0 & Ttrain.EDSSTOT_base<5.5;%[1-5],2251
id_highbase = Ttrain.EDSSTOT_base>5; %  2047
base_disability=[id_lowbase,id_medbase,id_highbase];

cat_label=cell(height(Ttrain),1); % gruppo di appartenenza cat_label(id_mediumbase)="medium";
labels = {'low' 'medium' 'high'};
for i=1:3
    cat_label(base_disability(:,i))=labels(i);
end
cat_label = categorical(cat_label,string(labels),'ordinal',true);

switch FOCUS
    case "Delta_EDSSTOT"
        names = {'safe' 'moderate' 'at risk'}; 
        L1=[0.1 1.5];L2=[0.1 1];L3=[0.1 0.5]; % discretize prende il limite dx
        risk_class(id_lowbase,1)=(discretize(Ttrain.(FOCUS)(id_lowbase),[-Inf L1 Inf],names)); % basso edss
        risk_class(id_medbase,1)=discretize(Ttrain.(FOCUS)(id_medbase),[-Inf L2 Inf],names); % medio edss
        risk_class(id_highbase,1)=discretize(Ttrain.(FOCUS)(id_highbase),[-Inf L3 Inf],names); % alto edss        
    case "Delta_FIM_TOT"
        names = {'at risk' 'moderate' 'safe'}; 
        quantili = quantile(Ttrain.(FOCUS)(Ttrain.(FOCUS)<0), .5);
        L = [quantili 0];        
        risk_class=discretize(Ttrain.(FOCUS),[-Inf L Inf],names); 
    case "Delta_LIFETOT"
        names = {'at risk' 'moderate' 'safe'}; 
        quantili = quantile(Ttrain.(FOCUS)(Ttrain.(FOCUS)<0), .5);
        L = [quantili 0]; % DELTA DI 10
        risk_class= discretize(Ttrain.(FOCUS),[-Inf L Inf],names); 
    case "Delta_MFISTOT"
        names = {'safe' 'moderate' 'at risk'}; 
        quantili = quantile(Ttrain.(FOCUS)(Ttrain.(FOCUS)>0), .5);
        L = [0.1 quantili];
        risk_class= discretize(Ttrain.(FOCUS),[-Inf L Inf],names); 
   end
% risk_class=categorical(risk_class,1:3,names,'Ordinal',true);
risk_class=categorical(risk_class,names,'ordinal',true);

%% histogram delle parti prese 
f=figure(1);clf;
switch FOCUS
    case "Delta_EDSSTOT"
        tiledlayout(2,1);
        ax(1) = nexttile;
        histogram(Ttrain.(FOCUS)(id_medbase),100)
        title('Classes division - patients with EDSS baseline [1-5]')
        xlabel(sprintf('%s from mean baseline value',FOCUS),'Interpreter','none')
        h=gca;
        hold on;
        x_lim = h.XLim;
        y2_lim = h.YLim;

        % Definisci le coordinate dei rettangoli
        x1 = [x_lim(1), L1(1), L1(1), x_lim(1)];
        x2 = [L1(1),  L1(2), L1(2),  L1(1)];
        x3 = [L1(2), x_lim(2), x_lim(2), L1(2)];
    
        % Colori per i rettangoli
        colors = [1, 0, 0; 0, 1, 0; 0, 0, 1]; % rosso, verde, blu
    
        % Disegna i rettangoli
        fill(x1, [y2_lim(1), y2_lim(1), y2_lim(2), y2_lim(2)], colors(1, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        fill(x2, [y2_lim(1), y2_lim(1), y2_lim(2), y2_lim(2)], colors(2, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        fill(x3, [y2_lim(1), y2_lim(1), y2_lim(2), y2_lim(2)], colors(3, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        
        text(mean(x1(1:2)), y2_lim(2)*0.95, 'Safe', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12, 'Color', 'k');
        text(mean(x2(1:2)), y2_lim(2)*0.95, 'Moderate', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12,  'Color', 'k');
        text(mean(x3(1:2)), y2_lim(2)*0.95, 'At Risk', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12,  'Color', 'k');
        
        hold off;

        ax(2) = nexttile;
        histogram(Ttrain.(FOCUS)(id_highbase),100)
        title(ax(2),'Classes division - sample with EDSS baseline => 5.5')
        xlabel(sprintf('%s from mean baseline value',FOCUS),'Interpreter','none')
        h1=gca;
        hold on;
        x_lim = h1.XLim;
        y_lim = h1.YLim;

        % Definisci le coordinate dei rettangoli
        x1 = [x_lim(1), L2(1), L2(1), x_lim(1)];
        x2 = [L2(1),  L2(2), L2(2),  L2(1)];
        x3 = [L2(2), x_lim(2), x_lim(2), L2(2)];
        % Colori per i rettangoli
        colors = [1, 0, 0; 0, 1, 0; 0, 0, 1]; % rosso, verde, blu

        % Disegna i rettangoli
        fill(x1, [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], colors(1, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        fill(x2, [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], colors(2, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        fill(x3, [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], colors(3, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        text(mean(x1(1:2)), y_lim(2)*0.95, 'Safe', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12,  'Color', 'k');
        text(mean(x2(1:2)), y_lim(2)*0.95, 'Moderate', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12,  'Color', 'k');
        text(mean(x3(1:2)), y_lim(2)*0.95, 'At Risk', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12,  'Color', 'k');

        hold off;
    otherwise
        histogram(Ttrain.(FOCUS),100)
        title('Classes division')
        xlabel(sprintf('%s from mean baseline value',FOCUS),'Interpreter','none')
        h=gca;
        hold on;
        x_lim = h.XLim;
        y3_lim = h.YLim;

    % Definisci le coordinate dei rettangoli
    x1 = [x_lim(1), L(1), L(1), x_lim(1)];
    x2 = [ L(1),  L(2), L(2),  L(1)];
    x3 = [L(2), x_lim(2), x_lim(2), L(2)];
    
    % Colori per i rettangoli
    colors = [1, 0, 0; 0, 1, 0; 0, 0, 1]; % rosso, verde, blu
    
    % Disegna i rettangoli
    fill(x1, [y3_lim(1), y3_lim(1), y3_lim(2), y3_lim(2)], colors(1, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    fill(x2, [y3_lim(1), y3_lim(1), y3_lim(2), y3_lim(2)], colors(2, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    fill(x3, [y3_lim(1), y3_lim(1), y3_lim(2), y3_lim(2)], colors(3, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    text(mean(x1(1:2)), y3_lim(2)*0.95, 'Safe', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12, 'Color', 'k');
    text(mean(x2(1:2)), y3_lim(2)*0.95, 'Moderate', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12,  'Color', 'k');
    text(mean(x3(1:2)), y3_lim(2)*0.95, 'At Risk', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12,  'Color', 'k');
    
    hold off;

end
% picturewidth=45;
% hw_ratio=0.5;
% set(findall(f,'-property','FontSize'),'FontSize',12)
% set(findall(f,'-property','Box'),'Box', 'off')
% % set(findall(ht,'-property','Interpreter'),'Interpreter','none')
% % set(findall(ht,'-property','TickLabelInterpreter'),'TickLabelInterpreter','none')
% set(f,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos=get(f,'Position');
% set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[pos(3) pos(4)])
% print(f, 'sample','-dpdf','-painters','-fillpage');
% print(f, 'sample','-dpng','-painters');
%% creo la matrice finale ad hoc per l'algortimo utilizzato



%% SCELTA del gruppo di EDSS baseline

TtrainRF = addvars(Ttrain, cat_label, risk_class);
TtrainRF.deltaT=years(TtrainRF.deltaT);

% TtrainNN = 
base_disability=[base_disability,true(height(TtrainRF),1)];
% subcategories = {all};
Class =3; 
% choice = base_disability(:,Class);

TtrainRF=TtrainRF(base_disability(:,Class),:);
Tsub=Ttrain(base_disability(:,Class),:);
Realclasses = risk_class(base_disability(:,Class));

vars=[problem demo independent{:} 'comorbidities' 'eqCOUR']; 
% vars=[problem demo independent{:} 'comorbidities' 'cat_label' 'eqCOUR']; % demo

%% grande k-fold sui pazienti

pat = unique(TtrainRF.Paziente);
% cv= cvpartition(numel(pat),'Kfold',20);
% cv=cvpartition(numel(pat),'Kfold',ceil(numel(pat)/(numel(pat)^(1/4))));
cv=cvpartition(numel(pat),'Leaveout');
predictions=Realclasses; 
oobvi=zeros(cv.NumTestSets,numel(vars));

% fare la prova di dare le label finali a caso e veder come predice 
scores=zeros(height(TtrainRF),numel(names));
[sc,pred]=deal(cell(cv.NumTestSets,1));
% misclass_mat = [0 1 1; 1.2 0 1; 1.4 1.2 0];                                       % missclassification cost matrix, moderate risk safe
leaf = 1;
ma = {};
ht=tic;
fprintf('starting CV ...\n');
tr = templateTree('MinLeafSize',leaf,'surrogate','on',"PredictorSelection","interaction-curvature","SplitCriterion","gdi");
parfor j=1:cv.NumTestSets
	wtrain = ismember(TtrainRF.Paziente,pat(training(cv,j)));
	wtest = ismember(TtrainRF.Paziente,pat(test(cv,j)));
    % wtrain=training(cv,j);
    % wtest = test(cv,j);
	fitmodel = fitcensemble(TtrainRF(wtrain,:),'risk_class',...
		'prior','uniform','method','RUSBoost',...
		'predictornames',vars,'learners',tr,'ClassNames',names,...
		'NumLearningCycles',150); %,'cost',misclass_mat, ,'LearnRate',0.8
    [pred{j},sc{j}]=predict(fitmodel,TtrainRF(wtest,:));                                                                                                                                                    
	[oobvi(j,:),ma{j}]=predictorImportance(fitmodel);
end

for j=1:cv.NumTestSets
	wtest = ismember(TtrainRF.Paziente,pat(test(cv,j)));
    % wtest = test(cv,j);
	predictions(wtest)=pred{j};
	scores(wtest,:)=sc{j};
end
dt=seconds(toc(ht));
fprintf('done in %s minutes\n',minutes(dt)); % pp 577 line, rr 1685
%% ADD ID
IDrf = NaN(size(TtrainRF,1),1);
TtrainRF = addvars(TtrainRF, IDrf, 'Before','Paziente');

for i=1:numel(pat)
    w = TtrainRF.Paziente==pat(i);
    Tw = TtrainRF(w,:);
    ID = T.IDmri(T.Paziente==unique(Tw.Paziente));
    if nnz(~isnan(ID))>0
       TtrainRF.IDrf(w) = ID(1);
    end
end
%% plot Variable Importance - boxchart

figure(2);clf;
ht=tiledlayout("flow");
[~,pos]=sort(mean(oobvi,1),'descend');
ax=nexttile(ht);
boxchart(oobvi(:,pos),'markerstyle','.');
ylabel("Predictor Importance Estimates","FontWeight","bold")
xlabel("Predictors","FontWeight","bold")
ax.XTick=categorical(1:numel(vars));
ax.XTickLabel=vars(pos);
ax.XTickLabelRotation = 45;
ax.TickLabelInterpreter = "none";
grid on;
%% plot della PDF
% plotto solo la dummy
id_dummy = find(string(vars) == "dummy");
% f=figure(2);clf;
hf=histfit(oobvi(:,id_dummy),15);                                           % 15 numero di bin arbitrario
pd=fitdist(oobvi(:,id_dummy),'normal');
cut=pd.mu+3*pd.sigma;
%% plot Variable Importance - histogram of average importance with bar errors
SD_importance= zeros(size(oobvi,2),1);mean_importance= zeros(size(oobvi,2),1);
for i=1:size(oobvi,2)
	[SD_importance(i),mean_importance(i)] = std(oobvi(:,i));
end
[POS,IDX]=sort(mean_importance,'descend');

f=figure(3);clf;
x= 1:size(oobvi,2);
bar(x,POS);
hold on
errorbar(x,POS,SD_importance(IDX),'k','linestyle','none','LineWidth',3)';
% hold on
% yline(cut,'-',{'Cut (3 \sigma)',sprintf('x = %1.1f',cut)}) 
% hold off
h = gca;
h.XTick=1:numel(vars);
h.XTickLabel = vars(IDX);
h.XTickLabelRotation = 45;
h.TickLabelInterpreter = "none";
ylabel("Predictor Importance Estimates","FontWeight","bold")
xlabel("Predictors","FontWeight","bold")
% 
% picturewidth=45;
% hw_ratio=0.5;
% set(findall(f,'-property','FontSize'),'FontSize',12)
% set(findall(f,'-property','Box'),'Box', 'off')
% % set(findall(ht,'-property','Interpreter'),'Interpreter','none')
% % set(findall(ht,'-property','TickLabelInterpreter'),'TickLabelInterpreter','none')
% set(f,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos=get(f,'Position');
% set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[pos(3) pos(4)])
% print(f, 'histo','-dpdf','-painters','-fillpage');
% print(f, 'histo','-dpng','-painters');
%% STUDIO come il modello predice - CONFUSION MATRIX
quantili = quantile(Tsub.deltaT, [.25 .5 .75 ]);
[DT,edg] = discretize(Tsub.deltaT, [0, quantili,max(Tsub.deltaT)]);
dedg = diff(edg);
mask = [1 0.6 0; 0.3 1 0.6; 0 0.3 1];

% confusion matrix predizione di classe

f=figure(4);clf;
tiledlayout('flow')
ax(1)=nexttile;
[C,order] = confusionmat(Realclasses,predictions,'order',names);
[hTOT, chiTOT, pTOT]=crosstab(Realclasses,predictions)
confusionchart(C,order)
Cw = C.*mask;
acc_tot(1) = sum(diag(C))/sum(sum(C)); % real accuracy
% acc_tot(2) = sum(sum(Cw))/sum(sum(C)); % weighted accuracy
ax(2) = nexttile;
rocObj = rocmetrics(Realclasses,scores,names);
plot(rocObj,ShowConfidenceIntervals=true)

picturewidth=45;
hw_ratio=0.5;
set(findall(f,'-property','FontSize'),'FontSize',12)
set(findall(f,'-property','Box'),'Box', 'off')
% set(findall(ht,'-property','Interpreter'),'Interpreter','none')
% set(findall(ht,'-property','TickLabelInterpreter'),'TickLabelInterpreter','none')
set(f,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos=get(f,'Position');
set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[pos(3) pos(4)])
print(f, 'roc','-dpdf','-painters','-fillpage');
print(f, 'roc','-dpng','-painters');
%%  in funzione di Delta_t + 
prec_foreachclass = {};sens_foreachclass = {};
% precision & sensitivity for each class
for i=1:numel(names)
    prec_foreachclass{1,i}=C(i,i)/(sum(C(:,i)));
    sens_foreachclass{1,i}=C(i,i)/(sum(C(i,:)));

end

f=figure(5);clf;
for j=1:max(DT)  
    
    % tiledlayout('flow');
	sel=(DT==j);
    ax(j)=nexttile;
    [Ci,Iorder] = confusionmat(Realclasses(sel),predictions(sel),'order',names);
    [hacc, chi(j), pACC(j)]=crosstab(Realclasses(sel),predictions(sel));
    acc_time(j,1) = sum(diag(Ci))/sum(sum(Ci));
    Ciw = Ci.*mask;
    acc_time(j,2) = sum(sum(Ciw))/sum(sum(Ci));
    for i=1:numel(names)
        prec_foreachclass{j+2,i}=Ciw(i,i)/(sum(Ciw(:,i)));
        sens_foreachclass{j+2,i}=Ciw(i,i)/(sum(Ciw(i,:)));
    end
    confusionchart(Ci,Iorder)
    title('Confusion matrix') 
end

picturewidth=45;
hw_ratio=0.5;
set(findall(f,'-property','FontSize'),'FontSize',12)
set(findall(f,'-property','Box'),'Box', 'off')
% set(findall(ht,'-property','Interpreter'),'Interpreter','none')
% set(findall(ht,'-property','TickLabelInterpreter'),'TickLabelInterpreter','none')
set(f,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos=get(f,'Position');
set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[pos(3) pos(4)])
print(f, 'cmtime','-dpdf','-painters','-fillpage');
print(f, 'cmtime','-dpng','-painters');
%% Accuracy/% in funzione del tempo
Pred=predictions==Realclasses; % si o no

f=figure(6);clf;
plot(cumsum(dedg),acc_time(:,1)*100,'s-','LineWidth',2); 
title('Percentage of accuracy (%)','FontWeight','bold');
xlabel('prediction ahead \Delta t [y]');
ylabel('Accuracy - percentage of correctly predicted classes (%)');
grid on;

 
% picturewidth=45;
% hw_ratio=0.5;
% set(findall(f,'-property','FontSize'),'FontSize',12)
% set(findall(f,'-property','Box'),'Box', 'off')
% % set(findall(ht,'-property','Interpreter'),'Interpreter','none')
% % set(findall(ht,'-property','TickLabelInterpreter'),'TickLabelInterpreter','none')
% set(f,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos=get(f,'Position');
% set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[pos(3) pos(4)])
% print(f, 'acc','-dpdf','-painters','-fillpage');
% print(f, 'acc','-dpng','-painters');

%% definizione di evento(i)
EDSSevent=TtrainRF.risk_class=='at risk'; % RISK CLASS 
MODELevent=predictions=='at risk';
pat=unique(TtrainRF.Paziente);

%% creazione tabella x SURVIVAL analysis
% if stratify
% TRFsub2 = TRFsub(TRFsub.eqCOUR=='PP',:);
% pat=unique(TRFsub.Paziente);

% follow
Tsurv=table('size',[numel(pat),5],'variabletypes',{'string','logical','double','logical','double'},'variablenames',{'Paziente','event','EventTime','model','ModelTime'});
id_patsurv = NaN(numel(pat),1);
% add index prima riga di ogni paziete , poi incolla a Tsurv. attenzione a
% T.Paziente; fitcox, pi hazardratio plotSurvival(coxmdl). e all abeta è il
% fattore di rischio , prossima volta ragazzi, aumento del rischio pe anno,
% aic per vedere diversi modelli
for j=1:numel(pat) % rullo sui pazienti
    w = TtrainRF.Paziente==pat(j); % trovo ogni paziente e la sua sottomatrice
    Tpat=TtrainRF(w,:); % sottomatrice del singolo paziente 
    id_patsurv(j)=find(w,1,'first'); % indice del paziente
	trueeventline=find(EDSSevent(w),1,'first'); % prendo il pezzo di vettore che ha info su eventi veri (edss at risk) per il singolo paziente e prendo la prima 
	modeleventline=find(MODELevent(w),1,'first'); % idem per il  modello
	truedateevent=Tpat.deltaT(trueeventline); % nella matrice del singolo paziente prendo la data dell'evento
	modeldateevent=Tpat.deltaT(modeleventline);
	trueevent=true;
	modelevent=true;
	if isempty(modeleventline)
		modeldateevent=Tpat.deltaT(end);
		modelevent=false;
	end
	if isempty(trueeventline)
		truedateevent=Tpat.deltaT(end);
		trueevent=false;
	end
	Tsurv(j,:)={pat(j) trueevent truedateevent modelevent modeldateevent};
end
Tbase = TtrainRF(id_patsurv,:);
Tsurv_complete = [Tsurv,Tbase(:,vars)];
% Tsurv.EventTime=duration(Tsurv.EventTime,'format','y');
% Tsurv.ModelTime=duration(Tsurv.ModelTime,'format','y');

%% KAPLAN MEIER analysis
% categoricalgroup=categorical(cat(1,repmat("true data",[height(Tsurv),1]),repmat("model",[height(Tsurv),1])));
categoricalgroup = categorical([repmat("true data",[height(Tsurv),1]);repmat("model",[height(Tsurv),1])]);
%t=cat(1,Tsurv.EventTime,Tsurv.ModelTime);
t = [Tsurv.EventTime;Tsurv.ModelTime];
% censored=cat(1,~Tsurv.event,~Tsurv.model);
censored=[~Tsurv.event;~Tsurv.model];
cmap=lines;
f=figure(7);clf;
KaplanMeier(categoricalgroup,t,censored,true,cmap);
xlabel(sprintf('# of years passed after %d yrs of monitoring',tstar));

ax = gca;
ax2=axes('Position',ax.Position,'XAxisLocation','top','YAxisLocation','right','Color','none','XLim',ax.XLim,'XTick',ax.XTick);
    ax2.YAxis.Visible='off';ax2.XAxis.Color='b';
    num=str2double(ax.XTickLabel);
    ev=zeros(1,numel(num));
    logicalgroup=dummyvar(categoricalgroup);
    logicalgroup=logical(logicalgroup(:,1));
    for j=2:numel(num)-1;ev(j)=nnz(t(logicalgroup)<=num(j));end    
    surv=nnz(logicalgroup)-ev;
    ax2.XTickLabel=cat(1,cellstr(num2str(surv(1:end-1)')),'n at risk',' ');
    title('Kaplan Meier curves')
% logrank([Tsurv.EventTime,~Tsurv.event],[Tsurv.ModelTime,~Tsurv.model]);
% 
% picturewidth=45;
% hw_ratio=0.5;
% set(findall(f,'-property','FontSize'),'FontSize',12)
% set(findall(f,'-property','Box'),'Box', 'off')
% % set(findall(ht,'-property','Interpreter'),'Interpreter','none')
% % set(findall(ht,'-property','TickLabelInterpreter'),'TickLabelInterpreter','none')
% set(f,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos=get(f,'Position');
% set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[pos(3) pos(4)])
% print(f, 'km_withsample','-dpdf','-painters','-fillpage');
% print(f, 'km_withsample','-dpng','-painters');
% % % 
%% vorrei visulaizzare quelli che hanno fup più lungo
patLONG = unique(TtrainRF.Paziente(TtrainRF.deltaT>7.3));
esempioRF = table();
f=figure(8);clf;
tiledlayout('flow')
for i=1:numel(patLONG)
    ax(i)=nexttile;
    w = TtrainRF.Paziente==patLONG(i);
    toplot = [TtrainRF(w,:),array2table(predictions(w))];
    esempioRF = cat(1,esempioRF,toplot);
    plot(toplot.deltaT,toplot.risk_class,'-o',toplot.deltaT,toplot.Var1,'-s')
    % hold on
    % xline(toplot.RELAPS,'-') 
    % hold off
    title(ax(i),sprintf('subject %s',patLONG(i)))
    % hold on
    % plot(toplot.deltaT,toplot.Var1)
end
l=legend('real class','predicted class')
l.Location='best';

% picturewidth=45;
% hw_ratio=0.5;
% set(findall(f,'-property','FontSize'),'FontSize',12)
% set(findall(f,'-property','Box'),'Box', 'off')
% % set(findall(ht,'-property','Interpreter'),'Interpreter','none')
% % set(findall(ht,'-property','TickLabelInterpreter'),'TickLabelInterpreter','none')
% set(f,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos=get(f,'Position');
% set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[pos(3) pos(4)])
% print(f, 'single','-dpdf','-painters','-fillpage');
% print(f, 'single','-dpng','-painters');

%% PREDETTI BENE E MALE
TtrainMED=[TtrainRF(base_disability(:,2),:),array2table(predictions(base_disability(:,2)))];
TtrainHIGH=[TtrainRF(base_disability(:,3),:),array2table(predictions(base_disability(:,3)))];

% TtrainHIGH=[TtrainRF,array2table(predictions)];

f=figure(9);clf;
tiledlayout('flow')
nexttile
histogram(TtrainHIGH.risk_class)
hold on
histogram(TtrainHIGH.Var1)
title('High EDSS (>5)')
legend('Real class', 'Predicted class')
nexttile
histogram(TtrainMED.risk_class)
hold on
histogram(TtrainMED.Var1)
title('Low EDSS [1-5]')
l=legend('Real class', 'Predicted class')
l.Location='best';

[Chigh,~] = confusionmat(TtrainHIGH.risk_class,TtrainHIGH.Var1,'order',names);
[Cmed,~] = confusionmat(TtrainMED.risk_class,TtrainMED.Var1,'order',names);

% [hTOT, chiTOT, pTOT]=crosstab(Realclasses,predictions)
nexttile
confusionchart(Chigh,order)
nexttile
confusionchart(Cmed,order)

% picturewidth=45;
% hw_ratio=0.5;
% set(findall(f,'-property','FontSize'),'FontSize',12)
% set(findall(f,'-property','Box'),'Box', 'off')
% % set(findall(ht,'-property','Interpreter'),'Interpreter','none')
% % set(findall(ht,'-property','TickLabelInterpreter'),'TickLabelInterpreter','none')
% set(f,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos=get(f,'Position');
% set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[pos(3) pos(4)])
% print(f, 'uniquemod','-dpdf','-painters','-fillpage');
% print(f, 'uniquemod','-dpng','-painters');

%% good vs bad

goodMED = TtrainMED(TtrainMED.risk_class==TtrainMED.Var1,:);
badMED = TtrainMED(~(TtrainMED.risk_class==TtrainMED.Var1),:);

goodHIGH = TtrainHIGH(TtrainHIGH.risk_class==TtrainHIGH.Var1,:);
badHIGH = TtrainHIGH(~(TtrainHIGH.risk_class==TtrainHIGH.Var1),:);

% medium
f=figure(10);clf;
ht=tiledlayout('flow');
for i=1:numel(independent{1})
    [~,idg]=sort(goodMED.(independent{1}{i}));
    [~,idb]=sort(badMED.(independent{1}{i}));
    nexttile
    yb=ksdensity(badMED.(independent{1}{i}),badMED.(independent{1}{i}));yg=ksdensity(goodMED.(independent{1}{i}),goodMED.(independent{1}{i}));
    plot(sort(badMED.(independent{1}{i})),yb(idb),sort(goodMED.(independent{1}{i})),yg(idg))
      xlabel('Score')
      ylabel('Probability Density Function')
    [~,p,ci,~] = ttest2(yb,yg);
    title(sprintf('Confer %s - p value %1.3f',independent{1}{i},p),'Interpreter','none')
end
legend('Bad predicted','Well predicted')
title(ht,'EDSS [1-5]','Fontweight','bold')


%high
f=figure(11);clf;
ht=tiledlayout('flow');
for i=1:numel(independent{1})
    [~,idg]=sort(goodHIGH.(independent{1}{i}));
    [~,idb]=sort(badHIGH.(independent{1}{i}));
    nexttile
    yb=ksdensity(badHIGH.(independent{1}{i}),badHIGH.(independent{1}{i}));yg=ksdensity(goodHIGH.(independent{1}{i}),goodHIGH.(independent{1}{i}));
    plot(sort(badHIGH.(independent{1}{i})),yb(idb),sort(goodHIGH.(independent{1}{i})),yg(idg))
    xlabel('Score')
    ylabel('Probability Density Function')

    [~,p,ci,~] = ttest2(yb,yg);
    title(sprintf('Confer %s - p value %1.3f',independent{1}{i},p),'Interpreter','none')
end
legend('Bad predicted','Well predicted')
title(ht,'EDSS [>5]','Fontweight','bold')

% picturewidth=45;
% hw_ratio=0.5;
% set(findall(f,'-property','FontSize'),'FontSize',12)
% set(findall(f,'-property','Box'),'Box', 'off')
% % set(findall(ht,'-property','Interpreter'),'Interpreter','none')
% % set(findall(ht,'-property','TickLabelInterpreter'),'TickLabelInterpreter','none')
% set(f,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos=get(f,'Position');
% set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[pos(3) pos(4)])
% print(f, 'hbase','-dpdf','-painters','-fillpage');
% print(f, 'hbase','-dpng','-painters');

% slopes

% medium
f=figure(12);clf;
ht=tiledlayout('flow');
for i=1:numel(independent{2})
    [~,idg]=sort(goodMED.(independent{2}{i}));
    [~,idb]=sort(badMED.(independent{2}{i}));
    nexttile
    yb=ksdensity(badMED.(independent{2}{i}),badMED.(independent{2}{i}));yg=ksdensity(goodMED.(independent{2}{i}),goodMED.(independent{2}{i}));
    plot(sort(badMED.(independent{2}{i})),yb(idb),sort(goodMED.(independent{2}{i})),yg(idg))
    ylabel('Probability Density Function')
    xlabel(sprintf('%i year-slope',tstar))
    [~,p,ci,~] = ttest2(yb,yg);
    title(sprintf('Confer %s - p value %1.3f',independent{2}{i},p),'Interpreter','none')
end
legend('Bad predicted','Well predicted')
title(ht,'EDSS [1-5]','Fontweight','bold')

% picturewidth=45;
% hw_ratio=0.5;
% set(findall(f,'-property','FontSize'),'FontSize',12)
% set(findall(f,'-property','Box'),'Box', 'off')
% % set(findall(ht,'-property','Interpreter'),'Interpreter','none')
% % set(findall(ht,'-property','TickLabelInterpreter'),'TickLabelInterpreter','none')
% set(f,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos=get(f,'Position');
% set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[pos(3) pos(4)])
% print(f, 'mslope','-dpdf','-painters','-fillpage');
% print(f, 'mslope','-dpng','-painters');


%high
f=figure(13);clf;
ht=tiledlayout('flow');
for i=1:numel(independent{2})
    [~,idg]=sort(goodHIGH.(independent{2}{i}));
    [~,idb]=sort(badHIGH.(independent{2}{i}));
    nexttile
    yb=ksdensity(badHIGH.(independent{2}{i}),badHIGH.(independent{2}{i}));yg=ksdensity(goodHIGH.(independent{2}{i}),goodHIGH.(independent{2}{i}));
    plot(sort(badHIGH.(independent{2}{i})),yb(idb),sort(goodHIGH.(independent{2}{i})),yg(idg))
   xlabel(sprintf('%i year-slope',tstar))
  ylabel('Probability Density Function')
    [~,p,ci,~] = ttest2(yb,yg);
    title(sprintf('Confer %s - p value %1.3f',independent{2}{i},p),'Interpreter','none')
end
l=legend('Bad predicted','Well predicted');
l.Location='best';
title(ht,'EDSS [>5]','Fontweight','bold')

%% partial dependence plot
fitmodel = fitcensemble(TtrainRF,'risk_class',...
		'prior','uniform','method','RUSBoost',...
		'predictornames',vars,'learners',tr,'ClassNames',names,...
		'NumLearningCycles',250);
rsLoss = resubLoss(fitmodel,'Mode','Cumulative');

f=figure(14);clf;
plot(rsLoss);
xlabel('Number of Learning Cycles');
ylabel('Resubstitution Loss');

%
%tiledlayout('flow')
% for i=1:numel(vars) 
%     f=figure(14+i);clf;    
%     plotPartialDependence(fitmodel,vars(i),fitmodel.ClassNames)
%     % picturewidth=45;
%     % hw_ratio=0.5;
%     % set(findall(f,'-property','FontSize'),'FontSize',12)
%     % set(findall(f,'-property','Box'),'Box', 'off')
%     % % set(findall(ht,'-property','Interpreter'),'Interpreter','none')
%     % % set(findall(ht,'-property','TickLabelInterpreter'),'TickLabelInterpreter','none')
%     % set(f,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
%     % pos=get(f,'Position');
%     % set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[pos(3) pos(4)])
%     % print(f, 'deltaT','-dpdf','-painters','-fillpage');
%     % print(f, 'deltaT','-dpng','-painters');
% end



%% focus su DRUG and REHAB TREAT
up2=unique(T.Paziente);
clear delta %= years(0) %NaN(height(T),1);
TtrainRF.Drug = string(NaN(height(TtrainRF),1));
TtrainRF.TREA01 = (NaN(height(TtrainRF),1));
TtrainRF.TREA02 = (NaN(height(TtrainRF),1));
TtrainRF.TREA03 = (NaN(height(TtrainRF),1));
TtrainRF.TREA04 = (NaN(height(TtrainRF),1));
TtrainRF.TREA05 = (NaN(height(TtrainRF),1));
TtrainRF.TREA06 = (NaN(height(TtrainRF),1));
TtrainRF.TREA07 = (NaN(height(TtrainRF),1));
TtrainRF.TREA08 = (NaN(height(TtrainRF),1));
TtrainRF.TREA09 = (NaN(height(TtrainRF),1));


for j=1:numel(up2)
	w=T.Paziente==up2(j);
    wf = find(w);
    Tpat = T(w,:);
    T.fup(w) = years(Tpat.Date-Tpat.Date(1));
    visitbeforetstar=(Tpat.Date-min(Tpat.Date))<=years(tstar);
    visittime=duration(Tpat.Date-min(Tpat.Date),'format','y');
    if range(Tpat.Date)>years(tstar)    
       delta(wf(~visitbeforetstar),1)=visittime(~visitbeforetstar)-years(tstar);
    end
end
T.deltaT = delta;

pat = unique(TtrainRF.Paziente);

for i=1:numel(pat)
    wT = T.Paziente == pat(i);
    w = TtrainRF.Paziente == pat(i);
    Tw = T(wT,:);
    Tpat = TtrainRF(w,:);
    for j=1:size(Tpat,1)
        l = Tw.deltaT == years(Tpat.deltaT(j)); 
        if nnz(l)>0
            Tpat.Drug(j) = Tw.DRUG(l);
            Tpat.TREA01(j) = Tw.TREA01(l);
            Tpat.TREA02(j) = Tw.TREA02(l);
            Tpat.TREA03(j) = Tw.TREA03(l);
            Tpat.TREA04(j) = Tw.TREA04(l);
            Tpat.TREA05(j) = Tw.TREA05(l);
            Tpat.TREA06(j) = Tw.TREA06(l);
            Tpat.TREA07(j) = Tw.TREA07(l);
            Tpat.TREA08(j) = Tw.TREA08(l);
            Tpat.TREA09(j) = Tw.TREA09(l);
        end
    end
    TtrainRF.Drug(w)=Tpat.Drug;
    TtrainRF.TREA01(w) = Tpat.TREA01;
    TtrainRF.TREA02(w) = Tpat.TREA02;
    TtrainRF.TREA03(w) = Tpat.TREA03;
    TtrainRF.TREA04(w) = Tpat.TREA04;
    TtrainRF.TREA05(w) = Tpat.TREA05;
    TtrainRF.TREA06(w) = Tpat.TREA06;
    TtrainRF.TREA07(w) = Tpat.TREA07;
    TtrainRF.TREA08(w) = Tpat.TREA08;
    TtrainRF.TREA09(w) = Tpat.TREA09;
end


%% ora considero ben predetti quelli che hanno almeno 50% di predizioni corrette
TtrainSUB = TtrainRF(~isnan(TtrainRF.IDrf),:);
clear corretti;
QC = [pat, repmat("to define",numel(pat),1)];
for i=1:numel(pat)
   w = TtrainRF.Paziente == pat(i); 
   Tpat = TtrainRF(w,:);
   predictions_QC = predictions(w);
   Real_QC = Realclasses(w);
   corretti(i,1) = sum(predictions_QC==Real_QC);
   corretti(i,2) = numel(predictions_QC)-corretti(i,1);
   if corretti(i,1) >= 0.5*(numel(predictions_QC))
       QC(i,2)= "good";
   else 
       QC(i,2)= "bad";
   end
end
 % good 331
 % bad 294
%% solo per quelli con MRI
pat = unique(TtrainSUB.Paziente);
clear corretti1;
QCsub = [unique(TtrainSUB.IDrf),pat, repmat("to define",numel(pat),1)];
for i=1:numel(pat)
   w = TtrainSUB.Paziente == pat(i); 
   Tpat = TtrainSUB(w,:);
   predictions_QC = predictions(w);
   Real_QC = Realclasses(w);
   corretti1(i,1) = sum(predictions_QC==Real_QC);
   % corretti1(i,2) = numel(predictions_QC)-corretti(i,1);
   if corretti1(i,1) >= 0.5*(numel(predictions_QC))
       QCsub(i,3)= "good";
   else 
       QCsub(i,3)= "bad";
   end
end
 % good
 % bad
 figure(16);clf;
 histogram(categorical(QCsub(:,2)))
 %% logistic regression con var indip carico lesionale 

 f2 = fullfile('C:','Users','federica.diantonio','Desktop','ME','DOTTORATO','MRI','DB_PROM_QCmri.xlsx');
 QC_MRI = readtable(f2);
 f1 = fullfile('C:','Users','federica.diantonio','Desktop','ME','DOTTORATO','MRI','IDs_name.xlsx');
 id_name = readtable(f1);

 QC_MRI_nomissing =QC_MRI(~(ismissing(QC_MRI.Les_V)),:);
 lesions = [string(QC_MRI_nomissing.name),QC_MRI_nomissing.Les_V,repmat(0,size(QC_MRI_nomissing.Les_V,1),1)];
 lesions = array2table(lesions);
 lesions.lesions3 = str2double(lesions.lesions3 );
 id_name.name = string(id_name.name);
 id_name.name = strrep(id_name.name,'_',' ');
 % add a column con ID 

 for i=1:size(lesions,1)
     w = id_name.name == lesions.lesions1(i);
     wf = find(w);
     if nnz(w)>0

        lesions.lesions3(i)=id_name.ID(w);
     end

 end
 lesions.Properties.VariableNames=["Paziente" "Carico lesionale" "ID"];
 lesions.("Carico lesionale")=str2double(lesions.("Carico lesionale"));
  % mi serve il carico lesionale dei 60 di cui ho good bad prediction 
QCsub = array2table(QCsub); QCsub.Properties.VariableNames=["ID" "Paziente" "Classification"];
QCsub.ID = str2double(QCsub.ID);

for i=1:size(QCsub,1)
    w = lesions.ID == QCsub.ID(i);
    wf = find(w);
    V_lesions(i)=lesions.("Carico lesionale")(wf)
end
QCsub = addvars(QCsub,V_lesions');
QCsub=renamevars(QCsub,'Var4','V_lesions');

%% logistic regression

boxplot(QCsub.V_lesions,categorical(QCsub.Classification))
[p,h]=ranksum(QCsub.V_lesions(QCsub.Classification=="good"),QCsub.V_lesions(QCsub.Classification=="bad"));


%%  GUARDIAMO I DMT DEI GOOD PREDICTED 
wg = QC(QC(:,2)=="good",1);Tgood = TtrainRF(1,:);
wb = QC(QC(:,2)=="bad",1);Tbad = TtrainRF(1,:);
for i=1:numel(wg)
    Tgood = [Tgood;TtrainRF(TtrainRF.Paziente == wg(i),:)];    
end
for i=1:numel(wb)
    Tbad = [Tbad;TtrainRF(TtrainRF.Paziente == wb(i),:)];
end
Tgood(1,:)=[];
Tbad(1,:)=[];
%% DMT 
figure(14);clf;
ht = tiledlayout('flow');
ax(1)=nexttile;
histogram(categorical(Tgood.Drug),'Normalization','pdf')
title("Well predicted")
grid on;

ax(2)=nexttile;
histogram(categorical(Tbad.Drug),'Normalization','pdf')
linkaxes([ax(1) ax(2)], 'y')
title("Bad predicted")
grid on;

%%
treat = ["TREA01" "TREA02" "TREA03" "TREA04" "TREA05" "TREA06" "TREA07" "TREA08" "TREA09"];
figure(15);clf;
tiledlayout('flow')
ax(1)=nexttile
for i=1:9
    bar(i,sum(Tgood.(treat(i)),'omitmissing')/size(Tgood,1))
    hold on;
end
hold off
title(ax(1),"Well predicted")
set(ax(1),'XTick',1:numel(treat),'XTickLabel',treat)
grid on;

ax(2)=nexttile;
for i=1:9
    bar(i,sum(Tbad.(treat(i)),'omitmissing')/size(Tbad,1))
    hold on;
end
hold off;
title("Bad predicted")
set(ax(2),'XTick',1:numel(treat),'XTickLabel',treat)
grid on;
linkaxes([ax(1) ax(2)], 'y')