clearvars;
set(0,'DefaultFigureWindowStyle','docked');
Namefile = 'data_promo2.mat';
FolderName = 'C:\Users\federica.diantonio\desktop\ME\DOTTORATO\';
File = fullfile(FolderName, Namefile);
load(File);

%% sistemazione variabili

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

%% un po' di nomi 
tstar=2;      % tempo di baseline                                                             
% nvisit=tstar*3;
% theorvisittime=years(linspace(0,tstar,nvisit));
Ddependent={'Delta_EDSSTOT' 'Delta_FIM_TOT' 'Delta_MFISTOT' 'Delta_LIFETOT'};
Class_base={'base_risk'};
Bdependent={'EDSSTOT' 'FIM_TOT' 'MFISTOT' 'LIFETOT'};
independent0={'EDSSTOT' 'FIM_TOT' 'BMI' 'ABILHTOT' 'HADSTOT' 'LIFETOT' 'MFISTOT' 'SDMTTOT' 'OAB_QTOT'};
problem={'deltaT'};
aux={'Paziente'};
external = {'eqCOUR' 'OtherDisease' 'RELAPS'};
comorb = {'metabolic' 'cardiovascular' 'psichic' 'onco/neuro' 'other'};
proof = {'dummy'};
demo = {'AGE' 'GEND' 'Duration_from_DD'};

clabel=["base" "slope"];
fupslabelfun=@(x,n) sprintf('%s_%s',x,clabel(n));

for j=1:2 % baseline e trend
	independent{j}=cellfun(@(x) fupslabelfun(x,j),independent0,'UniformOutput',false);
end

%% Ttrain con 1 y fup
Ttrain=table(); % creaimo matrice con info mean e slope del primo anno e solo le righe dei delta t >1 in cui metto la variazione di edss rispetto al mean score del primo anno
up=unique(T.Paziente);
for j=1:numel(up)
	w=T.Paziente==up(j);
	Tpat=T(w,:);
	[~,pos]=sort(Tpat.Date,'ascend');
	Tpat=Tpat(pos,:); % matrice del paziente ordinata temporalmente
	if range(Tpat.Date)>years(tstar) % se il range di fup del paziente è maggiore di 1 anno
		visitbeforetstar=(Tpat.Date-min(Tpat.Date))<=years(tstar); % rendo le visite in formato anno giusto       
		visittime=duration(Tpat.Date-min(Tpat.Date),'format','y');  % prendo le visite prima dell'anno		
		var=zeros(numel(independent0),2);% costruisco una matrice con 2 colonne (baseline e trend) e righe tante quante le variabili predittori
		for v=1:numel(independent0) % loppo sulle variabili predittrici
			tstarvisittime=visittime(visitbeforetstar);% prendo le date delle visite pre anno
			visit=Tpat.(independent0{v})(visitbeforetstar); % di ogni variabile prendo i valori pre anno
			visitbase=mean(visit,'omitnan'); % media degli score del primo anno
			nan=isnan(visit);% vedo i nan della variabile
			switch numel(visit) % per ogni numero di visite effettive che ha il paziente
				case 1 % se ha una visita sola
					var(v,:)=[NaN,NaN];% BASELINE e DERIVATA NULLA
				case 2 % se ha due visite 
					switch nnz(nan)  % numero di missin values
						case 0 % NESSUN MISSING VALUE
							lm=polyfit(years(tstarvisittime),visit,1);
							var(v,:)=[visitbase,lm(1)];
						case 1 % 1 MISSING VALUE
							var(v,:)=[NaN,NaN];
						case 2 % 2 MISSING VALUE
							var(v,:)=[NaN,NaN];
                    end % se sono tutti nan lasciamo i nan e poi il soggetto se ne andrà 
				otherwise
					if nnz(~nan)>=2
						p=polyfit(years(tstarvisittime(~nan)),visit(~nan),1);
						f=find(visit(~nan),1,'first');
						var(v,:)=[visitbase,p(1)];
					else
						var(v,:)=[NaN,NaN];
					end

			end
        end
        % deltaT=visittime(~visitbeforetstar)-years(tstar);  % GIA' COSTRUISCE MATRICE A PARTIRE DA TSTAR IN POI
        deltaT=visittime(~visitbeforetstar);
		Vdelta=zeros(nnz(~visitbeforetstar),numel(Bdependent));
        % Vclass=zeros(nnz(visitbeforetstar),numel(Bdependent));

		for v=1:numel(Bdependent)
			% Vdelta(:,v)=table2array(Tpat(~visitbeforetstar,Bdependent(v))-Tpat(find(visitbeforetstar,1,'last'),Bdependent(v)));
            Vdelta(:,v)=table2array( Tpat(~visitbeforetstar,Bdependent(v))-mean(Tpat(visitbeforetstar,Bdependent(v))) );
            % Vclass(:,v)=table2array( Tpat(visitbeforetstar,Bdependent(v))-Tpat(find(visitbeforetstar,1,"first"),Bdependent(v)));
            % Vclass(:,v)=table2array( Tpat(visitbeforetstar,Bdependent(v))-mean(Tpat(visitbeforetstar,Bdependent(v))));
        end        
		Cdelta=cumsum(Vdelta)./years(deltaT); %Vdelta;
        % Cclass=cumsum(Vclass)./years(visittime(visitbeforetstar));
		varbase=repmat(var(:)',[nnz(~visitbeforetstar),1]); %mette in colonna tutte base e poi tutte slope e poi trasversa e ripete pe rtutt le visite post 1 year
		% classbase = repmat(Cclass(end,:),[nnz(~visitbeforetstar),1]);
        Tpatbase=table(deltaT);
		Tpatbase=[Tpatbase,array2table([varbase,Cdelta],'variablenames',[independent{:},Ddependent])];
		Tpatbase=[Tpat(~visitbeforetstar,[aux,demo,external]),Tpatbase];
		Ttrain=[Ttrain;Tpatbase];
	end
end

%% pulizia dei missing
FOCUS="Delta_EDSSTOT";%"Delta_MFISTOT"; %"Delta_EDSSTOT";%"Delta_FIM_TOT";%"Delta_LIFETOT";

missinglines=any(ismissing(Ttrain),2);
Ttrain(missinglines,:)=[];
Ttrain(Ttrain.RELAPS>0,:)=[];
out_id = Ttrain.(FOCUS) < quantile(Ttrain.(FOCUS),0.05) | Ttrain.(FOCUS)> quantile(Ttrain.(FOCUS),0.99);
Ttrain(out_id,:)=[];
%% DIVIDO IN CLASSI EDSS SULLA BASELINE

id_lowbase = Ttrain.EDSSTOT_base==0; % 11
id_medbase = Ttrain.EDSSTOT_base>0 & Ttrain.EDSSTOT_base<=5;%[1-5],2251
id_highbase = Ttrain.EDSSTOT_base>5; %  2047
base_disability=[id_lowbase,id_medbase,id_highbase];
labels = ["basso" "medium" "high"];
cat_label = table('Size', [size(base_disability,1) 1], 'VariableTypes',"string");

% cat_label=table(height(Ttrain),1);
for i=1:3
    cat_label(base_disability(:,i),1)=array2table(labels(i));
end
cat_label = categorical(table2array(cat_label),labels,'ordinal',true);
% Ttrain=addvars(Ttrain,cat_label);


% names = ["negligible" "moderate" "severe"];
names = ["negligible" "severe"];  
% L1 = [0.1 1.5];L2 = [0.1 1];L3 = [0.1 0.5]; % discretize prende il limite dx
L1 = 0.1;L2 = 0.1;L3 = 0.1; % discretize prende il limite dx
risk_class(id_lowbase,1)=(discretize(Ttrain.(FOCUS)(id_lowbase),[-Inf L1 Inf],names)); % basso edss
risk_class(id_medbase,1)=discretize(Ttrain.(FOCUS)(id_medbase),[-Inf L2 Inf],names);   % medio edss
risk_class(id_highbase,1)=discretize(Ttrain.(FOCUS)(id_highbase),[-Inf L3 Inf],names); % alto edss

risk_class = categorical(risk_class,names,'ordinal',true);

% % Ttrain=addvars(Ttrain,risk_class);
% temp_risk=[risk_class{1},risk_class{2},risk_class{3}];
%risk = array2table(temp_risk);risk.Properties.VariableNames=["risk_classMFIS" "risk_classFIM" "risk_classLIFE"];
% Ttrain.dummy= rand(height(Ttrain),1);
%% tutto parte dalla scelta della variabile da predire

% comorbidità
% load comorbidities.mat;
cat_morbius=categorize_comorbidities(Ttrain.OtherDisease);
comorbidities=sum(cat_morbius,2);
Ttrain.comorbidities=categorical(comorbidities.sum,0:4,{'none' 'low' 'medium' 'high' 'severe'},'ordinal',true);
dummy= rand(height(Ttrain),1);  % aggiunta dummy var 
base_risk = NaN(height(Ttrain),1);

TtrainRF = addvars(Ttrain,base_risk, cat_label ,dummy, risk_class);
TtrainRF.deltaT=years(TtrainRF.deltaT);                                     % trasforma in double il numero di anni duration


%% da Ttrain costruisco la prima 

pat = unique(TtrainRF.Paziente);
t_range = 2.5; 
for i=1:numel(pat)
    w = TtrainRF.Paziente==pat(i);
    Tpat=TtrainRF(w,:);
    if Tpat.deltaT(1) <= t_range && Tpat.RELAPS(1) == 0 % delta T conta il tempo che apssa dopo il primo anno di fup (o tstar di fup)
        TtrainRF.base_risk(w)=Tpat.risk_class(1);
        % cat(1,Tclass,Tpat(1,:));
    end
end
TtrainRF(isnan(TtrainRF.base_risk),:)=[];

% edss 443 pat and 2815 clinical records
%% predizione a 1y
target_time=[3, 4, 5]; % visita piu vicina a 2 anni a meno di 6 mesi 
%target_time=2;
time_range=0.5; % scarto consentito dall'anno di visita = 6 mesi
clear TRFtime sub_risk_base;
TRFtime = cell(numel(target_time),1);
pat = unique(TtrainRF.Paziente);

for k=1:numel(target_time)
    for j=1:numel(pat)	
        w = TtrainRF.Paziente == pat(j);
        % if nnz(w)>0
        Tpat=TtrainRF(w,:);
	    [~,sel]=min(abs((Tpat.deltaT-target_time(k))));
	    % w=find(Tpat.deltaT<time_range,1,'first');
	    if abs(Tpat.deltaT(sel)-target_time(k))<time_range
            % sub_risk_base(j,:)=risk_base(pat_with_classbase==pat(j));
		    TRFtime{k,1}(j,:)=Tpat(sel,:);
        end
    end
    TRFtime{k}=rmmissing(TRFtime{k});
    [freq(k,:),names(k,:)]=histcounts(TRFtime{k}.risk_class);
end


%% primo step: assegnare la classe di rischio istantanea
out = {'risk_class'};
% base = [{'EDSSTOT_base'} {'FIM_TOT_base'}    {'BMI_base'}    {'ABILHTOT_base'}    {'HADSTOT_base'}    {'LIFETOT_base'}    {'MFISTOT_base'}    {'SDMTTOT_base'}    {'OAB_QTOT_base'}];
base = [{'EDSSTOT_base'} {'FIM_TOT_base'}];
slope = [{'EDSSTOT_slope'} {'FIM_TOT_slope'}];
% slope = [{'EDSSTOT_slope'} {'FIM_TOT_slope'}    {'BMI_slope'}    {'ABILHTOT_slope'}    {'HADSTOT_slope'}    {'LIFETOT_slope'}    {'MFISTOT_slope'}    {'SDMTTOT_slope'}    {'OAB_QTOT_slope'}];
vars = [base slope 'base_risk'];
leaf = 1;
%vars=[demo independent{:} 'base_risk' 'eqCOUR' 'comorbidities']; %cat_morbius.Properties.VariableNames proof
% vars = [demo base slope 'base_risk' 'eqCOUR' 'comorbidities' ];
Ttime = cell(numel(target_time),1);OOB = cell(numel(target_time),1);scores = cell(numel(target_time),1);predictions=cell(1,numel(target_time));
names = unique(TRFtime{k}.risk_class);
for k=1:numel(target_time)
    Ttime{k} = TRFtime{k}(:,[vars,out]);    
    % cv=cvpartition(height(TtrainRF),'Kfold',20);
    cv=cvpartition(height(Ttime{k}),'leaveout');
    oobvi=zeros(cv.NumTestSets,numel(vars));
    scores{k}=zeros(height(Ttime{k}),numel(names));
    [sc,pred]=deal(cell(cv.NumTestSets,1));
    % pred=deal(cell(cv.NumTestSets,1))
    ht=tic;
    fprintf('starting CV ...\n');
    tr = templateTree('MinLeafSize',leaf,'surrogate','off',"PredictorSelection","interaction-curvature","SplitCriterion","gdi");

    parfor j=1:cv.NumTestSets   
	    fitmodel = fitcensemble(Ttime{k}(training(cv,j),:),'risk_class',...
		'prior','uniform','method','RUSBoost',...
		'predictornames',vars,'learners',tr,'ClassNames',names,...
		'NumLearningCycles',200); 
    % fitmodel = fitrensemble(TtrainRF(training(cv,j),:),'Delta_EDSSTOT','method','Bag',...
	% 	'predictornames',vars,'learners',tr,...
	% 	'NumLearningCycles',200); 
		[pred{j},sc{j}]=predict(fitmodel,Ttime{k}(test(cv,j),vars)); % scores ha cinque colonne perchè ho 3 possibili classi                       
	    oobvi(j,:)=predictorImportance(fitmodel);
    end
    OOB{k}=oobvi;
    predictions{k}=Ttime{k}.risk_class;

    for j=1:cv.NumTestSets
	    wtest = test(cv,j);
	    predictions{k}(wtest)=pred{j};
	    scores{k}(wtest,:)=sc{j};
    end
    dt=seconds(toc(ht));
    fprintf('done in %s minutes\n',minutes(dt));
end
%% accuracy

figure(2);clf;
tiledlayout('flow');

stats = cell(numel(target_time),1);
% if p<0.05
%     disp('keep going!');
% end
for i=1:numel(target_time)
    ax(i) = nexttile;
    [C,order]=confusionmat(Ttime{i}.risk_class,predictions{i},'order',names);confusionchart(C,order);
    accu(i)=sum(diag(C))/sum(C,'all');
    [~,~,p(i)]=crosstab(Ttime{i}.risk_class,predictions{i});
    title(sprintf('Accuracy %1.1f %% (p value %1.3f) T=%d years',accu(i)*100,p(i),target_time(i)));
    ax = gca; 
    ax.FontSize = 12; 
    stats{i} = statsOfMeasure(C, 1);
end
writetable(stats{1},'acc_etc1.xlsx')
writetable(stats{2},'acc_etc2.xlsx')
writetable(stats{3},'acc_etc3.xlsx')
% writetable(stats{4},'acc_etc4.xlsx')
export_fig cm_all -png -jpg -pdf -transparent


%% plot Variable Importance - histogram of average importance with bar errors
SD_importance=cell(numel(target_time),1);mean_importance=cell(numel(target_time),1);
figure(4);clf;
ht=tiledlayout(1,numel(target_time))
for k=1:numel(target_time)
    % id_dummy = find(string(vars) == "dummy");
    % pd=fitdist(OOB{k}(:,id_dummy),'normal');
    % cut=pd.mu+3*pd.sigma;
    % SD_importance = zeros(size(OOB{k},2),1);mean_importance= zeros(size(OOB{k},2),1);
    for i=1:size(OOB{k},2)
	    [SD_importance{k}(i),mean_importance{k}(i)] = std(OOB{k}(:,i));
    end
    x= 1:size(OOB{k},2);
    if k==1
        [POS(k,:),IDX(k,:)]=sort(mean_importance{k},'descend');    
        ax(k) = nexttile;        
        barh(x,POS(k,:),'k','FaceAlpha',0.2);
        hold on
        errorbar(POS(k,:),x,SD_importance{k}(IDX(k,:)),'k','linestyle','none','LineWidth',3)';
        % hold on
        % xline(cut,'-',{'Gaussian estimated cut (3 \sigma)',sprintf('x = %1.2f',cut)},'LineWidth',3) 
        hold off
        title(ax(k),sprintf('T = %d years',target_time(k)))
        h = gca;
        h.YTick=1:numel(vars);
        h.YTickLabel = vars(IDX(k,:));
        h.TickLabelInterpreter = "none";
        % ax = gca; 
        h.FontSize = 12; 
    else 
        ax(k) = nexttile;
        % x= 1:size(OOB{k},2);
        barh(x,mean_importance{k}(IDX(1,:)),'k','FaceAlpha',0.2);
        hold on
        errorbar(mean_importance{k}(IDX(1,:)),x,SD_importance{k}(IDX(1,:)),'k','linestyle','none','LineWidth',3)';
        % hold on
        % xline(cut,'-',{'Gaussian estimated cut (3 \sigma)',sprintf('x = %1.2f',cut)},'LineWidth',3) 
        hold off
        title(ax(k),sprintf('T = %d years',target_time(k)))
        h = gca;
        h.YTick=1:numel(vars);
        h.YTickLabel = vars(IDX(1,:));
        % h.XTickLabelRotation = 45;
        h.TickLabelInterpreter = "none";
        % xlabel("Predictor Importance Estimates","FontWeight","bold")
        % ylabel("Predictors","FontWeight","bold")
        % ax = gca; 
        h.FontSize = 12; 

    end
    
end
%grid on;
linkaxes([ax(1) ax(2) ax(3)],'x')
xlabel(ht,"Predictor Importance Estimates","FontWeight","bold")
ylabel(ht,"Predictors","FontWeight","bold")
% export_fig feat_impo_dummy -png -jpg -pdf -transparent

% feat_impo=[POS,SD_importance(IDX)];feat_impo=array2table(feat_impo);
% feat_impo=addvars(feat_impo,vars(IDX)','before','feat_impo1');
% feat_impo.Properties.VariableNames=["predictor" "Mean importance" "SD"];
% writetable(feat_impo,'impo.xlsx')
export_fig histo -png -jpg -pdf -transparent

%% roc curves
% rocObj = rocmetrics(TRFtime.risk_class,scores,names);
% plot(rocObj,ShowConfidenceIntervals=true,AverageROCType="weighted") 
% [Xsafe,Ysafe,Tsafe,AUCsafe] = perfcurve(Ttime.risk_class,scores(:,1),'low');
% [Xmod,Ymod,Tmod,AUCmod] = perfcurve(TRFtime.risk_class,scores(:,2),'moderate');


figure(5);clf;
Tauc = cell(numel(target_time),3);
for k=1:numel(target_time)
    [Tauc{k,1}(:,1),Tauc{k,1}(:,2),Tauc{k,2},Tauc{k,3}] = perfcurve(Ttime{k}.risk_class,scores{k}(:,2),'severe'); %0.8234 PER t=3
end
% ax(2) = nexttile;
h(1) = plot(Tauc{1,1}(:,1),Tauc{1,1}(:,2),'LineWidth',2,'Color','black');
hold on;
h(2) = plot(Tauc{2,1}(:,1),Tauc{2,1}(:,2),'LineWidth',2,'LineStyle','--','Color','black');
hold on;
h(3) = plot(Tauc{3,1}(:,1),Tauc{3,1}(:,2),'LineWidth',2,'LineStyle',':','Color','black');
% hold on; 
% h(4) =plot(Tauc{4,1}(:,1),Tauc{4,1}(:,2),'LineWidth',2,'LineStyle','-.','Color','black');


xlabel('False positive rate'); ylabel('True positive rate');
title('ROC Curves - Reduced model')
% title(sprintf('AUC %1.2f',AUCrisk));
axis square;
hold on;
% rocObj = rocmetrics(Realclasses,scores,names);
% plot(rocObj,ShowConfidenceIntervals=true) 
x = linspace(0,1);
y = linspace(0,1);
h4=plot(x,y,'LineStyle','-.','Color','k');
hold off;
%(T=%d yrs)',target_time
%legend(h1,sprintf('T = %d',target_time(1)),sprintf('T = %d',target_time(2)),sprintf('T = %d',target_time(3)),'Location','southeast')
% legend([h1 h2 h3],'Low','Moderate','Severe','Location','southeast')
legend([h(1) h(2) h(3) ],'T=2yrs','T=3yrs','T=4yrs','T=5yrs','Location','southeast')

ax = gca; 
ax.FontSize = 12; 

export_fig roc_all_BW4 -png -jpg -pdf -transparent
mat = [Tauc{1,3},size(Ttime{1},1),freq(1,:);Tauc{2,3}, size(Ttime{2},1),freq(2,:);Tauc{3,3}, size(Ttime{3},1),freq(3,:)];
writematrix(mat,'auc_bis.xlsx') % ...size(Ttime{1},1);freq'],'auc.xlsx');; Tauc{4,3}, size(Ttime{4},1),freq(4,:)
%% definizione di evento(i)
EDSSevent=cell(numel(target_time),1);
MODELevent=cell(numel(target_time));
for k=1:numel(target_time)
    % EDSSevent=TRFsub.classes=='3'; % RISK CLASS
    EDSSevent{k}=Ttime{k}.risk_class=='severe'; % RISK CLASS 
    MODELevent{k}=predictions{k}=='severe';
end
%% creazione tabella x SURVIVAL analysis
Tsurv = cell(numel(target_time),1);Tsurv_complete = cell(numel(target_time),1);
for k=1:numel(target_time)
    pat=unique(TRFtime{k}.Paziente);
% pat=unique(TtrainRF.Paziente);
    Tsurv{k}=table('size',[numel(pat),5],'variabletypes',{'string','logical','double','logical','double'},'variablenames',{'Paziente','event','EventTime','model','ModelTime'});
    patindex=NaN(numel(pat),1);
    for j=1:numel(pat)
        w = TRFtime{k}.Paziente==pat(j);
        Tpat=TRFtime{k}(w,:);
	    patindex(j)=find(w,1,'first');
	    trueeventline=find(EDSSevent{k}(w),1,'first');
	    modeleventline=find(MODELevent{k}(w),1,'first');
	    truedateevent=Tpat.deltaT(trueeventline);
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
	    Tsurv{k}(j,:)={pat(j) trueevent truedateevent modelevent modeldateevent};
    end
    Tsurv{k}.EventTimeY=years(Tsurv{k}.EventTime);
    Tsurv{k}.ModelTimeY=years(Tsurv{k}.ModelTime);
    vars=setdiff(TRFtime{k}.Properties.VariableNames,"Paziente");
    Tsurv_complete{k}=[Tsurv{k},TRFtime{k}(patindex,vars)];
end

%% KAPLAN MEIER analysis

%ht=tiledlayout("flow");

% for k=1:numel(target_time)
    categoricalgroup = categorical([repmat("true data",[height(Tsurv{k}),1]);repmat("model",[height(Tsurv{k}),1])]);

    t = [Tsurv{k}.EventTimeY;Tsurv{k}.ModelTimeY];

    censored=[~Tsurv{k}.event;~Tsurv{k}.model];

    cmap=bone;
    f=figure(5+k);clf;
   % ax(k)=nexttile;
    KaplanMeier(categoricalgroup,years(t),censored,true,cmap);
    xlabel('Years');
    ax = gca;
    ax.FontSize = 12;

    ax2 = gca;
    ax2=axes('Position',ax.Position,'XAxisLocation','top','YAxisLocation','right','Color','none','XLim',ax.XLim,'XTick',ax.XTick);
    ax2.YAxis.Visible='off';ax2.XAxis.Color='b';
    num=str2double(ax.XTickLabel);
    ev=zeros(1,numel(num));
    logicalgroup=dummyvar(categoricalgroup);
    logicalgroup=logical(logicalgroup(:,1));
    for j=2:numel(num)-1;ev(j)=nnz(t(logicalgroup)<=years(num(j)));end    
    surv=nnz(logicalgroup)-ev;
    ax2.XTickLabel=cat(1,cellstr(num2str(surv(1:end-1)')),'n at risk',' ');
    %title(sprintf('%d year of baseline (T=%d yrs)',tstar, target_time(k)))
     title(sprintf('Kaplan Meier curves (T=%d yrs)', target_time(k)))
% title(sprintf('%d years of baseline - prediction at %d year',tstar, target_time))
% [truewblfit,trueexpfit]=singleKaplanMeier(ht,categoricalgroup=="true data",years(t),censored,[cmap(1,:);cmap(1,:).^2;cmap(1,:).^(0.5)]);
% [modelwblfit,modelexpfit]=singleKaplanMeier(ht,categoricalgroup=="model",years(t),censored,[cmap(2,:);cmap(2,:).^2;cmap(2,:).^(0.5)]);
% ax = gca; 
    ax2.FontSize = 12; 
    figure(8+k);clf;
    logrank([Tsurv{k}.EventTime,~Tsurv{k}.event],[Tsurv{k}.ModelTime,~Tsurv{k}.model]);
% end

% export_fig kmBW_3 -png -jpg -pdf -transparent

