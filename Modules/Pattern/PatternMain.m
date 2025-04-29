clear all
clc

fs = 10000;
[file, path] = uigetfile ('.mat', 'Select the NetworkBurst folder');
load(fullfile(path, file));
for i = 1:length(netBurstsPattern)
    el(i) = size(netBurstsPattern{i},1);
end

ElActive = max(el);
mkdir(path,'PatternAnalysis\ActivationMatrix');
NetTot = [];
%% Crea matrici asimmetrice A(i,j) per ogni NB 
% Calcolo del delay tra tutti gli elettrodi per ogni Network Burst
% NetTot è una matrice tridimensionale asimmetrica: sulle righe
% l'elettrodo, anche nelle colonne e i rispettivi ritardi, la terza
% dimensionalità identifica il numero di network burst identificati
for k = 1: length(netBurstsPattern)
    NetDelay = NaN*ones(120,120);
    for i = 1:length(netBurstsPattern{k})
        for t = i+1:size(netBurstsPattern{k},1)
            NetDelay(netBurstsPattern{k}(i,1), netBurstsPattern{k}(t,1)) = (netBurstsPattern{k}(i,2)- netBurstsPattern{k}(t,2))/fs*1000; %ms
            NetDelay(netBurstsPattern{k}(t,1), netBurstsPattern{k}(i,1)) = - NetDelay(netBurstsPattern{k}(i,1), netBurstsPattern{k}(t,1)); 
        end
    end
    NetDelay(logical(eye(120))) = zeros(1,120);
    NetTot = cat(3, NetTot, NetDelay);
    save(strcat(fullfile(path,'PatternAnalysis\ActivationMatrix'),'\DelayMatrix_NB=',string(k)),'NetDelay');
end

%% Calcolare similarità per ogni NB
% La similarità tra network burst viene effettuato confrontando ogni
% ritardo dell'iesimo elettrodo con il i'esimo elettrodo in tutti i network
% burst.
% Se in modulo la differenza dei ritardi dell'kesimo elettrodo con se
% stesso al iesimo e jesimo network bust presentano un delay inferiore a 50
% ms si conta 1, altrimenti 0. 
% I valori vengono sommati e divisi rispetto a numero di
% elettrodi attivi moltiplicati per il numero di elettrodi attivi -1 e
% inserito in una matrice di similarità quadrata di dimensione del numero
% dei network burst S(i,j)

SimTot = [];
mkdir(path,'PatternAnalysis\SimilarityMatrix');
for i = 1:size(NetTot,3)
    for j = i+1:size(NetTot,3)
        tmp = abs(NetTot(:,:,i)-NetTot(:,:,j));
        tmp(logical(eye(120))) = 100.*ones(1,120);
        check = tmp <= 50;
        S(i,j) = (1/(ElActive*(ElActive-1)))*sum(sum(check));
        S(j,i) = S(i,j);
    end
    save(strcat(fullfile(path,'PatternAnalysis\SimilarityMatrix'),'\SimilarityMatrix'),'S');
end

mkdir(strcat(fullfile(path,'PatternAnalysis','\Figures')));
similarity = figure();
imagesc(S);
ax = gca;
axis(ax,'off');
xlabel(ax,'Burst #')
ylabel(ax,'Burst #')
title('Similarity')
axis equal
colorbar
ax.XLabel.Visible = 'on';
ax.YLabel.Visible = 'on';
savefig(similarity, strcat(fullfile(path,'PatternAnalysis\Figures'),'\SimilarityMatrix'))
close(similarity)

%% Riorganizza la matrice S con il metodo agglomerato del dendrogramma
% Sreshaped è la matrice S riorganizzata secondo il metodo standard
% agglomerative dendrogram OrderedS
% SArray = squareform(S);
meanDistFull_dist = pdist(S,'euclidean');
tree = linkage(meanDistFull_dist,'average');
%tree = linkage(SArray,'euclidean');
[a,b,idx] = dendrogram(tree,length(S));

Sreshaped = S(idx,idx);

OrderedS = figure();
imagesc(Sreshaped);
ax = gca;
axis(ax,'off');
xlabel(ax,'Burst #')
ylabel(ax,'Burst #')
title('Similarity')
axis equal
colorbar
ax.XLabel.Visible = 'on';
ax.YLabel.Visible = 'on';
savefig(OrderedS, strcat(fullfile(path,'PatternAnalysis\Figures'),'\OrderedMatrix'))
close(OrderedS)

%% Creare i possibili set definiti come matrici quadrate di lunghezza sqrt(numero di NB)
% Per identificare i reali set:
% Riorganizzo la matrice dei delay A(i,j,k) in modo da avere le matrici
% rispettive ai NB (lungo k) ordinate secondo il metodo precedentemente
% individuato NetReshaped.

NetReshaped = [];
for k = 1:length(idx)
    NetReshaped(:,:,k) = NetTot(:,:,idx(k));
end

% Per ogni possibile set costituito da sqrt(numero di NB) calcolo la 
% matrice dei delay medio Aset e salvo quali network burst appartengono
% al set 
 
Aset = [];
%% Versione corretta
% numberNB = floor(sqrt(length(netBursts)));

%% Versione con imposizione del sqrt(valore minimo del numero di NB tra le configurazioni)
numberNB = 7;

%%
for k = 0:size(NetReshaped,3)-numberNB
    Aset(:,:,k+1) = mean(NetReshaped(:,:,[k+1:k+numberNB]),3);
    NBset{k+1}= [k+1:k+numberNB];
end


%% Calcolare la similarità all'interno di ogni set identificato
% Per ogni set calcolo il valore di similarità interno: prendo l i'esimo
% network burst che appartiene al set e calcolo la differenza tra la
% matrice A rispettiva e la media Aset. Se i valori sono inferiori a 50
% considero 1 0 altrimenti. Sommo tutto e divido per il numero NB
% appartenenti al set.
% SS è il valore di similarità interna di ogni NB nel set 
SSet = [];
SS = [];
mkdir(path,'PatternAnalysis\SimilarityMatrixSet');
count = 1;
for i = 1:size(Aset,3)
    for j = 0:numberNB-1
        start = count+j;
        tmp = abs(NetReshaped(:,:,start)-Aset(:,:,i));
        check = tmp <= 50;
        SS(j+1) = (1/(ElActive*(ElActive-1)))*sum(sum(check));
    end
    SSet(i) = mean(SS);
    count = count+1;
end
save(strcat(path,'PatternAnalysis\SimilarityMatrixSet\SSet'),'SSet');

%% Inserire la funzione per: 
% identificare il set con similarità maggiore
% Rimuovere Set che contengono NB del Set scelto
% Rimuovere i set che hanno un NB con un valore di similarità interno al
% set selezionato superiore al valore minimo di similarità dei NB interni
% al set.

[NBsetChoose, SsetChoose] = defineCluster(SSet, NBset, NetReshaped, Aset, ElActive);
mkdir(path,'PatternAnalysis\Modules')
numModules = length(NBsetChoose);
save(strcat(path,'PatternAnalysis\Modules\Analysis'),'NBsetChoose','SsetChoose','numModules')
close all

%% Identificare il set con similarità maggiore
% Seleziono l'indice setMag del set che ha similarità interna maggiore
% [value, setMag] = max(SSet);

% %% Rimuovi Set che contengono i NB del set con similarità maggiore
% % Elimino i set che contengono almeno un NB appartenente al set
% % identificato
% chooseNB = NBset{setMag};
% count = 1;
% for i = 1:length(NBset)
%     if ~isempty(intersect(NBset{i},chooseNB))
%         delete(count) = i;
%         count = count + 1;
%     end
% end
% 
% NBset(delete) = {0};
% SSet(delete) = 0;
% 
% %% Rimuovo i set che hanno un NB che ha un valore di similarità interno al set selezionato superiore al valore minimo di similarità dei NB interni al set.
% % Calcolare similarità interna dei NB nel set identificato
% tmpInt = [];
% SSInt = [];
% for i = 1:length(chooseNB)  
%     tmpInt = abs(NetReshaped(:,:,chooseNB(i)) - Aset(:,:,setMag));
%     check = tmpInt<= 50;
%     SSInt(i) = (1/(ElActive*(ElActive-1)))*sum(sum(check));
% end
% 
% valueCheck = min(SSInt);
% 
% % Calcolo il valore di similarità delle matrici dei delay dei NB esterni al
% % set selezionato con la matrice del delay medio del set
% tmpExt = [];
% SSExt = []; 
% NBCheck = setdiff([1:length(netBursts)], chooseNB);
% for i = 1:length(NBCheck)  
%     tmpExt(:,:) = abs(NetReshaped(:,:,NBCheck(i)) - Aset(:,:,setMag));
%     check = tmpExt <= 50;
%     SSExt(NBCheck(i)) = (1/(ElActive*(ElActive-1)))*sum(sum(check));
% end
% 
% 
% NBdelete = find(SSExt > valueCheck);
% 
% % Trova i set che contengono i NB da eliminare
% 
% count = 1;
% clear delete
% for i = 1:length(NBset)
%     if ~isempty(intersect(NBset{i},NBdelete))
%         delete(count) = i;
%         count = count + 1;
%     end
% end
% NBset(delete) = {0};
% SSet(delete) = 0;
% 
% 
% %Salva per il set identificato il valore di similarità interno di ogni NB
% % i NB appartenenti al set e i valori di NetReshaped dei rispettivi NB e la
% % matrice dei delay medio del set