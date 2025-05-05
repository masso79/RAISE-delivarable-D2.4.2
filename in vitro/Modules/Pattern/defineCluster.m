function [NBsetChoose, SsetChoose] = defineCluster(SSet, NBset, NetReshaped, Aset, ElActive)
%Seleziono l'indice setMag del set che ha similarità interna maggiore
cnt = 1;
while sum(SSet) ~=0
    clear delete
    [value, setMag] = max(SSet);

    % % Elimino i set che contengono almeno un NB appartenente al set
    % identificato

    chooseNB = NBset{setMag};
    NBsetChoose{cnt} = chooseNB;
    SsetChoose(cnt,:) = [SSet(setMag) setMag];
    
    count = 1;
    for i = 1:length(NBset)
        if ~isempty(intersect(NBset{i},chooseNB))
            delete(count) = i;
            count = count + 1;
        end
    end

    NBset(delete) = {0};
    SSet(delete) = 0;

    % Calcolare similarità interna dei NB nel set identificato
    tmpInt = [];
    SSInt = [];
    for i = 1:length(chooseNB)  
        tmpInt = abs(NetReshaped(:,:,chooseNB(i)) - Aset(:,:,setMag));
        check = tmpInt<= 50;
        SSInt(i) = (1/(ElActive*(ElActive-1)))*sum(sum(check));
    end

    valueCheck = min(SSInt);

    % Calcolo il valore di similarità delle matrici dei delay dei NB esterni al
    % set selezionato con la matrice del delay medio del set
    tmpExt = [];
    SSExt = [];
    NBCheck = setdiff(find(SSet), chooseNB);
    for i = 1:length(NBCheck)  
        tmpExt(:,:) = abs(NetReshaped(:,:,NBCheck(i)) - Aset(:,:,setMag));
        check = tmpExt <= 50;
        SSExt(NBCheck(i)) = (1/(ElActive*(ElActive-1)))*sum(sum(check));
    end

    NBdelete = find(SSExt > valueCheck);


    % Trova i set che contengono i NB da eliminare
    count = 1;
    clear delete
    for i = 1:length(NBset)
        if ~isempty(intersect(NBset{i},NBdelete))
            delete(count) = i;
            count = count + 1;
        end
    end
    
    exist delete var
    if ans ~=0
        NBset(delete) = {0};
        SSet(delete) = 0;
    end

    cnt = cnt +1;
    NBset(setMag) = {0};
    SSet(setMag) = 0;
end
end



