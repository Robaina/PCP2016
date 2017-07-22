function GEM=getRxnLocation(GEM)

    MeRxns=[];BSRxns=[];MeBSRxns=[];LeafRxns=[];
    BScomp={'c','i','m','p','a','t','v','x'};
    Mecomp={'d','j','n','q','b','u','w','f'};
    for i=1:length(GEM.rxns),
        MetLoc=GEM.mets(find(GEM.S(:,i)~=0));
        for k=1:length(MetLoc),
            L(1)=strfind(MetLoc{k},'[');
            L(2)=strfind(MetLoc{k},']');
            MetLoc{k}=MetLoc{k}(L(1)+1:L(2)-1);
        end
        if sum(ismember(MetLoc,Mecomp))>0 && sum(ismember(MetLoc,BScomp))==0,
            MeRxns=[MeRxns,i];
        elseif sum(ismember(MetLoc,BScomp))>0 && sum(ismember(MetLoc,Mecomp))==0,
            BSRxns=[BSRxns,i];
        elseif sum(ismember(MetLoc,Mecomp))>0 && sum(ismember(MetLoc,BScomp))>0,
            MeBSRxns=[MeBSRxns,i];
        elseif sum(ismember(MetLoc,'L'))>0 && sum(ismember(MetLoc,BScomp))==0 && sum(ismember(MetLoc,Mecomp))==0,
            LeafRxns=[LeafRxns,i];
        end
    end

    GEM.MeRxns=unique(MeRxns)';
    GEM.BSRxns=unique(BSRxns)';
    GEM.MeBSRxns=unique(MeBSRxns)';
    GEM.LeafRxns=unique(LeafRxns)';
    GEM.OtherRxns=setdiff(1:length(GEM.rxns),unique([MeRxns,BSRxns,MeBSRxns,LeafRxns]))';
    GEM.Mecompartments=Mecomp;
    GEM.BScompartments=BScomp;
 end
