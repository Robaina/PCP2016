function BIO=getBiologyNoData(GEM,fctable)

%Extract Coupled Groups and get cellular localization of reactions and coupled groups, as well as
%implicated metabolic subsystems

for b=1:2,
    ctype=[1,4];
    Cgroups={};
    n=1;
    for i=1:size(fctable,1),
        if length(find(fctable(:,i)==ctype(b)))>1,
            FC=unique([i;find(fctable(:,i)==ctype(b))]);
            Cgroups{n,1}=FC;
            n=n+1;
        end
    end
    %eliminate repeated groups
    CgroupsMat=zeros(length(Cgroups),size(fctable,1));
    for w=1:length(Cgroups),
        CgroupsMat(w,Cgroups{w})=1;
    end
    [~,Uidx]=unique(CgroupsMat,'rows');
    Cgroups=Cgroups(Uidx);
  
    %map GEM reactions to localizations in the maize leaf model
    if ~isfield(GEM,'MeRxns'),
        MeRxns=[];BSRxns=[];MeBSRxns=[];LeafRxns=[];
        for i=1:length(GEM.rxns),
            MetLoc=GEM.mets(find(GEM.S(:,i)~=0));
            for k=1:length(MetLoc),
                L(1)=strfind(MetLoc{k},'[');
                L(2)=strfind(MetLoc{k},']');
                MetLoc{k}=MetLoc{k}(L(1)+1:L(2)-1);
            end
            if sum(ismember(MetLoc,GEM.Mecompartments(:,2)))>0 && sum(ismember(MetLoc,GEM.BScompartments(:,2)))==0,
                MeRxns=[MeRxns,i];
            elseif sum(ismember(MetLoc,GEM.BScompartments(:,2)))>0 && sum(ismember(MetLoc,GEM.Mecompartments(:,2)))==0,
                BSRxns=[BSRxns,i];
            elseif sum(ismember(MetLoc,GEM.Mecompartments(:,2)))>0 && sum(ismember(MetLoc,GEM.BScompartments(:,2)))>0,
                MeBSRxns=[MeBSRxns,i];
            elseif sum(ismember(MetLoc,'L'))>0 && sum(ismember(MetLoc,GEM.BScompartments(:,2)))==0 && sum(ismember(MetLoc,GEM.Mecompartments(:,2)))==0,
                LeafRxns=[LeafRxns,i];
            end
        end

        GEM.MeRxns=unique(MeRxns)';
        GEM.BSRxns=unique(BSRxns)';
        GEM.MeBSRxns=unique(MeBSRxns)';
        GEM.LeafRxns=unique(LeafRxns)';
    end

    %Localize reactions in FCgroups and FCgroupsData:
       %Leaf=1
       %Me=2
       %BS=3
       %MeBS=4
       %Undetermined=0
    TotsubSys={};TotCoupledsubSysMeBS={};TotCoupledsubSysMe={};TotCoupledsubSysBS={};TotCoupledsubSysLeaf={};
    for l=1:length(Cgroups),
        GroupCard(l,1)=length(Cgroups{l});
        LocCgroups{l,1}=zeros(length(Cgroups{l}),1);
        LocCgroups{l,1}(ismember(Cgroups{l},GEM.LeafRxns))=1;
        LocCgroups{l,1}(ismember(Cgroups{l},GEM.MeRxns))=2;
        LocCgroups{l,1}(ismember(Cgroups{l},GEM.BSRxns))=3;
        LocCgroups{l,1}(ismember(Cgroups{l},GEM.MeBSRxns))=4;
        if ismember(4,LocCgroups{l,1})||(ismember(2,LocCgroups{l,1}) && ismember(3,LocCgroups{l,1})),
           CompartDist{l,1}='Me-BS';
        elseif ismember(2,LocCgroups{l,1}) && ~ismember(3,LocCgroups{l,1}) && ~ismember(4,LocCgroups{l,1}), 
           CompartDist{l,1}='Me';
        elseif ismember(3,LocCgroups{l,1}) && ~ismember(2,LocCgroups{l,1}) && ~ismember(4,LocCgroups{l,1}), 
            CompartDist{l,1}='BS';
        elseif ismember(1,LocCgroups{l,1}) && ~ismember(2,LocCgroups{l,1}) && ~ismember(3,LocCgroups{l,1}) && ~ismember(4,LocCgroups{l,1}),
           CompartDist{l,1}='Leaf';
        elseif sum(LocCgroups{l,1})==0,   
            CompartDist{l,1}='Undetermined';
        end
        
        %get subsystems in group and number of reactions per subsystem
        if ~isfield(GEM,'UniquesubSystems'),
            GEM.uniquesubSystems=unique(GEM.subSystems);
            for h=1:length(GEM.uniquesubSystems),
                GEM.uniquesubSystems{h,2}=length(find(strcmp(GEM.subSystems,GEM.uniquesubSystems{h})));
            end
        end
        S1=GEM.subSystems(Cgroups{l},1);
        S1(strcmp(S1,''))=[];
        subSystems{l,1}=unique(S1);
        TotsubSys=[TotsubSys;S1];
        if strcmp(CompartDist{l},'Me-BS'),
            S22=GEM.subSystems(Cgroups{l},1);
            S22(strcmp(S22,''))=[];
            TotCoupledsubSysMeBS=[TotCoupledsubSysMeBS;S22];
            UCoupledsubSysMeBS=unique(TotCoupledsubSysMeBS);
        elseif strcmp(CompartDist{l},'Me'),
            S32=GEM.subSystems(Cgroups{l},1);
            S32(strcmp(S32,''))=[];
            TotCoupledsubSysMe=[TotCoupledsubSysMe;S32];
            UCoupledsubSysMe=unique(TotCoupledsubSysMe);
        elseif strcmp(CompartDist{l},'BS'),
            S42=GEM.subSystems(Cgroups{l},1);
            S42(strcmp(S42,''))=[];
            TotCoupledsubSysBS=[TotCoupledsubSysBS;S42];
            UCoupledsubSysBS=unique(TotCoupledsubSysBS);
        elseif strcmp(CompartDist{l},'Leaf'),
            S52=GEM.subSystems(Cgroups{l},1);
            S52(strcmp(S52,''))=[];
            TotCoupledsubSysLeaf=[TotCoupledsubSysLeaf;S52];
            UCoupledsubSysLeaf=unique(TotCoupledsubSysLeaf);
        end
    end

    %Count reactions in subsystems
    A=[];B=[];C=[];D=[];
    if exist('UCoupledsubSysMe','var') && ~isempty(UCoupledsubSysMe),
        for a=1:length(UCoupledsubSysMe),
            A(a,1)=length(find(strcmp(TotCoupledsubSysMe,UCoupledsubSysMe{a})));
            A(a,2)=length(find(strcmp(TotCoupledsubSysMe,UCoupledsubSysMe{a})))/GEM.uniquesubSystems{find(strcmp(GEM.uniquesubSystems(:,1),UCoupledsubSysMe{a})),2};
        end
        [~,Idxa]=sort(A(:,2),'descend');
        UCoupledsubSysMe(:,1)=UCoupledsubSysMe(Idxa,1);
        UCoupledsubSysMe(:,2:3)=num2cell(A(Idxa,:));
    elseif ~exist('UCoupledsubSysMe','var'),
        UCoupledsubSysMe={};
    end
    if exist('UCoupledsubSysBS','var') && ~isempty(UCoupledsubSysBS),
        for a=1:length(UCoupledsubSysBS),
            B(a,1)=length(find(strcmp(TotCoupledsubSysBS,UCoupledsubSysBS{a})));
            B(a,2)=length(find(strcmp(TotCoupledsubSysBS,UCoupledsubSysBS{a})))/GEM.uniquesubSystems{find(strcmp(GEM.uniquesubSystems(:,1),UCoupledsubSysBS{a})),2};
        end
        [~,Idxb]=sort(B(:,2),'descend');
        UCoupledsubSysBS(:,1)=UCoupledsubSysBS(Idxb,1);
        UCoupledsubSysBS(:,2:3)=num2cell(B(Idxb,:));
    elseif ~exist('UCoupledsubSysBS','var'),
       UCoupledsubSysBS={};
    end
    if exist('UCoupledsubSysMeBS','var') && ~isempty(UCoupledsubSysMeBS),
        for a=1:length(UCoupledsubSysMeBS),
            C(a,1)=length(find(strcmp(TotCoupledsubSysMeBS,UCoupledsubSysMeBS{a})));
            C(a,2)=length(find(strcmp(TotCoupledsubSysMeBS,UCoupledsubSysMeBS{a})))/GEM.uniquesubSystems{find(strcmp(GEM.uniquesubSystems(:,1),UCoupledsubSysMeBS{a})),2};
        end
        [~,Idxc]=sort(C(:,2),'descend');
        UCoupledsubSysMeBS(:,1)=UCoupledsubSysMeBS(Idxc,1);
        UCoupledsubSysMeBS(:,2:3)=num2cell(C(Idxc,:));
    elseif ~exist('UCoupledsubSysMeBS','var'),
        UCoupledsubSysMeBS={};
    end
     if exist('UCoupledsubSysLeaf','var') && ~isempty(UCoupledsubSysLeaf),
        for a=1:length(UCoupledsubSysLeaf),
            D(a,1)=length(find(strcmp(TotCoupledsubSysLeaf,UCoupledsubSysLeaf{a})));
            D(a,2)=length(find(strcmp(TotCoupledsubSysLeaf,UCoupledsubSysLeaf{a})))/GEM.uniquesubSystems{find(strcmp(GEM.uniquesubSystems(:,1),UCoupledsubSysLeaf{a})),2};
        end
        [~,Idxd]=sort(D(:,2),'descend');
        UCoupledsubSysLeaf(:,1)=UCoupledsubSysLeaf(Idxd,1);
        UCoupledsubSysLeaf(:,2:3)=num2cell(D(Idxd,:));
    elseif ~exist('UCoupledsubSysLeaf','var'),
        UCoupledsubSysLeaf={};
     end
     
    %Analize coupled groups
    for i=1:length(Cgroups(:,1)),
        CgroupsRxns{i,1}=GEM.rxnNames(Cgroups{i,1});
        CgroupsSubSys{i,1}=GEM.subSystems(Cgroups{i,1},1);
        CgroupsMechanisms{i,1}=GEM.rxnMechanisms(Cgroups{i,1});
    end
    for i=1:length(Cgroups),
        CGidx(i,1)=Cgroups{i}(1);
        CGidx(i,2)=length(Cgroups{i});
    end
    
    %Count number of coupled reaction pairs per cellular localization
    
    %Localize groups and leading reaction in M/BS directionally coupled groups
        MeGroups=(strcmpi(CompartDist,'Me'));
        BSGroups=(strcmpi(CompartDist,'BS'));
        MeBSGroups=(strcmpi(CompartDist,'Me-BS'));
        LeafGroups=(strcmpi(CompartDist,'Leaf'));
        OtherGroups=setdiff(1:length(Cgroups),unique([find(MeGroups==1);find(BSGroups==1);find(MeBSGroups==1);find(LeafGroups==1)]));
    if ctype(b)==1,
        BIO.MeFCgroupsMedRxn=median(CGidx(MeGroups,2));
        BIO.BSFCgroupsMedRxn=median(CGidx(BSGroups,2));
        BIO.MeBSFCgroupsMedRxn=median(CGidx(MeBSGroups,2));
        BIO.LeafFCgroupsMedRxn=median(CGidx(LeafGroups,2));
        BIO.OtherFCgroupsMedRxn=median(CGidx(OtherGroups,2));
        BIO.MFCgroups=CGidx(MeGroups,1);
        BIO.BSFCgroups=CGidx(BSGroups,1);
        BIO.MBSFCgroups=CGidx(MeBSGroups,1);
        BIO.LeafFCgroups=CGidx(LeafGroups,1);
        BIO.OtherFCgroups=CGidx(OtherGroups,1);
    elseif ctype(b)==4,
        BIO.MeDCgroupsMedRxn=median(CGidx(MeGroups,2));
        BIO.BSDCgroupsMedRxn=median(CGidx(BSGroups,2));
        BIO.MeBSDCgroupsMedRxn=median(CGidx(MeBSGroups,2));
        BIO.LeafDCgroupsMedRxn=median(CGidx(LeafGroups,2));
        BIO.OtherDCgroupsMedRxn=median(CGidx(OtherGroups,2));
        BIO.MDCgroups=CGidx(MeGroups,1);
        BIO.BSDCgroups=CGidx(BSGroups,1);
        BIO.MBSDCgroups=CGidx(MeBSGroups,1);
        BIO.LeafDCgroups=CGidx(LeafGroups,1);
        BIO.OtherDCgroups=CGidx(OtherGroups,1);
        BIO.MLeading=find(ismember(BIO.MBSDCgroups,GEM.MeRxns)==1);
        BIO.BSLeading=find(ismember(BIO.MBSDCgroups,GEM.BSRxns)==1);
        BIO.MBSLeading=find(ismember(BIO.MBSDCgroups,GEM.MeBSRxns)==1);
        BIO.OtherLeading=setdiff(BIO.MBSDCgroups,[BIO.MLeading;BIO.BSLeading;BIO.MBSLeading]);
    end
    Tablecell=[num2cell(CGidx(:,1)),Cgroups,CompartDist,CgroupsRxns,CgroupsSubSys,CgroupsMechanisms];
    if ctype(b)==1,
        names2{1}='FC Group Idx';names2{2}='FC Group';names2{3}='Cell Type';names2{4}='Rxn. Names';names2{6}='Rxn. Mechanisms';names2{5}='SubSystems';
        FCtablecell=[names2;Tablecell];
        clearvars -except FCtablecell GEM fctable alpha PcTot CorSeq rhoCutoff permutations DataCorTot BIO MeanGroupCorFC MeanGroupCorDC
    elseif ctype(b)==4,

         names2{1}='DC Group Idx';names2{2}='DC Group';names2{3}='Cell Type';names2{4}='Rxn. Names';names2{6}='Rxn. Mechanisms';names2{5}='SubSystems';
         DCtablecell=[names2;Tablecell];

    end
end

%Obtain total FC, DC and uncoupled pairs of reactions pair cell type
fctable(tril(fctable==fctable))=-1;
[FCpairs(:,1),FCpairs(:,2)]=find(fctable==1);
[DCpairs(:,1),DCpairs(:,2)]=find(fctable==4);
[UnCpairs(:,1),UnCpairs(:,2)]=find(fctable==0);

BIO.FCtablecell=FCtablecell;
BIO.DCtablecell=DCtablecell;
BIO.UnPairs=UnCpairs;
BIO.FCPairs=FCpairs;
BIO.DCPairs=DCpairs;

end