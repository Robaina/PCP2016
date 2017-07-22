function Val=getSubSystems(GEM,BIO)

%gets ranked table of subsystems in FC and DC groups per localization
Mfcgroups=find(strcmpi(BIO.FCtablecell(:,3),'Me'));
BSfcgroups=find(strcmpi(BIO.FCtablecell(:,3),'BS'));
MBSfcgroups=find(strcmpi(BIO.FCtablecell(:,3),'Me-BS'));

Mdcgroups=find(strcmpi(BIO.DCtablecell(:,3),'Me'));
BSdcgroups=find(strcmpi(BIO.DCtablecell(:,3),'BS'));
MBSdcgroups=find(strcmpi(BIO.DCtablecell(:,3),'Me-BS'));

%gets ranked table of subsystems of highly correlated FC and DC groups 
Hfcgroups=find(cell2mat(BIO.FCtablecellData(2:end,6))==1)+1;
Hdcgroups=find(cell2mat(BIO.DCtablecellData(2:end,6))==1)+1;

GEMsubSys=unique(GEM.subSystems);
for i=1:length(GEMsubSys),
    N1(i,1)=length(find(strcmpi(GEM.subSystems,GEMsubSys{i})));
end

MSubSysFC=[];BSSubSysFC=[];MBSSubSysFC=[];
MSubSysDC=[];BSSubSysDC=[];MBSSubSysDC=[];
HSubSysFC=[];HSubSysDC=[];

%SubSystems of groups per localization
%FC
for i=1:length(Mfcgroups),
    MSubSysFC=[MSubSysFC;BIO.FCtablecell{Mfcgroups(i),5}(:,1)];
end
MSubSysUFC=unique(MSubSysFC);
for i=1:length(MSubSysUFC),
    N2(i,1)=length(find(strcmpi(MSubSysFC,MSubSysUFC{i})));
    N2(i,2)=N2(i,1)/N1(find(strcmpi(GEMsubSys,MSubSysUFC{i})));
end
[~,idx2]=sort(N2(:,2),'descend');
N2=N2(idx2,:);MSubSysUFC=MSubSysUFC(idx2);

for i=1:length(BSfcgroups),
    BSSubSysFC=[BSSubSysFC;BIO.FCtablecell{BSfcgroups(i),5}(:,1)];
end
BSSubSysUFC=unique(BSSubSysFC);
for i=1:length(BSSubSysUFC),
    N3(i,1)=length(find(strcmpi(BSSubSysFC,BSSubSysUFC{i})));
    N3(i,2)=N3(i,1)/N1(find(strcmpi(GEMsubSys,BSSubSysUFC{i})));
end
[~,idx3]=sort(N3(:,2),'descend');
N3=N3(idx3,:);BSSubSysUFC=BSSubSysUFC(idx3);

for i=1:length(MBSfcgroups),
    MBSSubSysFC=[MBSSubSysFC;BIO.FCtablecell{MBSfcgroups(i),5}(:,1)];
end
MBSSubSysUFC=unique(MBSSubSysFC);
try
    for i=1:length(MBSSubSysUFC),
        N4(i,1)=length(find(strcmpi(MBSSubSysFC,MBSSubSysUFC{i})));
        N4(i,2)=N4(i,1)/N1(find(strcmpi(GEMsubSys,MBSSubSysUFC{i})));
    end
    [~,idx4]=sort(N4(:,2),'descend');
    N4=N4(idx4,:);MBSSubSysUFC=MBSSubSysUFC(idx4);
catch
    MBSSubSysUFC=[];
end

%DC
for i=1:length(Mdcgroups),
    MSubSysDC=[MSubSysDC;BIO.DCtablecell{Mdcgroups(i),5}(:,1)];
end
MSubSysUDC=unique(MSubSysDC);
for i=1:length(MSubSysUDC),
    N5(i,1)=length(find(strcmpi(MSubSysDC,MSubSysUDC{i})));
    N5(i,2)=N5(i,1)/N1(find(strcmpi(GEMsubSys,MSubSysUDC{i})));
end
[~,idx5]=sort(N5(:,2),'descend');
N5=N5(idx5,:);MSubSysUDC=MSubSysUDC(idx5);

for i=1:length(BSdcgroups),
    BSSubSysDC=[BSSubSysDC;BIO.DCtablecell{BSdcgroups(i),5}(:,1)];
end
BSSubSysUDC=unique(BSSubSysDC);
for i=1:length(BSSubSysUDC),
    N6(i,1)=length(find(strcmpi(BSSubSysDC,BSSubSysUDC{i})));
    N6(i,2)=N6(i,1)/N1(find(strcmpi(GEMsubSys,BSSubSysUDC{i})));
end
[~,idx6]=sort(N6(:,2),'descend');
N6=N6(idx6,:);BSSubSysUDC=BSSubSysUDC(idx6);

for i=1:length(MBSdcgroups),
    MBSSubSysDC=[MBSSubSysDC;BIO.DCtablecell{MBSdcgroups(i),5}(:,1)];
end
MBSSubSysUDC=unique(MBSSubSysDC);
try
    for i=1:length(MBSSubSysUFC),
        N7(i,1)=length(find(strcmpi(MBSSubSysFC,MBSSubSysUFC{i})));
        N7(i,2)=N7(i,1)/N1(find(strcmpi(GEMsubSys,MBSSubSysUDC{i})));
    end
    [~,idx7]=sort(N7(:,2),'descend');
    N7=N7(idx7,:);MBSSubSysUDC=MBSSubSysUDC(idx7);
catch
    MBSSubSysUDC=[];
end

%SubSystems of highly correlated FC and DC groups
%FC
for i=1:length(Hfcgroups),
    HSubSysFC=[HSubSysFC;BIO.FCtablecellData{Hfcgroups(i),9}(:,1)];
end
HSubSysUFC=unique(HSubSysFC);
for i=1:length(HSubSysUFC),
    N8(i,1)=length(find(strcmpi(HSubSysFC,HSubSysUFC{i})));
    N8(i,2)=N8(i,1)/N1(find(strcmpi(GEMsubSys,HSubSysUFC{i})));
end
[~,idx8]=sort(N8(:,2),'descend');
N8=N8(idx8,:);HSubSysUFC=HSubSysUFC(idx8);

%DC
for i=1:length(Hdcgroups),
    HSubSysDC=[HSubSysDC;BIO.DCtablecellData{Hdcgroups(i),9}(:,1)];
end
HSubSysUDC=unique(HSubSysDC);
for i=1:length(HSubSysUDC),
    N9(i,1)=length(find(strcmpi(HSubSysDC,HSubSysUDC{i})));
    N9(i,2)=N9(i,1)/N1(find(strcmpi(GEMsubSys,HSubSysUDC{i})));
end
[~,idx9]=sort(N9(:,2),'descend');
N9=N9(idx9,:);HSubSysUDC=HSubSysUDC(idx9);

names={'M-subSystem','#Rxns','Proportion','BS-subSystem','#Rxns','Proportion','MBS-subSystem','#Rxns','Proportion','HC-subSystem','#Rxns','Proportion'};
Val.FCMsubsystems=[names(1:3);MSubSysUFC,num2cell(N2)];
Val.FCBSsubsystems=[names(4:6);BSSubSysUFC,num2cell(N3)];
if ~isempty(MBSSubSysUFC)
Val.FCMBSsubsystems=[names(7:9);MBSSubSysUFC,num2cell(N4)];
end
Val.DCMsubsystems=[names(1:3);MSubSysUDC,num2cell(N5)];
Val.DCBSsubsystems=[names(4:6);BSSubSysUDC,num2cell(N6)];
if ~isempty(MBSSubSysUDC)
Val.DCMBSsubsystems=[names(7:9);MBSSubSysUDC,num2cell(N7)];
end
Val.HFCsubsytems=[names(10:12);HSubSysUFC,num2cell(N8)];
Val.HDCsubsytems=[names(10:12);HSubSysUDC,num2cell(N9)];

[~,Idx]=sort(N1,'descend');
N1=N1(Idx);GEMsubSys=GEMsubSys(Idx);
Val.GEMsubSys=[{'subSystem','#Rxns'};GEMsubSys,num2cell(N1)];

end
