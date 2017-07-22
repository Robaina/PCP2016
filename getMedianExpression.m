function Sol=getMedianExpression(GEM,fctable,Data,permutations,FigureName)

%get total fully (partially) and directionally coupled reactions
fctable(tril(fctable==fctable))=-1;
[FCrxns(:,1),FCrxns(:,2)]=find(fctable==1);
[DCrxns(:,1),DCrxns(:,2)]=find(fctable==3);
FCrxns=unique([FCrxns;FCrxns]);
DCrxns=unique([DCrxns;DCrxns]);
UnCrxns=unique(setdiff(1:length(GEM.rxns),[FCrxns;DCrxns]));
Sol.FCrxns=FCrxns;
Sol.DCrxns=DCrxns;
Sol.UnCrxns=UnCrxns;
%get median and mad values of each group and generate samples from the
%total set of reactions
Data(isinf(Data))=0;
DataFC=Data(FCrxns,:);DataDC=Data(DCrxns,:);DataUnC=Data(UnCrxns,:);
Sol.DataFC=DataFC(find(DataFC(:,1)~=0));
Sol.DataDC=DataDC(find(DataDC(:,1)~=0));
Sol.DataUnC=DataUnC(find(DataUnC(:,1)~=0));
for i=1:size(Data,2),
    Totmedmad(i,1)=median(Data(find(Data(:,i)~=0),i));
    Totmedmad(i,2)=mad(Data(find(Data(:,i)~=0),i));
    FCmedmad(i,1)=median(DataFC(find(DataFC(:,i)~=0),i));
    FCmedmad(i,2)=mad(DataFC(find(DataFC(:,i)~=0),i));
    DCmedmad(i,1)=median(DataDC(find(DataDC(:,i)~=0),i));
    DCmedmad(i,2)=mad(DataDC(find(DataDC(:,i)~=0),i));
    UnCmedmad(i,1)=median(DataUnC(find(DataUnC(:,i)~=0),i));
    UnCmedmad(i,2)=mad(DataUnC(find(DataUnC(:,i)~=0),i));
end

for i=1:permutations,
    P1=randsample(1:length(FCrxns),length(UnCrxns),'true');
    DataP1=DataFC(P1,:);
    for j=1:size(Data,2),
        TotsamplesFCmed(i,j)=median(DataP1(find(DataP1(:,j)~=0),j));
        TotsamplesFCmad(i,j)=mad(DataP1(find(DataP1(:,j)~=0),j));
    end
    
    P2=randsample(1:length(DCrxns),length(UnCrxns),'true');
    DataP2=DataDC(P2,:);
    for j=1:size(Data,2),
        TotsamplesDCmed(i,j)=median(DataP2(find(DataP2(:,j)~=0),j));
        TotsamplesDCmad(i,j)=mad(DataP2(find(DataP2(:,j)~=0),j));
    end
end

FCmedPvalue=zeros(1,size(Data,2));
FCmadPvalue=FCmedPvalue;
DCmedPvalue=FCmedPvalue;
DCmadPvalue=FCmedPvalue;

%get proportion of samples with greater mean and mad value (p-value
%estimation)

for i=1:size(Data,2),
    FCmedPvalue(1,i)=length(find(TotsamplesFCmed(:,i)>=UnCmedmad(i,1)));
    FCmadPvalue(1,i)=length(find(TotsamplesFCmad(:,i)>=UnCmedmad(i,2)));
    DCmedPvalue(1,i)=length(find(TotsamplesDCmed(:,i)>=UnCmedmad(i,1)));
    DCmadPvalue(1,i)=length(find(TotsamplesDCmad(:,i)>=UnCmedmad(i,2)));
end

Sol.FCmedPvalue=FCmedPvalue/permutations;
Sol.FCmadPvalue=FCmadPvalue/permutations;
Sol.DCmedPvalue=DCmedPvalue/permutations;
Sol.DCmadPvalue=DCmadPvalue/permutations;
Sol.Totmedmad=Totmedmad';
Sol.FCmedmad=FCmedmad';
Sol.DCmedmad=DCmedmad';
Sol.UnCmedmad=UnCmedmad';
Sol.TotsamplesFCmad=TotsamplesFCmad;
Sol.TotsamplesDCmad=TotsamplesDCmad;

h(1)=figure('Color','w');
g=bar([FCmedmad(:,1),DCmedmad(:,1),Totmedmad(:,1)],'grouped');
set(g(3),'FaceColor',[0 0 0])
xlabel('leaf sections','FontSize',16);ylabel('median expression','FontSize',16);
set(gca,'fontsize',15)
h_legend=legend(g,{'FC pairs','DC pairs','UnC pairs'});
legend boxoff 
set(h_legend,'FontSize',12);
% print(h(1),[Directory,FigureName,'\','MedianExp'],'-dtiff')

h(2)=figure('Color','w');
g=bar([FCmedmad(:,2),DCmedmad(:,2),Totmedmad(:,2)],'grouped');
set(g(3),'FaceColor',[0 0 0])
xlabel('leaf sections','FontSize',16);ylabel('mad','FontSize',16);
set(gca,'fontsize',15)
% print(h(2),[Directory,FigureName,'\','MadExp'],'-dtiff')

h(3)=figure('Color','w');
g=bar([Sol.FCmedPvalue',Sol.DCmedPvalue'],'grouped');
% set(g(1),'FaceColor',[0.9 0.0 0.0])
set(g(2),'FaceColor',[0.0590 0.6838 0.7254])
xlabel('leaf sections','FontSize',16);ylabel('p-value(median)','FontSize',16);
set(gca,'fontsize',15)
% print(h(3),[Directory,FigureName,'\','PvalueMedian'],'-dtiff')

h(4)=figure('Color','w');
g=bar([Sol.FCmadPvalue',Sol.DCmadPvalue'],'grouped');
set(g(2),'FaceColor',[0.0590 0.6838 0.7254])
xlabel('leaf sections','FontSize',16);ylabel('p-value(mad)','FontSize',16);
set(gca,'fontsize',15)
% print(h(4),[Directory,FigureName,'\','PvalueMad'],'-dtiff')
savefig(h,[FigureName,'.fig'],'compact')
close all
end
