function Val=MaizeAnalysis(GEM,Data,fctable,alpha,rhoCutoff,permutations,KSalphavalue,printFigures,FigureName,biology)

%Main function: GEM: TS4red model Data: already mapped gene data
%(mapgene2rxn) or the privided Wang(Pick)Data matrix 
%fctable: coupling table returned by F2C2, provided in data files 
%alpha: the significance level for the evaluation (0.05 in this study) 
%rhoCutoff: a cutoff value to obtain highly correlated groups, it is meant to be a
%percentile in the correlation distribution(0.9 in this study)
%permutations: number of random samples to extract (10^4 in this study)

%An example would be:

% Leaf1=MaizeFC(TS4red,LeafData1,fctable,0.05,0.9,1e4);

%USAGE with the truncated model (TS4redTruncated), only the data files for
%Leaf Data 1 are provided (in this study the correlation in the truncated
%model was not evaluated, only the coupled groups). In this case, the set
%up would be:

% LeafTrunc=MaizeFC(TS4redTruncated,LeafData1Truncated,fctableTruncated,0.05,0.9,1e4)

%**************************************************************************
%**************************************************************************

%Take only reactions with full spatial data coverage and identify reaction
%pairs with the same data series
Data(isinf(Data))=0;
DataRxns=zeros(size(Data));DataRxns(abs(Data)>0)=1;DataRxns=sum(DataRxns,2);
DataRxns=find(DataRxns==size(Data,2));
EqualData=zeros(size(Data,1));
for i=1:size(Data,1),
    for j=1:size(Data,1),
        if Data(i,:)==Data(j,:),
            EqualData(i,j)=1;
        end
    end
end
%Compute correlation matrix

Cor=corr(Data(DataRxns,:)');
Cor(EqualData(DataRxns,DataRxns)==1)=0;

%Get correlation between all reaction pairs and get distribution
Datafc=fctable(DataRxns,DataRxns);
DataCorTot=Cor;
DataCorTot=triu(DataCorTot,1);DataCorTot(DataCorTot==0)=[];

%Compute correlation for adjacent reaction pairs and get
%distribution
S=GEM.S;DataS=S(:,DataRxns);DataS(sum(DataS,2)==0,:)=[];RevRxns=find(GEM.rev==1);m1=1;m2=1;
for i=1:size(DataS,1),
    Rxns=find(DataS(i,:)~=0);
    N=length(Rxns);
    if N>1,
      for j=1:N-1,
          for k=j+1:N,
              AdjPairs1(m1,:)=Rxns([j,k]);
              DataCorAdj1(m1)=Cor(Rxns(j),Rxns(k));
              m1=m1+1;
          end
      end
    end
end
[~,Idx1]=unique(sort(AdjPairs1,2),'rows');
DataCorAdj1=DataCorAdj1(Idx1);
DataCorAdj1(isnan(DataCorAdj1) | DataCorAdj1==0)=[];

%Compute correlation and extract distribution for the fully (partially) and
%directionally coupled pairs only
DataCorFC=triu(Cor,1);
DataCorFC(Datafc~=1 | DataCorFC==0)=[];
DataCorDC=Cor;
DataCorDC(Datafc~=4 | DataCorDC==0)=[];
KStestFCTotcor=zeros(permutations,1);
KStestDCTotcor=zeros(permutations,1);
KStestAdj1Totcor=zeros(permutations,1);
KStestFCAdj1cor=zeros(permutations,1);
KStestDCAdj1cor=zeros(permutations,1);
for i=1:permutations,
    %FC vs Tot
    P1=randsample(1:length(DataCorTot),length(DataCorFC),'true');
    DataCorTotSamplesFC(i,:)=DataCorTot(P1);
    KStestFCTotcor(i)=kstest2(DataCorFC,DataCorTotSamplesFC(i,:),'Alpha',KSalphavalue,'Tail','smaller');
    %DC vs Tot
    P2=randsample(1:length(DataCorTot),length(DataCorDC),'true');
    DataCorTotSamplesDC(i,:)=DataCorTot(P2);
    KStestDCTotcor(i)=kstest2(DataCorDC,DataCorTotSamplesDC(i,:),'Alpha',KSalphavalue,'Tail','smaller');
    %Adj1 vs Tot
    P31=randsample(1:length(DataCorTot),length(DataCorAdj1),'true');
    DataCorTotSamplesAdj1(i,:)=DataCorTot(P31);
    KStestAdj1Totcor(i)=kstest2(DataCorAdj1,DataCorTotSamplesAdj1(i,:),'Alpha',KSalphavalue,'Tail','smaller');   
    %FC vs Adj1
    P4=randsample(1:length(DataCorAdj1),length(DataCorFC),'true');
    DataCorAdj1SamplesFC(i,:)=DataCorAdj1(P4);
    KStestFCAdj1cor(i)=kstest2(DataCorFC,DataCorAdj1SamplesFC(i,:),'Alpha',KSalphavalue,'Tail','smaller');
    %DC vs Adj1
    P5=randsample(1:length(DataCorAdj1),length(DataCorDC),'true');
    DataCorAdj1SamplesDC(i,:)=DataCorAdj1(P5);
    KStestDCAdj1cor(i)=kstest2(DataCorDC,DataCorAdj1SamplesDC(i,:),'Alpha',KSalphavalue,'Tail','smaller');
end

%Compute cumulative distribution of correlation and distance values
n=1;
Cmin=min(DataCorTot);
CorSeq=Cmin:0.01:1;

for i=CorSeq;
    PcTot(n)=length(find(DataCorTot>i))/length(DataCorTot);
    PcFC(n)=length(find(DataCorFC>i))/length(DataCorFC);
    PcDC(n)=length(find(DataCorDC>i))/length(DataCorDC);
    PcAdj1(n)=length(find(DataCorAdj1>i))/length(DataCorAdj1);
    n=n+1;
end

%Plot DataCorTot samples
for i=1:permutations,
    n=1;
    for j=CorSeq;
      PTotSamplescorFC(i,n)=length(find(DataCorTotSamplesFC(i,:)>j))/length(DataCorFC);
      PTotSamplescorAdj1(i,n)=length(find(DataCorTotSamplesAdj1(i,:)>j))/length(DataCorAdj1);
      PTotSamplescorDC(i,n)=length(find(DataCorTotSamplesDC(i,:)>j))/length(DataCorDC);
      Adj1SamplescorFC(i,n)=length(find(DataCorAdj1SamplesFC(i,:)>j))/length(DataCorFC);
      n=n+1;
   end
end

if printFigures=='T',
    %plot total pairs sample vs FC
    h(1)=figure('Color','w');
    g1=hggroup;
    g2=hggroup;
    g3=hggroup;
    for i=1:permutations,
        hold on
        plot(CorSeq,PTotSamplescorFC(i,:),'Color',[0.7,0.7,0.7],'LineWidth',0.3,'Parent',g1)
    end
    hold on
    plot(CorSeq,PcFC,'Color',[0.2081 0.1663 0.5292],'LineWidth',2,'Parent',g3)
    hold on
    plot(CorSeq,PcTot,'k','LineWidth',1,'Parent',g2)
    xlabel('Pearson correlation','FontSize',16)
    ylabel('density','FontSize',16)
    h_legend=legend([g1,g2,g3],'Total pairs samples','Total pairs','FC pairs','Location','SouthWest');
    set(h_legend,'FontSize',14);
    set(gca,'fontsize',15)
    legend boxoff   
    axis([Cmin 1 0 1])
%     print(h(1),[Directory,FigureName,'\','FCtot'],'-dtiff')
    
    %plot total pairs sample vs DC
    h(2)=figure('Color','w');
    g1=hggroup;
    g2=hggroup;
    g3=hggroup;
    for i=1:permutations,
        hold on
        plot(CorSeq,PTotSamplescorDC(i,:),'Color',[0.7,0.7,0.7],'LineWidth',0.3,'Parent',g1)
    end
    hold on
    plot(CorSeq,PcDC,'Color',[0.0590 0.6838 0.7254],'LineWidth',2,'Parent',g3)
    hold on
    plot(CorSeq,PcTot,'k','LineWidth',1,'Parent',g2)
    xlabel('Pearson correlation','FontSize',16)
    ylabel('density','FontSize',16)
    h_legend=legend([g1,g2,g3],'Total pairs samples','Total pairs','DC pairs','Location','SouthWest');
    set(h_legend,'FontSize',14);
    set(gca,'fontsize',15)
    legend boxoff 
    axis([Cmin 1 0 1])
%     print(h(2),[Directory,FigureName,'\','DCtot'],'-dtiff')

    %plot total pairs sample vs adjacent1
    h(3)=figure('Color','w');
    g1=hggroup;
    g2=hggroup;
    g3=hggroup;
    for i=1:permutations,
        hold on
        plot(CorSeq,PTotSamplescorAdj1(i,:),'Color',[0.7,0.7,0.7],'LineWidth',0.3,'Parent',g1)
    end
    hold on
    plot(CorSeq,PcAdj1,'g','LineWidth',2,'Parent',g3)
    hold on
    plot(CorSeq,PcTot,'k','LineWidth',1,'Parent',g2)
    xlabel('Pearson correlation','FontSize',16)
    ylabel('density','FontSize',16)
    h_legend=legend([g1,g2,g3],'Total pairs samples','Total pairs','Adj pairs','Location','SouthWest');
    set(h_legend,'FontSize',14);
    set(gca,'fontsize',15)
    legend boxoff 
    axis([Cmin 1 0 1])
%     print(h(3),[Directory,FigureName,'\','AdjTot'],'-dtiff')

    %plot adjacent1 pairs sample vs FC
    h(4)=figure('Color','w');
    g1=hggroup;
    g2=hggroup;
    g3=hggroup;
    for i=1:permutations,
        hold on
        plot(CorSeq,Adj1SamplescorFC(i,:),'Color',[0.7,1,0.7],'LineWidth',0.3,'Parent',g1)
    end
    hold on
    plot(CorSeq,PcFC,'Color',[0.2081 0.1663 0.5292],'LineWidth',2,'Parent',g3)
    hold on
    plot(CorSeq,PcAdj1,'g','LineWidth',1,'Parent',g2)
    xlabel('Pearson correlation','FontSize',16)
    ylabel('density','FontSize',16)
    h_legend=legend([g1,g2,g3],'Adj pairs samples','Adj pairs','FC pairs','Location','SouthWest');
    set(h_legend,'FontSize',14);
    set(gca,'fontsize',15)
    legend boxoff   
    axis([Cmin 1 0 1])
%     print(h(4),[Directory,FigureName,'\','FCAdj'],'-dtiff')

    %plot adjacent1 pairs sample vs DC
    h(5)=figure('Color','w');
    g1=hggroup;
    g2=hggroup;
    g3=hggroup;
    for i=1:permutations,
        hold on
        plot(CorSeq,Adj1SamplescorFC(i,:),'Color',[0.7,1,0.7],'LineWidth',0.3,'Parent',g1)
    end
    hold on
    plot(CorSeq,PcDC,'Color',[0.0590 0.6838 0.7254],'LineWidth',2,'Parent',g3)
    hold on
    plot(CorSeq,PcAdj1,'g','LineWidth',1,'Parent',g2)
    xlabel('Pearson correlation','FontSize',16)
    ylabel('density','FontSize',16)
    h_legend=legend([g1,g2,g3],'Adj pairs samples','Adj pairs','DC pairs','Location','SouthWest');
    set(h_legend,'FontSize',14);
    set(gca,'fontsize',15)
    legend boxoff   
    axis([Cmin 1 0 1])
%     print(h(5),[Directory,FigureName,'\','DCAdj'],'-dtiff')
    
    savefig(h,[FigureName,'.fig'],'compact')
    close all
end

%Get coupled groups and tables
if biology=='T',
  BIO=getBiology(GEM,Data,DataRxns,fctable,alpha,rhoCutoff,PcTot,DataCorTot,CorSeq,permutations);
  %Get subsystems per coupled group and rank them according to the proportion of reactions
  GroupSubSystems=getSubSystems(GEM,BIO);
  Val.Biology=BIO;
  Val.GroupSubSystems=GroupSubSystems;
%   Val.NFCgroupsrho=length(find(cell2mat(BIO.FCtablecellData(2:end,6))==1));
%   Val.NDCgroupsrho=length(find(cell2mat(BIO.DCtablecellData(2:end,6))==1));
  Val.rho=rhoCutoff;
end

% Val.DataCorTot=DataCorTot;
% Val.DataCorFC=DataCorFC;
% Val.DataCorDC=DataCorDC;
% Val.AdjPairsCornoDir=DataCorAdj1;

% Val.PTot=PcTot;
% Val.PFC=PcFC;
% Val.PDC=PcDC;
% Val.PAdj=PcAdj1;

Val.PvalueFC_Totcor=1-length(find(KStestFCTotcor==1))/permutations;
Val.PvalueDC_Totcor=1-length(find(KStestDCTotcor==1))/permutations;
Val.PvalueAdj_TotcorNoDir=1-length(find(KStestAdj1Totcor==1))/permutations;
Val.PvalueFC_Adj1cor=1-length(find(KStestFCAdj1cor==1))/permutations;
Val.PvalueDC_Adj1cor=1-length(find(KStestDCAdj1cor==1))/permutations;
     
end
