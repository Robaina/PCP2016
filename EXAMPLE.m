function Sol=EXAMPLE(WangGeneNames,WangLeafData)

    %This is an example on how to generate the results presented in this study,
    %using the maize leaf model and Leaf Data 1 (Wang et al. ,2014). 
    %The MaizeModel1, corresponding to the maize

    %get COBRA structure (Depends on COBRA toolbox)
    MaizeModel1=readCbModel('TS4_Model.xml');

    %eliminate blocked reactions and get reduced model
    MaizeModel1=myreduceModel(MaizeModel1);

    %assign reaction localization (Me, BS)
    MaizeModel1=getRxnLocation(MaizeModel1);

    %get coupling relations among reactions of the model (Depends on F2C2
    %toolbox and glpk solver)

    MaizeModel1Network.stoichiometricMatrix=full(MaizeModel1.S);
    MaizeModel1Network.reversibilityVector=MaizeModel1.rev;
    MaizeModel1Network.Reactions=MaizeModel1.rxns;
    MaizeModel1Network.Metabolites=MaizeModel1.mets;

    FCAMaizeModel1=F2C2('glpk',MaizeModel1Network);

    %Map original expression data onto metabolic model
    Maize1LeafData1=zeros(length(MaizeModel1.rxns),size(WangLeafData,2));
    for i=1:size(WangLeafData,2),
        Maize1LeafData1(:,i)=mapgene2rxn(MaizeModel1,WangGeneNames,WangLeafData(:,i));
    end

    %Analyse coupling relations and correspondance to gene co-expression
    CorMaize1LeafData1=MaizeAnalysis(Maize1Model1,Maize1LeafData1,FCAMaizeModel1,0.05,90,1,0.05,'T','SolMaizeLeafData1','T');
    
    %Analyse gene expression values associated to reactions in coupled
    %relations
    ExpMaize1LeafData1=getMedianExpression(MaizeModel1,FCAMaizeModel1,Maize1LeafData1,1,'ExpMaize1LeafData1');

    %Analyse coupling relations in the truncated maize leaf model
    MaizeModel1Truncated=getTruncatedmodel(MaizeModel1);

    MaizeModel1TruncatedNetwork.stoichiometricMatrix=full(MaizeModel1Truncated.S);
    MaizeModel1TruncatedNetwork.reversibilityVector=MaizeModel1Truncated.rev;
    MaizeModel1TruncatedNetwork.Reactions=MaizeModel1Truncated.rxns;
    MaizeModel1TruncatedNetwork.Metabolites=MaizeModel1Truncated.mets;

    FCAMaizeModel1Truncated=F2C2('glpk',MaizeModel1TruncatedNetwork);

    BioMaizeModel1Truncated=getBiologyNoData(MaizeModel1Truncated,FCAMaizeModel1Truncated);

    Sol.MaizeModel1=MaizeModel1;
    Sol.MaizeModel1Truncated=MaizeModel1Truncated;
    Sol.CorMaize1LeafData1=CorMaize1LeafData1;
    Sol.ExpMaize1LeafData1=ExpMaize1LeafData1;
    Sol.BioMaizeModel1Truncated=BioMaizeModel1Truncated;

end