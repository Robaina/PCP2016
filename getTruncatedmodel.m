function GEMtrunc=getTruncatedmodel(GEM)
    GEMtrunc=contextmodel2COBRA(setdiff(1:length(GEM.rxns),GEM.MeBSRxns),GEM,'GEMtrunc',1);
    GEMtrunc=myreduceModel(GEMtrunc);
    GEMtrunc=getRxnLocation(GEMtrunc);
end