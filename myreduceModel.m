function [reducedModel,Vrange]=myreduceModel(GEM,L,eps)

%Eliminate blocked reactions and extract reduced model
%Semidan, march 2014

if nargin<2,
   L=1e3;
end
if isempty(L),
   L=1e3;
end
if nargin<3,
   eps=1e-6;
end
if isempty(eps),
   eps=1e-6;
end

%perform FVA on GEM
[Vrange(:,1),Vrange(:,2)]=FVA(GEM,L,'gurobi');

%identify blocked reactions
Blocked=find(abs(Vrange(:,1))<=eps & abs(Vrange(:,2))<=eps);
if isempty(Blocked),
   Blocked=0;
end

%update reversibility vector and stoichiometric matrix
for i=1:length(GEM.rxns),
    if abs(Vrange(i,1))>eps && abs(Vrange(i,2))<eps,
        GEM.S(:,i)=-GEM.S(:,i);
    end
end
GEM.rev=zeros(length(GEM.rxns),1);
GEM.rev(abs(Vrange(:,1))>eps & abs(Vrange(:,2))>eps)=1;

%extract reduced model
reducedModel=contextmodel2COBRA(setdiff(1:length(GEM.rxns),Blocked),GEM,'reduced GEM',1);
reducedModel.Blocked=Blocked;

end