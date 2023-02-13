function [rxnList, stoichiometries] = yuv_findRxnsFromMets(model, mets)
% Get all Metabolites from a set of reactions.
%
% USAGE:
%
%    [metList, stoichiometries] = findMetsFromRxns(model, reactions)
%    [metList] = findMetsFromRxns(model, reactions)
%
% INPUTS:
%    model:             COBRA model structure
%    reactions:         Reaction IDs (Cell Array) or positions (Double array) to 
%                       find the corresponding metabolites for
%
% OUTPUT:
%    metList:           If only one output is requested, returns an ordered set of all
%                       metabolites involved in the provided reactions.
%                       Otherwise, this is a Cell Array of cell arrays
%                       containing the metabolites involved in each of the
%                       provided reactions.
%    stoichiometries:   this is a Cell array of double arrays of the
%                       stoichiometric coefficients corresponding to the
%                       reactions in the order of provided reaction ids.
%                       If reactions not in the model are provided, those
%                       will be represented by empty arrays.
% .. Author: - Thomas Pfau Jan 2018

if ~isnumeric(mets)
    metInd = findMetIDs(model, mets);
else
    metInd = mets;
end

metNotInModel = (metInd == 0);

if any(metNotInModel)
    warning('The following metabolities are not in the model:\n%s',strjoin(mets(metNotInModel),'; '));
end

metInd = metInd(~metNotInModel);
metStoich = model.S(metInd,:);

%if only reactions are requested
if nargout < 2
    rxnList = model.rxns(sum(abs(metStoich),1) > 0);
    return
else
    %Initialize the outputs.
    rxnList = cell(numel(metInd),1);
    stoichiometries = cell(numel(metInd),1);
    num_Rxnlist = cell(length(metInd),1); 
    %Init the relevant reactions.
    relRxnList = cell(numel(metInd),1);
    relStoichiometries = cell(numel(metInd),1);
    for i = 1:numel(metInd)
        relpos = metStoich(i,:) ~= 0;
        relRxnList{i} = model.rxns(relpos);
        relStoichiometries{i} = metStoich(i,relpos);

        num_Rxnlist{i} = find(relpos == 1); 

    end
    rxnList(~metNotInModel) = relRxnList;
    rxnList(metNotInModel) = {{}};
    rxnList{:,2} = cell2mat(num_Rxnlist); 
    stoichiometries(~metNotInModel) = relStoichiometries;
end

