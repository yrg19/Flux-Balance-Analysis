function [plotting_struct] = yuv_dynamicFBA(model, substrateRxns, initConcentrations, initBiomass, timeStep, nSteps, plotRxns, exclUptakeRxns,forced_to_zero, plot_yesno, flux_del)
% Modified by Yuval Gal Shohet, 2022 
% Performs dynamic FBA simulation using the static optimization approach
%
% USAGE:
%
%    [concentrationMatrix, excRxnNames, timeVec, biomassVec] = dynamicFBA(model, substrateRxns, initConcentrations, initBiomass, timeStep, nSteps, plotRxns, exclUptakeRxns)
%
% INPUTS:
%    model:                  COBRA model structure
%    substrateRxns:          List of exchange reaction names for substrates
%                            initially in the media that may change (e.g. not
%                            h2o or co2)
%    initConcentrations:     Initial concentrations of substrates (in the same
%                            structure as `substrateRxns`)
%    initBiomass:            Initial biomass (must be non zero)
%    timeStep:               Time step size
%    nSteps:                 Maximum number of time steps
%    plt:                    Plotting 
%    forced_to_zero:         force concentration to be zero, despite initial flux 
%    flux_del:               force flux to be zero after every iteration.
%
% OPTIONAL INPUTS:
%    plotRxns:               Reactions to be plotted (Default = {'EX_glc(e)', 'EX_ac(e)', 'EX_for(e)'})
%    exclUptakeRxns:         List of uptake reactions whose substrate concentrations do not change
%                            (Default = {'EX_co2(e)', 'EX_o2(e)', 'EX_h2o(e)', 'EX_h(e)'})
%
% OUTPUTS:
%    concentrationMatrix:    Matrix of extracellular metabolite concentrations
%    excRxnNames:            Names of exchange reactions for the EC metabolites
%    timeVec:                Vector of time points
%    biomassVec:             Vector of biomass values
%
% If no initial concentration is given for a substrate that has an open
% uptake in the model (i.e. `model.lb < 0`) the concentration is assumed to
% be high enough to not be limiting. If the uptake rate for a nutrient is
% calculated to exceed the maximum uptake rate for that nutrient specified
% in the model and the max uptake rate specified is > 0, the maximum uptake
% rate specified in the model is used instead of the calculated uptake
% rate.
%
% NOTE:
%
%    The dynamic FBA method implemented in this function is essentially
%    the same as the method described in
%    [`Varma, A., and B. O. Palsson. Appl. Environ. Microbiol. 60:3724 (1994)`].
%    This function does not implement the dynamic FBA using dynamic optimization approach
%    described in [`Mahadevan, R. et al. Biophys J, 83:1331-1340 (2003)`].
%
% .. Author: - Markus Herrgard 8/22/06

global WAITBAR_TYPE

if (nargin < 7) || isempty(plotRxns)
    plotRxns = {'EX_glc(e)','EX_ac(e)','EX_for(e)'};
end

% Uptake reactions whose substrate concentrations do not change
if (nargin < 8) || isempty(exclUptakeRxns)
    exclUptakeRxns = {'EX_co2(e)','EX_h2o(e)','EX_h(e)'};
end

%%
% Find exchange rxns
excInd = findExcRxns(model,false);
excInd = excInd & ~ismember(model.rxns,exclUptakeRxns);
excRxnNames = model.rxns(excInd);

% Figure out if substrate reactions are correct
missingInd = find(~ismember(substrateRxns,excRxnNames));
if (~isempty(missingInd))
    for i = 1:length(missingInd)
        fprintf('%s\n',substrateRxns{missingInd(i)});
    end
    error('Invalid substrate uptake reaction!');
end

% Initialize concentrations
[~, substrateMatchInd] = ismember(substrateRxns,excRxnNames);
concentrations = zeros(length(excRxnNames),1);
concentrations(substrateMatchInd) = initConcentrations;

% Deal with reactions for which there are no initial concentrations
originalBound = -model.lb(excInd);
noInitConcentration = (concentrations == 0 & originalBound > 0);
concentrations(noInitConcentration) = 1000;

biomass = initBiomass;

% Initialize bounds
uptakeBound =  concentrations/(biomass*timeStep);

% Make sure bounds are not higher than what are specified in the model
aboveOriginal = (uptakeBound > originalBound) & (originalBound > 0);
uptakeBound(aboveOriginal) = originalBound(aboveOriginal);

model.lb(excInd) = -uptakeBound;

% flux reaction deletion 
if ~isempty(flux_del)
    if isnumeric(flux_del)
        forcedInd = flux_del;
    else
        [~, forcedInd] = ismember(flux_del,model.rxns);
    end
    model.lb(forcedInd) = 0;
end

% if forced to zero make zero 
if ~isempty(forced_to_zero)
    [~, forcedInd] = ismember(forced_to_zero,excRxnNames); 
    concentrations(forcedInd) = 0; 
end 



%%
concentrationMatrix = sparse(concentrations);
biomassVec = biomass;
timeVec(1) = 0;

real_indxs = [find(contains(model.rxns',substrateRxns(1))),find(contains(model.rxns',substrateRxns(2)))];

fprintf('Step number\tBiomass\n');
showprogress(0,'Dynamic FBA analysis in progress ...');
uptakevec = model.lb(real_indxs);
muvec = []; 
concvec = concentrations(substrateMatchInd); 
fluxes_per_time(:,1) = model.lb; 
fluxes_per_met(:,1) = zeros(length(model.mets), 1); 
for stepNo = 1:nSteps
    % Run FBA
    sol = optimizeCbModel(model,'max','one');
    mu = sol.f; % growth rate after solving for model 
    if (sol.stat ~= 1 || mu == 0)
        fprintf('\nNo feasible solution - nutrients exhausted. Biomass:\t %f\n', biomass);
       break;
    end
    muvec = [muvec, mu]; 
    uptakeFlux = sol.x(excInd);
    uptakevec = [uptakevec, sol.x(real_indxs)];
    biomass = biomass*exp(mu*timeStep);
    %biomass = biomass*(1+mu*timeStep);
    biomassVec(end+1) = biomass;

    % Update concentrations
    concentrations = concentrations - uptakeFlux/mu*biomass*(1-exp(mu*timeStep));
    concvec = [concvec, concentrations(substrateMatchInd)];
    concentrations(concentrations <= 0) = 0;
    concentrationMatrix(:,end+1) = sparse(concentrations);

    % Update bounds for uptake reactions
    uptakeBound =  concentrations/(biomass*timeStep);
    % This is to avoid any numerical issues
    uptakeBound(uptakeBound > 1000) = 1000;
    % Figure out if the computed bounds were above the original bounds
    aboveOriginal = (uptakeBound > originalBound) & (originalBound > 0);
    % Revert to original bounds if the rate was too high
    uptakeBound(aboveOriginal) = originalBound(aboveOriginal);
    uptakeBound(abs(uptakeBound) < 1e-9) = 0;
    
    model.lb(excInd) = -uptakeBound;

        % flux reaction deletion 
    if ~isempty(flux_del)
        if isnumeric(flux_del)
            forcedInd = flux_del; 
        else 
            [~, forcedInd] = ismember(flux_del,model.rxns); 
        end 
        model.lb(forcedInd) = 0; 
    end 

    if WAITBAR_TYPE ~= 1
        fprintf('%d\t%f\n',stepNo,biomass);
    end
    showprogress(stepNo/nSteps);
    timeVec(stepNo+1) = stepNo*timeStep;

    fluxes_per_time(:,stepNo+1) = sol.v; 
    fluxes_per_met(:,stepNo+1) = sol.y; 
end

plotting_struct.concentrationMatrix = concentrationMatrix;
plotting_struct.excRxnNames = excRxnNames; 
plotting_struct.timeVec = timeVec; 
plotting_struct.biomassVec = biomassVec; 
plotting_struct.nSteps = nSteps; 
plotting_struct.substrateMatchInd = substrateMatchInd; 
plotting_struct.initConcentrations = initConcentrations; 
plotting_struct.uptakevec = uptakevec; 
plotting_struct.concvec = concvec; 
plotting_struct.fluxes = fluxes_per_time; 
plotting_struct.metabolites = fluxes_per_met; 

if plot_yesno == 1
    plot_dfba(plotting_struct)
end 



% xylosevec = concentrationMatrix(substrateMatchInd(1),:); 
% xylosevec = (xylosevec/initConcentrations(1)) .*100; 
% glucosevec = concentrationMatrix(substrateMatchInd(2),:); 
% glucosevec = (glucosevec/initConcentrations(2)) .*100; 
% 
% subplot(2,1,1)
% plot(timeVec, xylosevec)
% ylim([0,100])
% title('Xylose v time')
% xlabel('Time, hr')
% ylabel('Percentage of Sugar Utilised')
% subplot(2,1,2)
% plot(timeVec, glucosevec)
% ylim([0,100])
% title('Glucose v time')
% xlabel('Time, hr')
% ylabel('Percentage of Sugar Utilised')



% %%
% selNonZero = any(concentrationMatrix>0,2);
% concentrationMatrix = concentrationMatrix(selNonZero,:);
% excRxnNames = excRxnNames(selNonZero);
% selPlot = ismember(excRxnNames,plotRxns);

% %% Plot concentrations as a function of time
% clf
% subplot(1,2,1);
% plot(timeVec,biomassVec);
% axis tight
% title('Biomass');
% subplot(1,2,2);
% plot(timeVec,concentrationMatrix(selPlot,:));
% axis tight
% legend(strrep(excRxnNames(selPlot),'EX_',''));
