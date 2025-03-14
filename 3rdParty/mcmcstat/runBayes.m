function output = runBayes(loop,nsimu,burnin,adaptint,params,problem)

problemStruct = problem{1};
controls = problem{2};

display = controls.display;

%data = problem.data;
% Arrange the data into the format mcmcstat requires 
% We have the same number of 'batches' as contrasts.
% Also need to pass problem in order to access this
% from the subfunctions.
numberOfContrasts = problemStruct.numberOfContrasts;

% Pre-allocate data to keep the compiler happy
% data = cell(1,numberOfContrasts);
% structDefineType = struct('ydata',[],'problem',problem);
% for i = 1:numberOfContrasts
%     data{i} =  structDefineType;
% end
data.data = cell(1,numberOfContrasts);
data.problem = problem;
for i = 1:numberOfContrasts
    thisData = problemStruct.data{i};
    if ~isempty(thisData)
        data.data{i} = [thisData(:,:)];
    end
    %data{i}.problem = problem;
end

% Make qcov based on the ranges of the parameters
qcov = [];
fitPars = problemStruct.fitParams;
fitConstr = problemStruct.fitLimits;
nParams = length(fitPars);

for i = 1:nParams
    thisConstr = fitConstr(i,:);
    qcov(i) = (abs(thisConstr(2) - thisConstr(1))*0.0001)^2;
end
qcov = diag(qcov);

% Define model and method options.
model.modelfun      = 'refModel';                  % will return a reflectivity curve
model.ssfun         = 'reflectivity_fitModel';     % will return chi squared
model.nbatch        = numberOfContrasts;
model.N = 1;

options.method      = 'dram';          % adaptation method (mh, am, dr, dram)
options.nsimu       = nsimu;           % no of simulation%
options.qcov        = qcov;            % proposal covariance (not sure why 11....)

options.adaptint    = adaptint;        % adaptation interval
options.printint    = 200;             % how often to show info of acceptance ratios
options.verbosity   = 1;               % how much to show output
options.waitbar     = 1;               % show graphical waitbar
options.updatesigma = 0;               % update error variance
options.stats       = 1;               % save extra statistics in result
options.burnintime  = burnin;          % burn in time..
options.ntry = 1;
options.drscale = [3 2 1];
%options.adascale = 2.4 / sqrt(nParams) * 0.01;


results = [];
loop = int32(loop);
for i = 1:loop
    if ~strcmpi(display, coderEnums.displayOptions.Off)
        fprintf('Running loop %d of %d ',i,loop);
    end

    %[results,chain,s2chain,sschain] = mcmcrun_compile(model, data, problem, params, options, results, display);
    [results,chain,s2chain,sschain] = mcmcrun(model, data, params, options, results);
    
    if ~strcmpi(display, coderEnums.displayOptions.Off)
        fprintf('\n');
    end
end

output.results = results;
output.chain = chain;
output.s2chain = s2chain;
output.sschain = sschain;
output.bestParams = results.mean;
output.data = data;

% out = mcmcpred_compile(results,chain,[],data,problem,500);
% outSld = mcmcpred_compile_sld(results,chain,[],data,problem,500);
% 
% problemStruct.fitParams = output.bestParams;
% problemStruct = unpackParams(problemStruct);
% result = reflectivityCalculation(problemStruct,controls);
% 
% output.bestFits = result.reflectivity;
% output.shiftedData = problemStruct.data;
% output.predictionIntervals = out;

end


%-------------------------------------------------------------------------

function makePredPlot(out,data,problem,bestFit,scalefac)


% figure(refFig);
% clf;
% hold on;

%data = problem.data;

thisData = data;

%figure(refFig); clf; hold on
errorbar(thisData(:,1),thisData(:,2),thisData(:,3),'.');
set(gca,'YScale','log','XScale','log');
hold on
plot(bestFit(:,1),bestFit(:,2).*scalefac,'k-');

numRanges = 1;
for i = 1:numRanges
    thisMin = out(2,:);
    thisMax = out(8,:);
    
    numPoints = size(thisData,1);
    
    plot(thisData(:,1),thisMin(1:numPoints),'-')
    plot(thisData(:,1),thisMax(1:numPoints),'-')
end

end




