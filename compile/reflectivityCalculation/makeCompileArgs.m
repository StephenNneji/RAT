function ARGS = makeCompileArgs()

% Define the arguments for compiling reflectivityCalculation
% using codegen.

%% Define argument types for entry-point 'reflectivityCalculation'.
maxArraySize = 10000;
maxDataSize = 10000;

ARGS = cell(1,1);
ARGS{1} = cell(2,1);
ARGS_1_1 = struct;
ARGS_1_1.TF = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_1.resample = coder.typeof(0,[1 maxArraySize],[0 1]);
ARG = coder.typeof(0,[maxDataSize 6],[1 0]);
ARGS_1_1.data = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARGS_1_1.dataPresent = coder.typeof(0,[1 maxArraySize],[0 1]);
ARG = coder.typeof(0,[1 2]);
ARGS_1_1.dataLimits = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARG = coder.typeof(0,[1 2]);
ARGS_1_1.simulationLimits = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARGS_1_1.numberOfContrasts = coder.typeof(0);
ARGS_1_1.geometry = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_1.useImaginary = coder.typeof(true);
ARGS_1_1.repeatLayers = coder.typeof(0, [1 maxArraySize],[0 1]);
ARG = coder.typeof(0,[1 5],[0 1]);
ARGS_1_1.contrastBackgroundParams = coder.typeof({ARG},[1 maxArraySize],[0 1]);
ARG = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_1.contrastBackgroundTypes = coder.typeof({ARG},[1 maxArraySize],[0 1]);
ARGS_1_1.contrastBackgroundActions = coder.typeof({ARG},[1 maxArraySize],[0 1]);
ARGS_1_1.contrastScalefactors = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.contrastBulkIns = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.contrastBulkOuts = coder.typeof(0,[1 maxArraySize],[0 1]);
ARG = coder.typeof(0,[1 5],[1 1]);
ARGS_1_1.contrastResolutionParams = coder.typeof({ARG},[1 maxArraySize],[0 1]);
ARG = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_1.contrastResolutionTypes = coder.typeof({ARG},[1 maxArraySize],[0 1]);
ARGS_1_1.backgroundParams = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.scalefactors = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.bulkIns = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.bulkOuts = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.resolutionParams = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.params = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.numberOfLayers = coder.typeof(0);
ARG = coder.typeof(0,[1 maxArraySize],[1 1]);
ARGS_1_1.contrastLayers = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARG = coder.typeof(0,[1 10],[1 1]);
ARGS_1_1.layersDetails = coder.typeof({ARG}, [maxArraySize 1],[1 1]);
ARG = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_1.customFiles = coder.typeof({ARG}, [1 maxArraySize], [0 1]);
ARGS_1_1.modelType = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_1.contrastCustomFiles = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.contrastDomainRatios = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.domainRatios = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.numberOfDomainContrasts = coder.typeof(0);
ARG = coder.typeof(0,[1 maxArraySize],[1 1]);
ARGS_1_1.domainContrastLayers = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARGS_1_1.fitParams = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.fitLimits = coder.typeof(0,[maxArraySize 2],[1 0]);
ARG = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_1.priorNames = coder.typeof({ARG}, [maxArraySize 1],[1 0]);
ARGS_1_1.priorValues = coder.typeof(0, [maxArraySize 3], [1 0]);
ARGS_1_1_names = struct;
ARG = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_1_names.params = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARGS_1_1_names.backgroundParams = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARGS_1_1_names.scalefactors = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARGS_1_1_names.bulkIns = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARGS_1_1_names.bulkOuts = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARGS_1_1_names.resolutionParams = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARGS_1_1_names.domainRatios = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARGS_1_1_names.contrasts = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARGS_1_1.names = coder.typeof(ARGS_1_1_names);
ARGS_1_1_checks = struct;
ARGS_1_1_checks.params = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1_checks.backgroundParams = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1_checks.scalefactors = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1_checks.bulkIns = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1_checks.bulkOuts = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1_checks.resolutionParams = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1_checks.domainRatios = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.checks = coder.typeof(ARGS_1_1_checks);
ARGS{1}{1} = coder.typeof(ARGS_1_1);
ARGS_1_2 = struct;
ARGS_1_2.procedure = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_2.parallel = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_2.calcSldDuringFit = coder.typeof(true);
ARGS_1_2.resampleMinAngle = coder.typeof(0);
ARGS_1_2.resampleNPoints = coder.typeof(0);
ARGS_1_2.display = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_2.xTolerance = coder.typeof(0);
ARGS_1_2.funcTolerance = coder.typeof(0);
ARGS_1_2.maxFuncEvals = coder.typeof(0);
ARGS_1_2.maxIterations = coder.typeof(0);
ARGS_1_2.updateFreq = coder.typeof(0);
ARGS_1_2.updatePlotFreq = coder.typeof(0);
ARGS_1_2.populationSize = coder.typeof(0);
ARGS_1_2.fWeight = coder.typeof(0);
ARGS_1_2.crossoverProbability = coder.typeof(0);
ARGS_1_2.strategy = coder.typeof(0);
ARGS_1_2.targetValue = coder.typeof(0);
ARGS_1_2.numGenerations = coder.typeof(0);
ARGS_1_2.nLive = coder.typeof(0);
ARGS_1_2.nMCMC = coder.typeof(0);
ARGS_1_2.propScale = coder.typeof(0);
ARGS_1_2.nsTolerance = coder.typeof(0);
ARGS_1_2.nSamples = coder.typeof(0);
ARGS_1_2.nChains = coder.typeof(0);
ARGS_1_2.jumpProbability = coder.typeof(0);
ARGS_1_2.pUnitGamma = coder.typeof(0);
ARGS_1_2.boundHandling = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_2.adaptPCR = coder.typeof(true);
ARGS_1_2.IPCFilePath = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS{1}{2} = coder.typeof(ARGS_1_2);

end
