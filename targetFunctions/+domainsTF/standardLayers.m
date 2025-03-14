function [qzshifts,scalefactors,bulkIns,bulkOuts,chis,reflectivity,...
    simulation,shiftedData,backgrounds,resolutions,domainLayerSlds,...
    domainSldProfiles,domainResampledLayers,subRoughs] = standardLayers(problemStruct,controls)
    % This is the main reflectivity calculation of the standard layers
    % calculation type. It extracts the required parameters for the contrasts
    % from the input arrays, then passes the main calculation to
    % 'coreLayersCalculation', which carries out the calculation itself. 
    % The core calculation is common for both standard and custom layers.
    
    % Extract parameters from problemStruct
    [numberOfContrasts, geometry, contrastBackgroundIndices, contrastQzshiftIndices,...
     contrastScalefactorIndices, contrastBulkInIndices, contrastBulkOutIndices,...
     contrastResolutionIndices, contrastDomainRatioIndices, backgroundParamArray,...
     qzshiftArray, scalefactorArray, bulkInArray, bulkOutArray, resolutionParamArray,...
     domainRatioArray, dataPresent, nParams, params, ~, resample,...
     contrastBackgroundTypes, contrastBackgroundActions, contrastResolutionTypes,...
     ~, useImaginary, repeatLayers, data, dataLimits, simLimits, contrastLayers,...
     layersDetails, customFiles, domainContrastLayers] = extractProblemParams(problemStruct);
    
    calcSld = controls.calcSldDuringFit;
    parallel = controls.parallel;
    resampleMinAngle = controls.resampleMinAngle;
    resampleNPoints = controls.resampleNPoints;
    
    % Allocate the memory for the output arrays before the main loop
    qzshifts = zeros(numberOfContrasts,1);
    scalefactors = zeros(numberOfContrasts,1);
    bulkIns = zeros(numberOfContrasts,1);
    bulkOuts = zeros(numberOfContrasts,1);
    subRoughs = zeros(numberOfContrasts,1);
    chis = zeros(numberOfContrasts,1);
    
    reflectivity = cell(numberOfContrasts,1);   
    simulation = cell(numberOfContrasts,1);
    shiftedData = cell(numberOfContrasts,1);
    backgrounds = cell(numberOfContrasts,1);
    resolutions = cell(numberOfContrasts,1);
    domainLayerSlds = cell(numberOfContrasts,2);
    domainSldProfiles = cell(numberOfContrasts,2);
    domainResampledLayers = cell(numberOfContrasts,2);

    layerSlds = cell(numberOfContrasts,1);
    sldProfiles = cell(numberOfContrasts,1);    
    resampledLayers = cell(numberOfContrasts,1);    
    
    domainContrastLayers1 = cell(numberOfContrasts,1);   
    domainContrastLayers2 = cell(numberOfContrasts,1);
    
    % First we need to allocate the absolute values of the input
    % parameters to all the layers in the layers list. This only needs
    % to be done once, and so is done outside the contrasts loop
    outParameterisedLayers = allocateParamsToLayers(params, layersDetails);
    
    % Substrate roughness is always first parameter for standard layers
    for i = 1:numberOfContrasts
        subRoughs(i) = params(1);
    end

    for i = 1:numberOfContrasts
        domainContrastLayers1{i} = domainContrastLayers{contrastLayers{i}(1)};
        domainContrastLayers2{i} = domainContrastLayers{contrastLayers{i}(2)};
    end

    if strcmpi(parallel, coderEnums.parallelOptions.Contrasts)
    
        % Loop over all the contrasts
        parfor i = 1:numberOfContrasts
        
            [qzshifts(i),scalefactors(i),bulkIns(i),bulkOuts(i),chis(i),...
             reflectivity{i},simulation{i},shiftedData{i},backgrounds{i},...
             resolutions{i},layerSlds{i},sldProfiles{i},resampledLayers{i}...
             ] = contrastCalculation(contrastBackgroundIndices{i},...
             contrastQzshiftIndices(i),contrastScalefactorIndices(i),...
             contrastBulkInIndices(i),contrastBulkOutIndices(i),...
             contrastResolutionIndices{i},contrastDomainRatioIndices(i),...
             backgroundParamArray,qzshiftArray,scalefactorArray,bulkInArray,...
             bulkOutArray,resolutionParamArray,domainRatioArray,dataPresent(i),...
             data{i},dataLimits{i},simLimits{i},repeatLayers{i},...
             contrastBackgroundTypes{i},contrastBackgroundActions{i},...
             contrastResolutionTypes{i},customFiles,nParams,parallel,...
             resampleMinAngle,resampleNPoints,resample(i),geometry,...
             subRoughs(i),calcSld,domainContrastLayers1{i},...
             domainContrastLayers2{i},outParameterisedLayers);
    
        end
    
    else
    
        % Loop over all the contrasts
        for i = 1:numberOfContrasts

            [qzshifts(i),scalefactors(i),bulkIns(i),bulkOuts(i),chis(i),...
             reflectivity{i},simulation{i},shiftedData{i},backgrounds{i},...
             resolutions{i},layerSlds{i},sldProfiles{i},resampledLayers{i}...
             ] = contrastCalculation(contrastBackgroundIndices{i},...
             contrastQzshiftIndices(i),contrastScalefactorIndices(i),...
             contrastBulkInIndices(i),contrastBulkOutIndices(i),...
             contrastResolutionIndices{i},contrastDomainRatioIndices(i),...
             backgroundParamArray,qzshiftArray,scalefactorArray,bulkInArray,...
             bulkOutArray,resolutionParamArray,domainRatioArray,dataPresent(i),...
             data{i},dataLimits{i},simLimits{i},repeatLayers{i},...
             contrastBackgroundTypes{i},contrastBackgroundActions{i},...
             contrastResolutionTypes{i},customFiles,nParams,parallel,...
             resampleMinAngle,resampleNPoints,resample(i),geometry,...
             subRoughs(i),calcSld,domainContrastLayers1{i},...
             domainContrastLayers2{i},outParameterisedLayers);
    
        end
    
    end
    
    for i = 1:numberOfContrasts
    
        contrastSldProfiles = sldProfiles{i};
        domainSldProfiles{i,1} = contrastSldProfiles{1};
        domainSldProfiles{i,2} = contrastSldProfiles{2};
    
        contrastLayerSlds = layerSlds{i};
        domainLayerSlds{i,1} = contrastLayerSlds{1};
        domainLayerSlds{i,2} = contrastLayerSlds{2};
    
        contrastResampledLayers = resampledLayers{i};
        domainResampledLayers{i,1} = contrastResampledLayers{1};
        domainResampledLayers{i,2} = contrastResampledLayers{2};
    
    end

    % Remove dummy imaginary column if present
    if ~useImaginary
        for i=1:numberOfContrasts
            domainLayerSlds{i,1}(:,3) = [];
            domainLayerSlds{i,2}(:,3) = [];
            domainResampledLayers{i,1}(:,3) = [];
            domainResampledLayers{i,2}(:,3) = [];
        end
    end

end


function [qzshiftValue,scalefactorValue,bulkInValue,bulkOutValue,chi,...
    reflectivity,simulation,shiftedData,background,resolution,layerSld,...
    sldProfile,resampledLayer] = contrastCalculation(backgroundParamIndex,...
    qzshiftIndex,scalefactorIndex,bulkInIndex,bulkOutIndex,resolutionParamIndex,...
    domainRatioIndex,backgroundParams,qzshifts,scalefactors,bulkIns,bulkOuts,...
    resolutionParams,domainRatios,dataPresent,data,dataLimits,simLimits,...
    repeatLayers,backgroundType,backgroundAction,resolutionType,...
    customFiles,nParams,parallel,resampleMinAngle,resampleNPoints,...
    resample,geometry,roughness,calcSld,domainContrastLayers1,...
    domainContrastLayers2,outParameterisedLayers)

    % Get the domain ratio for this contrast
    if isempty(domainRatioIndex)
        domainRatioIndex = 1;
    end
    domainRatio = domainRatios(domainRatioIndex);

    % Extract the relevant parameter values for this contrast
    % from the input arrays.
    % First need to decide which values of the backgrounds, scalefactors
    % data shifts and bulk contrasts are associated with this contrast
    [qzshiftValue,scalefactorValue,bulkInValue,bulkOutValue] = backSort(qzshiftIndex,...
     scalefactorIndex,bulkInIndex,bulkOutIndex,qzshifts,scalefactors,bulkIns,bulkOuts);

    % Apply scale factors and q shifts to the data
    shiftedData = shiftData(scalefactorValue,qzshiftValue,dataPresent,data,dataLimits,simLimits);
    [simulationXData, dataIndices] = makeSimulationRange(shiftedData, simLimits);

    background = constructBackground(backgroundType,backgroundParamIndex,...
        shiftedData,customFiles,backgroundParams,simulationXData,dataIndices);
    resolution = constructResolution(resolutionType,resolutionParamIndex,...
        shiftedData,customFiles,resolutionParams,simulationXData,dataIndices);

    % Also need to determine which layers from the overall layers list
    % are required for this contrast, and put them in the correct order 
    % according to geometry. We run it twice, once for each domain...
    thisContrastLayers1 = allocateLayersForContrast(domainContrastLayers1,outParameterisedLayers);
    thisContrastLayers2 = allocateLayersForContrast(domainContrastLayers2,outParameterisedLayers);
       
    % Call the core layers calculation - need to do this once for each
    % domain
    [sldProfile1,reflect1,simul1,~,layerSld1,resampledLayer1,~] = normalTF.coreLayersCalculation(thisContrastLayers1,roughness,...
     geometry,bulkInValue,bulkOutValue,resample,calcSld,shiftedData,simulationXData,dataIndices,repeatLayers,...
     resolution,background,backgroundAction,nParams,parallel,resampleMinAngle,resampleNPoints);

    [sldProfile2,reflect2,simul2,shiftedData,layerSld2,resampledLayer2,~] = normalTF.coreLayersCalculation(thisContrastLayers2,roughness,...
     geometry,bulkInValue,bulkOutValue,resample,calcSld,shiftedData,simulationXData,dataIndices,repeatLayers,...
     resolution,background,backgroundAction,nParams,parallel,resampleMinAngle,resampleNPoints);
    
    % Calculate the average reflectivities....
    [reflectivity,simulation] = domainsTF.averageReflectivity(reflect1,reflect2,simul1,simul2,domainRatio);

    % Get an overall chi-squared for the new averaged curve..
    chi = chiSquared(shiftedData,reflectivity,nParams);

    % Store returned values for this contrast in the output arrays.
    sldProfile = {sldProfile1, sldProfile2};
    layerSld = {layerSld1, layerSld2};
    resampledLayer = {resampledLayer1, resampledLayer2};

end
