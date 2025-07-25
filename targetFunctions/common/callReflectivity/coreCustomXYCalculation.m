function [reflectivity,simulation,shiftedData,sldProfile,layers,...
    resampledLayers,chi] = coreCustomXYCalculation(bulkIn,bulkOut,...
    shiftedData,simulationXData,dataIndices,background,resolution,...
    backgroundAction,parallel,resampleMinAngle,resampleNPoints,...
    roughness,contrastSld,nParams)

    % Resample the layers - always required for a custom XY calculation
    sldProfile = contrastSld(:,[1,2]);
    sldProfileIm = contrastSld(:,[1,3]);
    resampledLayers = resampleLayers(sldProfile,sldProfileIm,...
        resampleMinAngle,resampleNPoints);
    
    layers = resampledLayers;

    reflectivityType = 'standardAbeles';
    [reflectivity,simulation] = callReflectivity(bulkIn,bulkOut,...
     simulationXData,dataIndices,1,layers,roughness,resolution,parallel,...
     reflectivityType);

    [reflectivity,simulation,shiftedData] = applyBackgroundCorrection(...
        reflectivity,simulation,shiftedData,background,backgroundAction);

    chi = chiSquared(shiftedData,reflectivity,nParams);

end
