function [sldProfile,reflect,simulation,shiftedData,theseLayers,resamLayers,chiSq] = ...
    coreLayersCalculation(layers,rough,geometry,bulkIn,bulkOut,resample,...
    calcSld,shiftedData,simulationXData,dataIndices,repeatLayers,resolution,...
    background,backgroundAction,params,parallelPoints,resampleMinAngle,...
    resampleNPoints)

%   This is the main reflectivity calculation for all layers models in the 
%   normal target function.
%
%   The function first arranges the roughness' in the correct order
%   according to geometry. Then, if required it calculates the SLD profile
%   and if requested resamples this into a number of zero-roughness layers
%   (roughness resampling). It the applies any scalefactors and qz shifts
%   the data requires. This is done before calculating the reflectivity to
%   ensure that the reflectivity is calculated over the same range in qz as
%   the shifted datapoints. The main reflectivity calculation is then
%   called, including the resolution function. The calculation outputs two
%   profiles - 'reflect' which is the same range as the points, and
%   'simulation' which can be a different range to allow extrapolation.
%   The background correction is the applied, and finally chi-squared is 
%   calculated.

% Pre-definition for Coder
sldProfile = [0 0];
sldProfileIm = [0 0];

% Build up the layers matrix for this contrast
[theseLayers, ssubs] = groupLayersMod(layers,rough,geometry,bulkIn,bulkOut);

% Make the SLD profiles.
% If resampling is needed, then enforce the calcSLD flag, so as to catch
% the error af trying to resample a non-existent profile.
if (resample == 1 && ~calcSld)
    calcSld = true;
end

% If calc SLD flag is set, then calculate the SLD profile
if calcSld

    % We process real and imaginary parts of the SLD separately
    thisSldLays = [theseLayers(:,1:2) theseLayers(:,4:end)];
    thisSldLaysIm = [theseLayers(:,1) theseLayers(:,3:end)];

    % Note bulkIn and bulkOut = 0 since there is never any imaginary part
    % for the bulk phases.
    sldProfile = makeSLDProfiles(bulkIn,bulkOut,thisSldLays,ssubs,repeatLayers);
    sldProfileIm = makeSLDProfiles(0,0,thisSldLaysIm,ssubs,repeatLayers);

end

% If required, then resample the SLD
if resample == 1
    layerSld = resampleLayers(sldProfile,sldProfileIm,resampleMinAngle,resampleNPoints);
    resamLayers = layerSld;
else
    layerSld = theseLayers;
    resamLayers = [0 0 0 0];
end

% Calculate the reflectivity
reflectivityType = 'standardAbeles';
[reflect,simulation] = callReflectivity(bulkIn,bulkOut,simulationXData,dataIndices,repeatLayers,layerSld,ssubs,resolution,parallelPoints,reflectivityType);

% Apply background correction
[reflect,simulation,shiftedData] = applyBackgroundCorrection(reflect,simulation,shiftedData,background,backgroundAction);

% Calculate chi squared.
chiSq = chiSquared(shiftedData,reflect,params);

end
