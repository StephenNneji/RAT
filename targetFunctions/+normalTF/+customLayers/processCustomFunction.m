function [resampledLayers,subRoughs] = processCustomFunction(contrastBulkIns,contrastBulkOuts,...
    bulkInArray,bulkOutArray,cCustFiles,numberOfContrasts,customFiles,params,useImaginary)

    % Top-level function for processing custom layers for all the
    % contrasts.
    resampledLayers = cell(numberOfContrasts,1);
    subRoughs = zeros(numberOfContrasts,1);

    bulkOuts = bulkOutArray(contrastBulkOuts);

    for i = 1:numberOfContrasts     % TODO - the ambition is for parfor here, but would fail for Matlab and Python CM's..

        % Choose which custom file is associated with this contrast
        functionHandle = customFiles{cCustFiles(i)};

        % Find values of 'bulkIn' and 'bulkOut' for this
        % contrast
        thisBulkIn = bulkInArray(contrastBulkIns(i));
        thisBulkOut = bulkOuts(i);

        thisContrastLayers = [1 1 1 1]; % typeDef
        coder.varsize('thisContrastLayers',[10000 6],[1 1]);
        if isnan(str2double(functionHandle))
            [thisContrastLayers, subRoughs(i)] = callMatlabFunction(functionHandle, params, thisBulkIn, bulkOuts, i, 0);
        else
            [thisContrastLayers, subRoughs(i)] = callCppFunction(functionHandle, params, thisBulkIn, bulkOuts, i-1, -1);
        end

        % If SLD is real, add dummy imaginary column
        contrastLayersSize = size(thisContrastLayers);
        if ~useImaginary
            thisContrastLayers = [thisContrastLayers(:,1:2) zeros(contrastLayersSize(1), 1) thisContrastLayers(:,3:end)];
        end

        % If the output layers has 6 columns, then we need to do
        % the hydration correction (the user has not done it in the
        % custom function).
        thisContrastLayers = applyHydration(thisContrastLayers,thisBulkIn,thisBulkOut);
        resampledLayers{i} = thisContrastLayers;

    end

end
