function problemStruct = checkRange(problemStruct, limits)
    % Check the range of fitted parameters, remove parameter from fit if range is too small.
    %
    % Parameters
    % ----------
    % problemStruct : struct
    %     The project struct.
    % limits : struct
    %     The limits for each parameter.
    %
    % Returns
    % -------
    % problemStruct : struct
    %     The project struct with fit information.
    fields = {"params", "backgroundParams", "scalefactors", "bulkIns",...
        "bulkOuts", "resolutionParams", "domainRatios"};
    titles = {"Parameter", "Background parameter", "Scalefactor", "Bulk in",...
        "Bulk out", "Resolution parameter", "Domain ratio"};

    for i = 1:length(fields)
        fitIndices = find(problemStruct.checks.(fields{i}));

        for j = 1:length(fitIndices)
            lower = limits.(fields{i})(fitIndices(j),1);
            upper = limits.(fields{i})(fitIndices(j),2);
            
            if (upper - lower) < 1e-10
                paramName = problemStruct.names.(fields{i}){fitIndices(j)};
                warning('%s "%s" was removed from the fit because its range is too small (< 1e-10).', titles{i}, paramName);
                problemStruct.checks.(fields{i})(fitIndices(j)) = 0;
            end
        end
    end

end

