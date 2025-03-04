classdef parametersClass < tableUtilities
    % ``parametersClass`` manages the parameters for the project. It provides methods to add, update and remove parameters.
    % Each parameter is stored as a row in a table and consist of a name, a value with its minimum and maximum limits, a flag indicating 
    % if the parameter should be fitted in the calculation, and information for Bayesian calculations such as whether the prior is uniform 
    % or gaussian, then mu and sigma for gaussian priors. ``parametersClass`` will be initialised with a default first parameter if no 
    % arguments are provided otherwise the provided arguments will be used to create the first parameter.  
    %
    % Examples
    % --------
    % Default values are used when adding the parameter if no arguments are provided.
    % >>> params = parametersClass();
    % for 2 inputs, the min value provided would be used to set the value and max value.
    % >>> params = parametersClass('Tails', 10);
    % The code above is equivalent to 
    % >>> params = parametersClass('Tails', 10, 10, 10);
    % Other ways of initialisation
    % >>> params = parametersClass('Tails', 10, 20, 30, true);
    % >>> params = parametersClass('Tails', 10, 20, 30, true, priorTypes.Uniform.value, 0, Inf);
    %
    % Parameters
    % ----------
    % name : string or char array, default: auto-generated name
    %     The name of the first parameter. 
    % min : double, default: 0.0
    %     The minimum value that the first parameter could take when fitted.
    % value : double, default: 0.0
    %     The value of the parameter
    % max : double, default: 0.0
    %     The maximum value that the first parameter could take when fitted.
    % fit : logical, default: false
    %     Whether the first parameter should be fitted in a calculation.
    % priorType : PriorTypes, default: PriorTypes.Uniform 
    %     For Bayesian calculations, whether the prior likelihood is assumed to be ‘uniform’ or ‘gaussian’.
    % mu : double, default: 0
    %     If the prior type is Gaussian, the mu and sigma values describing the Gaussian function for the prior likelihood.
    % sigma : double, default: Inf
    %     If the prior type is Gaussian, the mu and sigma values describing the Gaussian function for the prior likelihood.
    %
    % Attributes
    % ----------
    % varTable : table
    %     The table which contains the properties for each parameter. 

    methods
        function obj = parametersClass(varargin)
            sz = [0, 8];
            varTypes = {'string','double','double','double','logical','string','double','double'};
            varNames = {'Name','Min','Value','Max','Fit?','Prior Type','mu','sigma'};
            obj.varTable = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
            if isempty(varargin)
                obj.addParameter();
            else
                obj.addParameter(varargin{:});
            end
        end
        
        function obj = addParameter(obj, name, min, value, max, fit, priorType, mu, sigma)
            % Adds an new parameter to the parameters table. 
            %
            % Examples
            % --------
            % To add a new parameter with no properties and an autogenerated name.
            % >>> params.addParameter();
            % To add a parameter with all available properties.
            % >>> params.addParameter('Tails', 20, 50, 60, true, 'gaussian', 1, 5);
            % Other examples of adding parameters with a subset of properties.
            % >>> params.addParameter('Tails');  % Parameter name only with others set to default
            % >>> params.addParameter('Tails', 23);  % Parameter name and min only. Value and max will be set to 23 to keep limits valid
            % >>> params.addParameter('Tails', 23, 24, 25, true);  % priors will be default
            %
            % Parameters
            % ----------
            % name : string or char array, default: auto-generated name
            %     The name of the parameter. 
            % min : double, default: 0.0
            %     The minimum value that the parameter could take when fitted.
            % value : double, default: 0.0
            %     The value of the parameter, default will be equal to ``min`` if this is not set.
            % max : double, default: 0.0
            %     The maximum value that the parameter could take when fitted, default will be equal to ``value`` if this is not set.
            % fit : logical, default: false
            %     Whether the parameter should be fitted in a calculation.
            % priorType : PriorTypes, default: PriorTypes.Uniform 
            %     For Bayesian calculations, whether the prior likelihood is assumed to be ‘uniform’ or ‘gaussian’.
            % mu : double, default: 0.0
            %     If the prior type is Gaussian, the mu and sigma values describing the Gaussian function for the prior likelihood.
            % sigma : double, default: Inf
            %     If the prior type is Gaussian, the mu and sigma values describing the Gaussian function for the prior likelihood.
            arguments
                obj
                name {mustBeTextScalar} = ''
                min {mustBeNumeric, isscalar} = 0.0
                value {mustBeScalarOrEmpty, mustBeNumeric} = []
                max {mustBeScalarOrEmpty, mustBeNumeric} = []
                fit {mustBeA(fit, 'logical')} = false
                priorType = priorTypes.Uniform
                mu {mustBeNumeric, isscalar} = 0.0
                sigma {mustBeNumeric, isscalar} = Inf
            end
            
            if isempty(value) && isempty(max)
                max = min;
                value = min;
            elseif ~isempty(value) && isempty(max)
                max = value;
            end
            values = [min, value, max];
            priorValues = [mu, sigma];

            if isempty(name)
                % No input parameter - create name and add defaults
                name = sprintf('new parameter %d',obj.autoNameCounter);
                newRow = {name, values(1), values(2), values(3), fit, priorType.value, priorValues(1), priorValues(2)};
            else                                  
                obj.validateLimits(values(1), values(2), values(3));
                priors = obj.validatePriors(priorType, priorValues(1), priorValues(2));
                newRow = {name, values(1), values(2), values(3), fit, priors{1}, priors{2}, priors{3}};    
            end
            obj.addRow(newRow{:});
        end
        
        function obj = removeParameter(obj, row)
            % Removes a parameter from the parameters table. 
            %
            % Examples
            % --------
            % To remove the second parameter in the table (parameter in row 2).  
            % >>> params.removeParameter(2);
            % To remove parameter with a specific name.
            % >>> params.removeParameter('Tails');
            %
            % Parameters
            % ----------
            % row : string or char array or integer
            %     If ``row`` is an integer, it is the row number of the parameter to remove. If it is text, 
            %     it is the name of the parameter to remove.
            obj.removeRow(row);
        end
        
        function obj = setParameter(obj, row, options)
            % General purpose method for updating properties of an existing parameter. 
            % Any unset property will remain unchanged.
            %
            % Examples
            % --------
            % To change the name and value of the second parameter in the table (parameter in row 2).
            % >>> params.setParameter(2, 'name', 'Heads', 'value', 50);
            % To change the all properties of a parameter called 'Tails'.
            % >>> params.setParameter('Tails', 'name', 'Heads', 'min', 20, 'value', 50, 'max', 60, ...
            % >>>                     'fit', true, 'priorType', 'gaussian', 'mu', 1, 'sigma', 5);
            % 
            % Parameters
            % ----------
            % row : string or char array or integer
            %     If ``row`` is an integer, it is the row number of the parameter to update. If it is text, 
            %     it is the name of the parameter to update.
            % options
            %    Keyword/value pair to properties to update for the specific parameter.
            %       * name (char array or string, default: '') the new name of the parameter.
            %       * min (double, default: []) the new minimum value of the parameter.
            %       * value (double, default: []) the new value of the parameter.
            %       * max (double, default: []) the new maximum value of the parameter.
            %       * fit (logical, default: logical.empty()) the new fit flag of the parameter.
            %       * priorTypes (priorTypes, priorTypes.empty()) the new prior type of the parameter.            
            %       * mu (double, default: []) the new mu value describing the Gaussian function for the prior.            
            %       * sigma (double, default: []) The new sigma value describing the Gaussian function for the prior.            
            arguments
                obj
                row
                options.name {mustBeTextScalar} = ''
                options.min {mustBeScalarOrEmpty, mustBeNumeric} = []
                options.value {mustBeScalarOrEmpty, mustBeNumeric} = []
                options.max {mustBeScalarOrEmpty, mustBeNumeric} = []
                options.fit {mustBeScalarOrEmptyLogical} = logical.empty()
                options.priorType = priorTypes.empty()
                options.mu {mustBeScalarOrEmpty, mustBeNumeric} = []
                options.sigma {mustBeScalarOrEmpty, mustBeNumeric} = []
            end
            row = obj.getValidRow(row);
            if isempty(options.name)
                options.name = obj.varTable{row, 1}{:};
            end
            
            if isempty(options.min)
                options.min = obj.varTable{row, 2};
            end
            
            if isempty(options.value)
                options.value = obj.varTable{row, 3};
            end
            
            if isempty(options.max)
                options.max = obj.varTable{row, 4};
            end
            
            if isempty(options.fit)
                options.fit = obj.varTable{row, 5};
            end
            
            if isempty(options.priorType)
                options.priorType = obj.varTable{row, 6};
            end
            
            if isempty(options.mu)
                options.mu = obj.varTable{row, 7};
            end
            
            if isempty(options.sigma)
                options.sigma = obj.varTable{row, 8};
            end

            obj.validateLimits(options.min, options.value, options.max);
            % Apply values
            obj.setName(row, options.name);
            obj.varTable{row, 2} = options.min;
            obj.varTable{row, 3} = options.value;
            obj.varTable{row, 4} = options.max;
            obj.setFit(row, options.fit);
            obj.setPrior(row, options.priorType, options.mu, options.sigma);
        end
        
        function obj = setPrior(obj, row, priorType, mu, sigma)
            % Sets the prior information of an existing parameter.
            %
            % Examples
            % --------
            % To change the prior of the second parameter in the table (parameter in row 2)
            % >>> params.setPrior(2, priorTypes.Gaussian, 1, 2);
            % To change the prior of a parameter called 'Tails'
            % >>> params.setPrior('Tails', 'uniform');
            %
            % Parameters
            % ----------
            % row : string or char array or integer
            %     If ``row`` is an integer, it is the row number of the parameter to update. If it is text, 
            %     it is the name of the parameter to update.
            % priorType : PriorTypes 
            %     The new prior type of the parameter.
            % mu : double, default: []
            %     The new mu value describing the Gaussian function for the prior likelihood, the value is not changed if empty.
            % sigma : double, default: []
            %     The new sigma value describing the Gaussian function for the prior likelihood, the value is not changed if empty.
            arguments
                obj
                row
                priorType
                mu {mustBeScalarOrEmpty, mustBeNumeric} = []
                sigma {mustBeScalarOrEmpty, mustBeNumeric} = []
            end
            
            row = obj.getValidRow(row);
            if isempty(mu)
                mu = obj.varTable{row, 7};
            end
            if isempty(sigma)
                sigma = obj.varTable{row, 8};
            end
            
            priors = obj.validatePriors(priorType, mu, sigma);
            obj.varTable{row, 6} = priors(1);
            obj.varTable{row, 7} = priors{2};
            obj.varTable{row, 8} = priors{3};
        end
        
        function obj = setValue(obj, row, value)
            % Sets the value of an existing parameter.
            %
            % Examples
            % --------
            % To change the value of the second parameter in the table (parameter in row 2)
            % >>> params.setValue(2, 3.4);
            % To change the value of a parameter called 'Tails'
            % >>> params.setValue('Tails', 3.4);
            %
            % Parameters
            % ----------
            % row : string or char array or integer
            %     If ``row`` is an integer, it is the row number of the parameter to update. If it is text, 
            %     it is the name of the parameter to update.
            % value : double
            %     The new value of the parameter.
            arguments
                obj
                row
                value {isscalar, mustBeNumeric}
            end
            row = obj.getValidRow(row);
            obj.validateLimits(obj.varTable{row, 2}, value, obj.varTable{row, 4});
            obj.varTable{row, 3} = value;
        end
        
        function obj = setName(obj, row, name)
            % Sets the name of an existing parameter.
            %
            % Examples
            % --------
            % To change the name of the second parameter in the table (parameter in row 2)
            % >>> params.setName(2, 'new name');
            % To change the name of a parameter called 'Tails' to 'Heads'
            % >>> params.setName('Tails', 'Heads');
            %
            % Parameters
            % ----------
            % row : string or char array or integer
            %     If ``row`` is an integer, it is the row number of the parameter to update. If it is text, 
            %     it is the name of the parameter to update.
            % name : string or char array
            %     The new name of the parameter.
            arguments
                obj
                row
                name {mustBeTextScalar}
            end
            row = obj.getValidRow(row);
            obj.varTable{row, 1} = {name};
        end
        
        function obj = setLimits(obj, row, min, max)
            % Sets the limits of an existing parameter. 
            %
            % Examples
            % --------
            % To change the limits of the second parameter in the table (parameter in row 2)
            % >>> params.setLimits(2, 0, 100);
            % To change the limits of a parameter with name 'Tails'
            % >>> params.setLimits('Tails', 0, 100);
            %
            % Parameters
            % ----------
            % row : string or char array or integer
            %     If ``row`` is an integer, it is the row number of the parameter to update. If it is text, 
            %     it is the name of the parameter to update.
            % min : double
            %     The new minimum value of the parameter.
            % max : double
            %     The new maximum value of the parameter.
            arguments
                obj
                row
                min {isscalar, mustBeNumeric}
                max {isscalar, mustBeNumeric}
            end
            row = obj.getValidRow(row);
            obj.validateLimits(min, obj.varTable{row, 3}, max);
       
            obj.varTable{row, 2} = min;
            obj.varTable{row, 4} = max;
        end
                
        function obj = setFit(obj, row, fit)
            % Sets the fit to off or on for an existing parameter.
            % 
            % Examples
            % --------
            % To change the fit flag of the second parameter in the table (parameter in row 2)
            % >>> params.setFit(2, 0, 100);
            % To change the fit flag of a parameter with name 'Tails'
            % >>> params.setFit('Tails', 0, 100);
            %
            % Parameters
            % ----------
            % row : string or char array or integer
            %     If ``row`` is an integer, it is the row number of the parameter to update. If it is text, 
            %     it is the name of the parameter to update.
            % fit : logical
            %     The new fit flag of the parameter.
            arguments
                obj
                row
                fit {mustBeA(fit, 'logical')}
            end
            row = obj.getValidRow(row);
           
            obj.varTable{row, 5} = fit;
        end
        
        function displayTable(obj, showPriors)
            % Prints the parameter table to the console.
            %
            % Examples
            % --------
            % To print the table with the prior information.
            % >>> params.displayTable(true);
            %
            % Parameters
            % ----------
            % showPriors : logical, default: false
            %     Indicates if the prior type, mu, and sigma columns should be displayed.
            arguments
                obj
                showPriors {mustBeA(showPriors, 'logical')} = false
            end

            numParams = height(obj.varTable);
            dim = [1, width(obj.varTable)];
                        
            if ~showPriors
                dim(2) = 5;
            end
            
            if numParams == 0    
                varNames = obj.varTable.Properties.VariableNames(1:dim(2));
                array = table('Size', dim, 'VariableTypes', repmat({'string'}, dim), 'VariableNames', varNames);
                array(1, :) = repmat({''}, 1, dim(2));
            else
                array = obj.varTable;
                p = 1:height(array);
                p = p(:);
                p = table(p);
                array = [p array(:, 1:dim(2))];
            end         
            disp(array);
        end
        
        function outStruct = toStruct(obj)
            % Converts the class parameters into a structure array.
            %
            % Returns
            % -------
            % outStruct : struct
            %     A struct which contains the properties for all the parameters.
            names = table2cell(obj.varTable(:,1));
            
            % Want these to be class 'char' rather than 'string'
            for n = 1:length(names)
                names{n} = char(names{n});
            end
            outStruct.names = names;
            
            mins = obj.varTable{:,2};
            maxs = obj.varTable{:,4};
            limits = cell([1, length(mins)]);
            for i = 1:length(mins)
                limits{i} = [mins(i) maxs(i)];
            end
            
            outStruct.limits = limits;
            
            outStruct.values = obj.varTable{:,3};
            
            outStruct.fit = double(obj.varTable{:,5});
            
            priors = table2cell(obj.varTable(:,6:8));
            priors = [outStruct.names priors];
            
            % Group each row into one cell. Should be a way of doing this
            % without a loop but I can't quite see it right now...
            pp = cell([1, size(priors, 1)]);
            for i = 1:size(priors, 1)
                thisPrior = priors(i,:);
                val2 = char(thisPrior{2});
                thisPrior{2} = val2;
                pp(i) = {thisPrior};
            end
            
            outStruct.priors = pp(:);
            
            
            % Need to force some of the outputs
            % to be row vectors, so transpose them
            outStruct.names = outStruct.names';
            outStruct.values = outStruct.values';
            outStruct.fit = outStruct.fit';
            
            % Fields order needs to be...
            
            % names
            % limits
            % values
            % fit
            % priors    
        end
        
    end
    
    methods (Access = protected)
        
        function index = getValidRow(obj, row)
            % Gets a valid row index, if ``row`` is a parameter name or validates the index, if ``row`` is an integer.
            % 
            % Examples
            % --------
            % To validate a row index is in table bounds.
            % >>> index = obj.getValidRow(2);
            % To find the index for a name.
            % >>> index = obj.getValidRow('Tails');
            %
            % Parameters
            % ----------
            % row : string or char array or integer
            %     If ``row`` is an integer, it is the row number of the parameter to update. If it is text, 
            %     it is the name of the parameter to update.
            % 
            % Returns
            % -------
            % index : integer
            %     A valid row index.
            %
            % Raises
            % ------
            % indexOutOfRange
            %     If ``row`` is an integer and it is less than 1 or more than the number of parameters.
            % nameNotRecognised
            %     If ``row`` is a parameter name and it cannot be found in the table.
            if isText(row)
                index = obj.findRowIndex(row, obj.varTable{:,1}, 'Unrecognised row name');
            else
                index = row;
                if (index < 1) || (index > obj.rowCount)
                    throw(exceptions.indexOutOfRange(sprintf('Row index out out of range 1 - %d', obj.rowCount)));
                end     
            end
        end

        function priors = validatePriors(~, priorType, mu, sigma)
            % Validates the prior types, mu and sigma variables. If mu and sigma are set on a uniform prior, 
            % they would be changed to 0 and Inf respectively and a warning issues.
            % 
            % Examples
            % --------
            % To validate a uniform prior.
            % >>> prior = obj.validatePriors('uniform', 2, 3); 
            % >>> prior
            % ans =
            %    1×3 cell array
            %      {'uniform'}    {[0]}    {[Inf]}
            %
            % Parameters
            % ----------
            % priorType : PriorTypes 
            %     The prior type to validate.
            % mu : double, default: []
            %     The mu value to validate.
            % sigma : double, default: []
            %     The sigma value to validate.
            %
            % Returns
            % -------
            % priors : cell
            %     A cell array containing the valid prior type, mu and sigma.
            %
            % Warns
            % -----
            % generic warning
            %     If ``mu`` or ``sigma`` is set on a uniform prior.
            arguments
                ~
                priorType
                mu {isscalar, mustBeNumeric}
                sigma {isscalar, mustBeNumeric}
            end
            invalidPriorsMessage = sprintf('Prior type must be a priorTypes enum or one of the following strings (%s)', ...
                                           strjoin(priorTypes.values(), ', '));
            priorType = validateOption(priorType, 'priorTypes', invalidPriorsMessage).value;

            if strcmp(priorType, priorTypes.Uniform.value)
                if mu ~= 0 
                    warning('mu cannot be %d when the prior types is uniform - resetting to 0', mu);
                end
                if sigma ~= Inf
                    warning('sigma cannot be %d when the prior types is uniform - resetting to Inf', sigma);
                end
                mu = 0;
                sigma = Inf;
            end
            priors = {priorType, mu, sigma};
        end

        function validateLimits(~, min, value, max)
            % Validates the value, minimum and maximum limit variables.
            % 
            % Examples
            % --------
            % Incorrect bounds throws exception.
            % >>> obj.validateLimits(2, 4, 1);  % throws exception as max is less than min 
            %
            % Parameters
            % ----------
            % min : double
            %     The minimum value to validate.
            % value : double
            %     The value to validate.
            % max : double
            %     The maximum value to validate.
            %
            % Returns
            % -------
            % index : integer
            %     A valid row index.
            %
            % Raises
            % ------
            % invalidValue
            %     If ``min`` is greater than ``max`` or ``value`` is outside the limits.
            arguments
                ~
                min {isscalar, mustBeNumeric}
                value {isscalar, mustBeNumeric}
                max {isscalar, mustBeNumeric}
            end
                
            if min > max
                throw(exceptions.invalidValue(sprintf('min limit %d must be less than or equal to max limit %d', min, max)));
            end

            if value < min || value > max
                throw(exceptions.invalidValue(sprintf('Parameter value %d must be within the limits %d to %d', value, min, max)));
            end
        end

    end

end
