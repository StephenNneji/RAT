function [allLayers, allRoughs] = loopMatalbCustlayWrapper_CustLaypoints(cBacks,cShifts,cScales,cNbas,cNbss,cRes,backs,...
     shifts,sf,nba,nbs,res,cCustFiles,numberOfContrasts,customFiles,params)
 
 %#codegen
 
 % This is a dummy function to allow source compile to proceed
 tempAllLayers = cell(numberOfContrasts,1);
 tempAllRoughs = zeros(numberOfContrasts,1);
 allLayers = cell(numberOfContrasts,1);
 allRoughs = zeros(numberOfContrasts,1);
 
 for i = 1:numberOfContrasts
     allLayers{i} = [1 , 1];    % Type def as double (size not important)
     tempAllLayers{i} = [0 0 0 0 0];
 end
 coder.varsize('tempAllLayers{:}',[10000 5],[1 1]);
 
% Dummy values to allow code generation to proceed....
for i = 1:numberOfContrasts
    tempAllLayers{i} = [0 0 0];
    tempAllRoughs(i) = 1;
end
 
 % All the following is intended to be casting from mxArray's to doubles.
 % I'm not sure if all of this is necessary, but it compiles...
 for i = 1:numberOfContrasts
     tempOut = tempAllLayers{i};
     n = [0 0];
     n = size(tempOut);
     newOut = zeros(n);
     newOut = tempOut;
     allLayers{i} = newOut;
 end
 
 %allLayers = tempAllLayers;
 allRoughs = tempAllRoughs;
 
end
 