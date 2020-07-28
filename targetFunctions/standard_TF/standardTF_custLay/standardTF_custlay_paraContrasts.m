function         [outSsubs,...
            backgs,...
            qshifts,...
            sfs,...
            nbas,...
            nbss,...
            resols,...
            chis,...
            reflectivity,...
            Simulation,...
            shifted_data,...
            layerSlds,...
            sldProfiles,...
            allLayers,...
            allRoughs] = standardTF_custlay_paraContrasts(resample, ...
            numberOfContrasts, ...
            geometry, ...
            repeatLayers , ...
            cBacks , ...
            cShifts , ...
            cScales , ...
            cNbas , ...
            cNbss, ...
            cRes , ...
            backs , ...
            shifts , ...
            sf, ...
            nba , ...
            nbs , ...
            res , ...
            dataPresent , ...
            allData , ...
            dataLimits , ...
            simLimits , ...
            nParams , ...
            params , ...
            contrastLayers , ...
            numberOfLayers , ...
            layersDetails, ...
            problemDef_limits, ...
            fname, ...
            lang, ...
            path,...
            backsType,...
            calcSLD);

            
%coder.extrinsic('callModules');            
            
%Pre-Allocation...
backgs = zeros(numberOfContrasts,1);
qshifts = zeros(numberOfContrasts,1);
sfs = zeros(numberOfContrasts,1);
nbas = zeros(numberOfContrasts,1);
nbss = zeros(numberOfContrasts,1);
resols = zeros(numberOfContrasts,1);
allRoughs = zeros(numberOfContrasts,1);
outSsubs = zeros(numberOfContrasts,1);
chis =  zeros(numberOfContrasts,1);
allLayers = cell(numberOfContrasts,1); 
layerSlds = cell(numberOfContrasts,1);
sldProfiles = cell(numberOfContrasts,1);
shifted_data = cell(numberOfContrasts,1);

reflectivity = cell(numberOfContrasts,1);
for i = 1:numberOfContrasts
    reflectivity{i} = [1 1 ; 1 1];
end
coder.varsize('reflectivity{:}',[10000 2],[1 0]);

Simulation = cell(numberOfContrasts,1);
for i = 1:numberOfContrasts
    Simulation{i} = [1 1 ; 1 1];
end
coder.varsize('Simulation{:}',[10000 2],[1 0]);

allLayers = cell(numberOfContrasts,1);
for i = 1:numberOfContrasts
    allLayers{i} = [1 ; 1];
end
coder.varsize('allLayers{:}',[10000 3],[1 1]);


parfor i = 1:numberOfContrasts
    [backgs(i),qshifts(i),sfs(i),nbas(i),nbss(i),resols(i)] = backSort(cBacks(i),cShifts(i),cScales(i),cNbas(i),cNbss(i),cRes(i),backs,shifts,sf,nba,nbs,res);
    
    [outLayers,allRoughs(i)] = call_customLayers(params,i,fname,path,nbas(i),nbss(i));
    allLayers{i} = outLayers;
    
    [outLayers, outSsubs(i)] = groupLayers_Mod(outLayers,allRoughs(i),geometry,nbas(i),nbss(i));
    
    sldProfile = makeSLDProfiles(nbas(i),nbss(i),outLayers,outSsubs(i),repeatLayers{i});
    sldProfiles{i} = sldProfile;

    if resample(i) == 1
        layerSld = resampleLayers(sldProfile);
        layerSlds{i} = layerSld;
    else
        layerSld = outLayers;
        layerSlds{i} = layerSld;
    end
   
    shifted_dat = shiftdata(sfs(i),qshifts(i),dataPresent(i),allData{i},dataLimits{i});
    shifted_data{i} = shifted_dat;
    
    [reflect,Simul] = callReflectivity(nbas(i),nbss(i),simLimits{i},repeatLayers{i},shifted_dat,layerSld,outSsubs(i),backgs(i),resols(i),'single');
    
    [reflect,Simul,shifted_dat] = applyBackgroundCorrection(reflect,Simul,shifted_dat,backgs(i),backsType(i));
    
    reflectivity{i} = reflect;
    Simulation{i} = Simul;
    
    chis(i) = chiSquared(shifted_dat,reflect,nParams);
end


end
