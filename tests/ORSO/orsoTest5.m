function out = orsoTest5()
% ORSO validation Test5 - Reflectivity plus resolution...

layers = dlmread('test1.layers');

% Change the units to Å
layers(:,2) = layers(:,2) .* 1e-6;
layers(:,3) = layers(:,3) .* 1e-6;

% Read in the data.....
data = dlmread('test5.dat');

% Group the Layers
thick = layers(:,1);
sld = complex(layers(:,2),layers(:,3));
rough = layers(:,4);

% Calculate reflectivity....
[~,argmin] = min(data(:, 1));
[~, argmax] = max(data(:, 1));
q = linspace(data(argmin, 1) - 3.5 * data(argmin, 4),...
             data(argmax, 1) + 3.5 * data(argmax, 4),...
             10001);
N = size(layers,1);
ref = abelesSingle(q,N,thick,sld,rough);

% Apply resolution....
FWHM = 2 * sqrt(2 * log(2)); % FWHM for Gaussian function
resol =  (0.05 / FWHM) .* q;
ref = gaussianConvolution(q, ref, q, resol);

% Plot the comparison....
figure(1); clf
semilogy(q,ref,'k-','LineWidth',2)
hold on
plot(data(:,1),data(:,2),'r.')

% Calculate the output....
ref = interp1(q, ref, data(:, 1));
out = allClose(ref, data(:, 2), rtol=0.033);
