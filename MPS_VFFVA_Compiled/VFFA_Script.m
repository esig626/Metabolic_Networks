%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VFFVA TEST!
% Load the model you need
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preload Parameters
% Change these parameters according to your CPU: 
nCores=3;
nThreads=2;
ex='';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create an .mps Model
path = '/Path/To/Your/Files';
cd(path)
modelmps = model;
output_filename = 'model';
convertProblem(modelmps, 0, output_filename);
modelmps = 'model.mps';
cd ..

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run VFFVA!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Hold on to your seat, I am initialising the system ...')
%Set install folder of MPI
setenv('PATH', ['/opt/local/bin:', getenv('PATH')]);
[status, cmdout] = system('mpiexec --version');
%Set schedule and chunk size parameters
setenv('OMP_SCHEDUELE', ['dynamic' ',' num2str(50)])
% Get the .mps formatted model 
cd(path); 
clc
% Run VFFVA!
disp('Running Very Fast Flux Variability Analysis!')
[status,cmdout]=system(['mpirun -np ' num2str(nCores)...
' --bind-to ' num2str('none') ' -x OMP_NUM_THREADS=' num2str(nThreads)...
' ./veryfastFVA ' modelmps ' ' num2str(90)  ' ' num2str(0) ' ' ex]);
clc
disp('I am done! I will now fetch your results, hang tight...')
resultFile=[modelmps(1:end-4) 'output.csv'];
results=readtable(resultFile);
opt=[results.minFlux,results.maxFlux];
disp('All good!')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except opt model
