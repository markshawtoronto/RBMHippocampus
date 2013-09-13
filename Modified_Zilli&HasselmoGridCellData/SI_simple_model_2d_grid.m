% 2d grid simulation using izhikevich simple model neuron network
% eric zilli - july 2, 2009 using simple model MATLAB code from
% Izhikevich's book modified to be a network
%
% release version 1.0. check modeldb for updates.
%
% this source code is released into the public domain
%
% This script uses experimentally recorded rat trajectory (given by the
% hafting_trajectory function) as input to multiple VCOs which can
% interfere to produce spatial firing of a postsynaptic cell which,
% hopefully, is arrayed in the classical hexagonal grid pattern.
%
% Running this script requires that an FI curve of the network has been
% generated (e.g. by SI_simple_model_FI_relation.m) for a cell or network
% with identical parameters to those used here.
%
% This is generally slow: behavioral timescales are big and neural
% timescales are tiny. Hints: using pcon=1 lets the script bypass the
% connectivity matrix which really speeds things up. If the resonant postsynaptic
% cell starts spiking too much the script will slow way down because the
% vector of spiketimes of the postsynaptic cell is a pre-allocated sparse
% vector. When adjusting postsynaptic cell parameters, start with
% short simulation to act as a sort of sanity check.
% For 2 VCO, 320 s simulations, my slow laptop takes about
% 1 hour. 160 s simulations at 30 minutes wall time is a good place to
% tweak parameters. 4 s simulations at 1.5 min to ballpark params.
%
% Presets are given for figures in the paper which all Type 2 excitable,
% but Type 1 excitability is as easy as setting the variable and providing
% the proper FI curve. Postsynaptic parameters might need to be adjusted.
% Note that parameters will need to change if the baseline frequency
% changes.
%
% Note: this is a cleaned up version of the script used to generate the
% paper figures--it is possible errors were introduced during the cleaning!
% Feel free to contact me if there is any difficulty reproducing any
% results from the manuscript.

% clear all


pcon=1; % Tweak to speed things up as recommended

simdur = 1*320e3; % ms
trajdt = 0.1; % sampling rate for returned trajectory (will interp down to dt later); ms
if simdur>4e3
  [dxs,dys,fdxs,fdys] = hafting_trajectory(simdur,trajdt);
else
  [dxs,dys,fdxs,fdys] = hafting_trajectory(4e3,trajdt);
end

% % straight line trajectory:
% fdxs = 4e-5*ones(size(fdxs));
% fdys = 0*fdys;

% % bent line trajectory:
% speedmult = 0.1;
% fdxs(1:end/2) = speedmult*.0005*dt;
% fdys(1:end/2) = 0;
% fdxs(end/2:end) = -speedmult*-.0005*dt;
% fdys(end/2:end) = 0.00001*speedmult*.0004*dt;

% % elliptical trajectory:
% tv=(0:dt:simdur)*2*2*pi/8e3; a=0.001/4; b=0.00025*1.5;
% fdxs = -dt*(- a*cos(pi/2)*sin(tv) - b*cos(tv)*sin(pi/2));
% fdys = dt*(- b*cos(pi/2)*cos(tv) - a*sin(pi/2)*sin(tv));

% if you're trying to find good postsynaptic cell parameters, you're gonna
% keep some sanity by seeding the random number generator to compare runs
% with different params:
rand('seed',7);
randn('seed',7);

% type=1 (Figures S5, S9) 1 VCO, Class 2 excitable, n=1, noise=0, filt
% type=2 (Figure 3) 2 VCOs, Class 2 excitable, n=1, noise=0, filt
% type=3 (Figure S1) 2 VCOs, Class 2 excitable, n=1, noise=0, unfilt
% type=4 (Figure 4) 2 VCOs, Class 2 excitable, n=1, noise=full, filt
% type=5 (Figures 6, 7) 2 VCOs, Class 2 excitable, synaptic, n=250, noise=full, filt
% type=6 (Figures S6, S7) 2 VCOs, Class 2 excitable, gap-junction, n=250, noise=full, filt (the preset params for this one could use work!)
% type=7 3 VCOs, Class 2 excitable, synaptic, n=250, noise=full, filt
% type=8 3 VCOs, Class 2 excitable, gap-junction, n=250, noise=full, filt (doesn't work well)
% type=9 (Figure 9) 2 VCOs, Class 2 excitable inhibitory, gap-junction, n=250, noise=full, filt
% type=10 2 VCOs, Class 2 excitable, gap-junction, n=250, noise=full, filt, nPostIn=1
% type=11 2 VCOs, Class 2 excitable, synaptic, n=5000, p=0.01
% type=12 6 VCOs, Class 2 excitable, synaptic, n=250, noise=full, filt, cheatprecession
% type=13 6 VCOs, Class 2 excitable, n=1, noise=0, filt, cheatprecession
% type=14 6 VCOs, Class 2 excitable, gap-junction, n=250, noise=full, filt, nPostIn=25, cheatprecession (doesn't work, ninhib messes up FI curve)
% type=15 sandbox
for type=[7]
  dt = .1; % ms

  nruns = 1;
  nVCOs = 3; % might be overridden below

  % generally you will want to run both models so that you have the
  % respective VCO phase differences as a visual performance measure (in
  % addition to the spatial firing itself)
  runAbstract = 1;
  runNetwork = 1;

  % generally leave these = 1, unless debugging one or the other
  runNetworkVCOs = 1;
  runNetworkBaseline = 1;

  % if true, figures will be saved to disk if a figure of the same name does
  % not already exist
  saveFigures = 0;

  % plotType = 0 means gated LIF post on top, first quarter time network or
  % abstract on bottom (Figures 3, 4, S1)
  % plotType = 1 means gated LIF post on top, non-gated res on bottom
  % (Figures 6, 7, S6, S7)
  % plotType = 2 means all three post plus abstract spatial firing for
  % nVCOs=1 (overrides abstractThr) (Figure S5, S9)
  % plotType = 3 means all three post plus xcorrs plus abstract (useful for setting params)
  plotType = 1;

  % if true, plots traces of cell #1 from each VCO and the postsynaptic
  % cells (each one gets its own figure)
  % this makes Figures 7, 12, 14)
  plotVCOsAndPosts = 1;

  % if journalChargesALotForColor, gray-scale plots will be produced
  journalChargesALotForColor=1;

  % if 1, plots the velocity component on the VCO 1 phase diff plot for
  % plotType==0
  plotVelocity = 0;

  % for plotType = 0, this controls whether the bottom plot is the abstract
  % model's output or the first 1/4 of the simulation for the network model
  showAbstractFig = 0;

  % if =1 the postsynaptic cells (LIF) are quiescent at rest (I=0), if =0 the
  % first postsynaptic cell (simple model) actively (tonically) spikes at
  % I=0 and its inputs become inhibitory
  stablePost = 1;

  %% baseline synapses onto the post cell are this many times stronger than
  %% VCO syns (not counting that they are also nVCO times stronger)
  baselineMult = 1;

  % live plot of network spatial activity
  runningPlot = 0;

  % number of cells in each VCO the postsynaptic cell receives input from
  nPostIn = 1;
  
  % if cheatHDRectification=1, VCO inputs to the grid cell are
  % suppressed if the current head direction is not within 180 degrees of
  % the VCO's preferred direction
  cheatHDRectification=0;

  % if opponentVCOInhibition=1, cells 1:nInhib of each VCO get an input
  % current of magnitude opponentInhibAmp if the animal's direction is not
  % within 90 degrees of its preferred direction
  opponentVCOInhibition = 0;
  nInhib = 0;
  opponentInhibAmp = 0;
  
  % if non-empty, we'll save the VCO spike times to disk (for efficiently
  % testing postsynaptic parameters, can just generate the VCOs themselves
  % one and re-use the spikes)
    saveVCOSpikesFileName = '';
%   saveVCOSpikesFileName = '15s_res2_fi_spikes_dur4s.txt';

  % if non-empty, we'll load VCO spike times from disk from this file
  loadVCOSpikesFileName = '';
%     loadVCOSpikesFileName = 'testspikes.txt';

  % connectivity matrix
  % if not-empty, this filed will be loaded to provide the connectivity
  % matrix C
  loadC = '';
  %   % to make a C:
  %   C = (rand(ncells)<pcon)>0;
  %   C = sparse(C-eye(size(C))>0); % remove autoconnections
  %   save C_simple_nX_pX.mat C

  if runningPlot
    runh = figure;
  end

  %%% ANY OF THESE may be replaced in the type-unique parameters below

  % if baselineLIFGating = 1, the LIF only gets active VCO inputs and only during the
  % basegateDur milliseconds after each baseline VCO spike
  baselineLIFGating = 1;
  basegateDur = 10; % ms
  abstractThr=3;
  simprefix = '';
  ncells = 1;
  pcon = 1;
  useNoise = 0;
  commonNoiseSTD = 0;
  useFilteredTrajectory = 1;

  if type==1
    simprefix = sprintf('filthaftingtraj_n1_0noise_1vco_%s',datestr(now,'mmmmdd'));
    excitationClass = 2;
    useGapJunctions = 1;
    showAbstractFig = 1;
    plotType = 2;
    ncells = 1;
    nPostIn = ncells;
    pcon = 1;
    useNoise = 0;
    useFilteredTrajectory = 0;
    load simple_model_RS2_FI_Jan09_n1.mat;
    uniqueNoiseSTD = 100*useNoise;
    g = 0;
    nVCOs = 1;
    abstractThr=1.65*nVCOs;
    % shared grid variables: freq = baselineFreq+beta*speed*(cos(prefHD-curHD));
    baselineFreq = freqs(-1+round(length(freqs)/2)); % Hz
    beta = 2; % Hz/(m/s)

    % non-gated lif params:
    tau = 5; % ms
    membraneDecay = exp(-dt/tau);
    postWeight=0.91;
    baseWeight = 1*postWeight;

    % gated lif params:
    basegateDur = 25; % ms
    tau = 20; % (msec)
    gatedmembraneDecay = exp(-dt/tau);
    weightMult = 1.3;
    gatedpostWeight = weightMult/ncells/nVCOs;

    % non-gated res params:
    weightMult = 1;
    respostWeight=2;
    resbaseWeight = 6;
    postc = -0.01;
    postw = baselineFreq*2*pi/1000; % Hz

  elseif type==2
    simprefix = sprintf('filthaftingtraj_n1_0noise_01res%s',datestr(now,'mmmmdd'));
    excitationClass = 2;
    useGapJunctions = 0;
    showAbstractFig = 1;
    plotType = 0;
    plotVelocity = 1;
    ncells = 1;
    nPostIn = ncells;
    pcon = 0;
    useNoise = 0;
    uniqueNoiseSTD = 100*useNoise;
    useFilteredTrajectory = 1;
    load simple_model_RS2_FI_Jan09_n1.mat;
    g = 0;
    % shared grid variables: freq = baselineFreq+beta*speed*(cos(prefHD-curHD));
    baselineFreq = freqs(round(length(freqs)/2)); % Hz
    beta = 2; % Hz/(m/s)

    % non-gated lif params:
    tau = 40*7.8989/baselineFreq; 35; % ms
    membraneDecay = exp(-dt/tau);
    postWeight = 0.14;
    baseWeight = .8;

    % gated lif params:
    basegateDur = 40; % ms
    tau = 25; % (msec)
    gatedmembraneDecay = exp(-dt/tau);
    weightMult = 1.4;
    gatedpostWeight = weightMult/ncells/nVCOs;

    % non-gated res params:
    weightMult = 1;
    respostWeight=1.4; weightMult*1/ncells/(nVCOs);
    resbaseWeight = 6; 5; 2*postWeight;
    postc = -0.01;
    postw = baselineFreq*2*pi/1000; % Hz (times 2pi for angular freq, /1000 to convert to /s
  elseif type==3
    simprefix = sprintf('unfilthaftingtraj_n1_0noise_01res_%s',datestr(now,'mmmmdd'));
    excitationClass = 2;
    useGapJunctions = 0;
    showAbstractFig = 0;
    plotType = 0;
    ncells = 1;
    nPostIn = ncells;
    pcon = 0;
    useNoise = 0;
    useFilteredTrajectory = 0;
    load simple_model_RS2_FI_Jan09_n1.mat;
    uniqueNoiseSTD = 100*useNoise;
    g = 0;
    % shared grid variables: freq = baselineFreq+beta*speed*(cos(prefHD-curHD));
    baselineFreq = freqs(round(length(freqs)/2)); % Hz
    beta = 2; % Hz/(m/s)
    % non-gated lif params:
    tau = 35; % ms
    membraneDecay = exp(-dt/tau);
    postWeight = 0.16;
    baseWeight = .8;

    % gated lif params:
    basegateDur = 40; % ms
    tau = 25; % (msec)
    gatedmembraneDecay = exp(-dt/tau);
    weightMult = 1.4;
    gatedpostWeight = weightMult/ncells/nVCOs;

    % non-gated res params:
    weightMult = 1;
    respostWeight=1.5; weightMult*1/ncells/(nVCOs);
    resbaseWeight = 5; 2*postWeight;
    postc = -0.0075;
    postw = baselineFreq*2*pi/1000; % Hz
  elseif type==4
    simprefix = sprintf('filthaftingtraj_n1_1noise_01res_%s',datestr(now,'mmmmdd'));
    excitationClass = 2;
    useGapJunctions = 0;
    showAbstractFig = 0;
    plotType = 0;
    ncells = 1;
    nPostIn = ncells;
    pcon = 0;
    useNoise = 1;
    useFilteredTrajectory = 1;
    basegateDur = 40;
    load simple_model_RS2n_FI_Jan09_n1;
    uniqueNoiseSTD = 100*useNoise;
    % shared grid variables: freq = baselineFreq+beta*speed*(cos(prefHD-curHD));
    baselineFreq = freqs(round(length(freqs)/2)); % Hz
    beta = 2; % Hz/(m/s)
    %% Use whatever parameters you like here, you won't get a grid cell!
    % non-gated lif params:
    tau = 35; % ms
    membraneDecay = exp(-dt/tau);
    postWeight = 0.16;
    baseWeight = .8;

    % gated lif params:
    basegateDur = 40; % ms
    tau = 25; % (msec)
    gatedmembraneDecay = exp(-dt/tau);
    weightMult = 1.4;
    gatedpostWeight = weightMult/ncells/nVCOs;

    % non-gated res params:
    weightMult = 1;
    respostWeight=1.5; weightMult*1/ncells/(nVCOs);
    resbaseWeight = 5; 2*postWeight;
    postc = -0.0075;
    postw = baselineFreq*2*pi/1000; % Hz
  elseif type==5
    simprefix = sprintf('filthaftingtraj_n250_1noise_01res_%s',datestr(now,'mmmmdd'));
    excitationClass = 2;
    useGapJunctions = 0;
    showAbstractFig = 0;
    plotType = 1;
    ncells = 250;
    nPostIn = ncells;
    pcon = 1;
    useNoise = 1;
    useFilteredTrajectory = 1;
    basegateDur = 1;
    load simple_model_RS2sn_FI_Jan09_n250.mat;
    uniqueNoiseSTD = 100*useNoise;
    g = 256;
    baselineFreq = freqs(round(2+length(freqs)/2)); % Hz
    beta = 2; % Hz/(m/s)

    % non-gated lif params:
    tau = 5; % ms
    weightMult = 0.7;
    membraneDecay = exp(-dt/tau);
    postWeight=weightMult*1/ncells/(2+nVCOs);
    baseWeight = 2*postWeight;

    % gated lif params:
    basegateDur = 1; 50; % ms
    tau = 5; % (msec)
    gatedmembraneDecay = exp(-dt/tau);
    weightMult = 1.5; 1.3;
    gatedpostWeight = 0.0016;

    % non-gated res params:
    weightMult = 1;
    respostWeight=0.002;
    resbaseWeight = 0.0024;
    postc = -0.01;
    postw = baselineFreq*2*pi/1000; % Hz
  elseif type==6
    simprefix = sprintf('filthaftingtraj_n250_1noise_02res_%s',datestr(now,'mmmmdd'));
    excitationClass = 2;
    useGapJunctions = 1;
    showAbstractFig = 0;
    plotType = 1;
    ncells = 250;
    nPostIn = ncells;
    pcon = 1;
    useNoise = 1;
    useFilteredTrajectory = 1;
    load simple_model_RS2gn_FI_Jan09_n250.mat;
    uniqueNoiseSTD = 100*useNoise;
    g = 0.1;

    % shared grid variables: freq = baselineFreq+beta*speed*(cos(prefHD-curHD));
    baselineFreq = freqs(1+round(length(freqs)/2)); % Hz
    beta = 2; % Hz/(m/s)

    % non-gated lif params:
    tau = 20; % ms
    membraneDecay = exp(-dt/tau);
    postWeight=.19/250;
    baseWeight = .8/250;

    % gated lif params:
    % I&F params
    basegateDur = 50; % ms
    tau = 50; % (msec)
    gatedmembraneDecay = exp(-dt/tau);
    weightMult = 1.5;
    gatedpostWeight = weightMult/ncells/nVCOs;

    % non-gated res params:
    weightMult = 8;
    respostWeight=weightMult*1/ncells/(nVCOs);
    resbaseWeight = 16*postWeight;
    postc = -0.01;
    postw = baselineFreq*2*pi/1000; % Hz
  elseif type==7
    rand('seed',1);
    randn('seed',1);
    simprefix = sprintf('filthaftingtraj_n250_1noise_3vco_%s',datestr(now,'mmmmdd'));
    excitationClass = 2;
    useGapJunctions = 0;
    showAbstractFig = 1;
   % showAbstractFig = 0;
   % plotType = 2;
    nVCOs = 3;
    ncells = 250;
    nPostIn = ncells;
    pcon = 1;
    useNoise = 1;
    useFilteredTrajectory = 1;
    load simple_model_RS2sn_FI_Jan09_n250.mat;
    uniqueNoiseSTD = 100*useNoise;
    g = 256;
    baselineFreq = freqs(round(-1+length(freqs)/2)); % Hz
    beta = 2; % Hz/(m/s)

    % non-gated lif params:
    tau = 5; % ms
    weightMult = 1; 0.8;
    membraneDecay = exp(-dt/tau);
    postWeight=weightMult*1/ncells/(2+nVCOs);
    baseWeight = 1.5*postWeight; 2*postWeight;

    % gated lif params:
    basegateDur = .05; 1; 50; % ms
    tau = 5; % (msec)
    gatedmembraneDecay = exp(-dt/tau);
    gatedpostWeight = 0.0014;

    % non-gated res params:
    weightMult = 7;
    respostWeight=0.0015; 0.002;
    resbaseWeight = 0.003; 0.0024;
    postc = -0.01;
    postw = baselineFreq*2*pi/1000; % Hz
  elseif type==8
    simprefix = sprintf('filthaftingtraj_n250_1noise_3vco_%s',datestr(now,'mmmmdd'));
    excitationClass = 2;
    useGapJunctions = 1;
    showAbstractFig = 0;
    plotType = 3;
    nVCOs = 3;
    ncells = 250;
    nPostIn = ncells;
    pcon = 1;
    useNoise = 1;
    useFilteredTrajectory = 1;
    load simple_model_RS2gn_FI_Jan09_n250.mat;
    uniqueNoiseSTD = 100*useNoise;
    g = 0.1;
    baselineFreq = freqs(1+round(length(freqs)/2)); % Hz
    beta = 2; % Hz/(m/s)

    % non-gated lif params:
    tau = 30; % ms
    membraneDecay = exp(-dt/tau);
    postWeight=.15/250;
    baseWeight = .8/250;

    % gated lif params:
    % I&F params
    basegateDur = 50; % ms
    tau = 50; % (msec)
    gatedmembraneDecay = exp(-dt/tau);
    weightMult = 1.5;
    gatedpostWeight = weightMult/ncells/nVCOs;

    % non-gated res params:
    weightMult = 9;
    respostWeight=weightMult*1/ncells/(nVCOs);
    resbaseWeight = 16*postWeight;
    postc = -0.01;
    postw = baselineFreq*2*pi/1000; % Hz
  elseif type==9
    simprefix = sprintf('filthaftingtraj_n250_inhib_1noise_02res_%s',datestr(now,'mmmmdd'));
    excitationClass = 2;
    useGapJunctions = 1;
    showAbstractFig = 0;
    plotType = 3;
    ncells = 250;
    nPostIn = ncells;
    stablePost = 0;
    pcon = 1;
    useNoise = 1;
    useFilteredTrajectory = 1;
    load simple_model_RS2gn_FI_Jan09_n250.mat;
    uniqueNoiseSTD = 100*useNoise;
    g = 0.1;

    % shared grid variables: freq = baselineFreq+beta*speed*(cos(prefHD-curHD));
    baselineFreq = freqs(-2+round(length(freqs)/2)); % Hz
    beta = 2; % Hz/(m/s)

    % simple model postsynaptic params:
    tau = 20; % ms
    membraneDecay = exp(-dt/tau);
    postWeight=1000*-.19/250;
    baseWeight = 1000*-.8/250;

    % gated lif params:
    % I&F params
    basegateDur = 50; % ms
    tau = 50; % (msec)
    gatedmembraneDecay = exp(-dt/tau);
    weightMult = 1.5;
    gatedpostWeight = weightMult/ncells/nVCOs;

    % non-gated res params:
    weightMult = 8;
    respostWeight=weightMult*1/ncells/(nVCOs);
    resbaseWeight = 16*postWeight;
    postc = -0.01;
    postw = baselineFreq*2*pi/1000; % Hz
  elseif type==10
    simprefix = sprintf('filthaftingtraj_n250_1noise_02res_%s',datestr(now,'mmmmdd'));
    excitationClass = 2;
    useGapJunctions = 1;
    showAbstractFig = 0;
    plotType = 1;
    ncells = 250;
    nPostIn = 1;
    pcon = 1;
    useNoise = 1;
    useFilteredTrajectory = 1;
    load simple_model_RS2gn_FI_Jan09_n250.mat;
    uniqueNoiseSTD = 100*useNoise;
    g = 0.1;

    % shared grid variables: freq = baselineFreq+beta*speed*(cos(prefHD-curHD));
    baselineFreq = freqs(2+round(length(freqs)/2)); % Hz
    beta = 2; % Hz/(m/s)

    % non-gated lif params:
    tau = 40; % ms
    membraneDecay = exp(-dt/tau);
    postWeight = 0.14;
    baseWeight = .8;

    % gated lif params:
    basegateDur = 30; % ms
    tau = 20; % (msec)
    gatedmembraneDecay = exp(-dt/tau);
    %     weightMult = 1.4;
    gatedpostWeight = 0.9;% weightMult/nVCOs;

    % non-gated res params:
    weightMult = 1;
    respostWeight=1.5;
    resbaseWeight = 5.5;
    postc = -0.01;
    postw = baselineFreq*2*pi/1000; % Hz (times 2pi for angular freq, /1000 to convert to /s
  elseif type==11
    simprefix = sprintf('filthaftingtraj_n5000_p01_1noise_02res_%s',datestr(now,'mmmmdd'));
    loadVCOSpikesFileName = 'fig_RS2sn_filthaftingtraj_n5000_p01_1noise_02res_July13_2D_spikes_0a.txt';
    excitationClass = 2;
    useGapJunctions = 0;
    showAbstractFig = 0;
    plotType = 3;
    ncells = 5000;
    nPostIn = 5000;
    pcon = 0.01;
    useNoise = 1;
    useFilteredTrajectory = 1;
    load simple_model_RS2sn_FI_Jul11_n5000.mat
    uniqueNoiseSTD = 100*useNoise;
    g = 256;
    loadC = 'fig_RS2sn_filthaftingtraj_n5000_p01_1noise_02res_July13_C_0a.mat';

    % shared grid variables: freq = baselineFreq+beta*speed*(cos(prefHD-curHD));
    baselineFreq = freqs(2+round(length(freqs)/2)); % Hz
    beta = 2; % Hz/(m/s)

    % non-gated lif params:
    tau = 12; % ms
    membraneDecay = exp(-dt/tau);
    postWeight = 0.00018;
    baseWeight =  .00019;

    % gated lif params:
    basegateDur = 3; % ms
    tau = 5; % (msec)
    gatedmembraneDecay = exp(-dt/tau);
    %     weightMult = 1.4;
    gatedpostWeight = 0.00035;

    % non-gated res params:
    weightMult = .001;
    respostWeight=1e-3;
    resbaseWeight = 1e-3;
    postc = -0.01;
    postw = baselineFreq*2*pi/1000; % Hz (times 2pi for angular freq, /1000 to convert to /s
  elseif type==12
    rand('seed',1);
    randn('seed',1);
    simprefix = sprintf('filthaftingtraj_6vcoprec_n250_1noise_6vco_prec_%s',datestr(now,'mmmmdd'));
    excitationClass = 2;
    useGapJunctions = 0;
    showAbstractFig = 0;
    plotType = 2;
    nVCOs = 6;
    ncells = 250;
    nPostIn = ncells;
    cheatHDRectification = 1;
    pcon = 1;
    useNoise = 1;
    useFilteredTrajectory = 1;
    load simple_model_RS2sn_FI_Jan09_n250.mat;
    uniqueNoiseSTD = 100*useNoise;
    g = 256;
    baselineFreq = freqs(round(-1+length(freqs)/2)); % Hz
    beta = 2; % Hz/(m/s)
    
    phaseLocked = 0;
    if phaseLocked==0
      % post-stronger (precessing) params:
      % non-gated lif params:
      tau = 5; % ms
      weightMult = 1.5; 1; 0.8;
      membraneDecay = exp(-dt/tau);
      postWeight=1.4*weightMult*1/ncells/(2+nVCOs);
      baseWeight = 0.5*postWeight; 1.5*postWeight; 2*postWeight;

      % gated lif params:
      basegateDur = .05; 1; 50; % ms
      tau = 1.2*5; % (msec)
      gatedmembraneDecay = exp(-dt/tau);
      gatedpostWeight = 0.0015; 0.0014;

      % non-gated res params:
      weightMult = 7;
      respostWeight= .0035; 1.5*0.0015; 0.002;
      resbaseWeight = .002; 0.003; 0.0024;
      postc = -0.01;
      postw = baselineFreq*2*pi/1000; % Hz
    elseif phaseLocked==1
      % baseline-stronger (non-precessing) params:
      % non-gated lif params:
      tau = 5; % ms
      weightMult = 1.25; 1; 0.8;
      membraneDecay = exp(-dt/tau);
      postWeight=1*weightMult*1/ncells/(2+nVCOs);
      baseWeight = 2.5*postWeight;

      % gated lif params:
      basegateDur = .05; 1; 50; % ms
      tau = 1.2*5; % (msec)
      gatedmembraneDecay = exp(-dt/tau);
      gatedpostWeight = 1.2*0.0014;

      % non-gated res params:
      weightMult = 7;
      respostWeight= .3*0.003; 0.002;
      resbaseWeight = 1.6*0.003; 0.0024;
      postc = -0.01;
      postw = baselineFreq*2*pi/1000; % Hz
    end

    % % params for no-cheat test:
%     % non-gated lif params:
%     tau = 5; % ms
%     weightMult = .5*1.5; 1; 0.8;
%     membraneDecay = exp(-dt/tau);
%     postWeight=weightMult*1/ncells/(2+nVCOs);
%     baseWeight = 1.5*postWeight; 2*postWeight;
% 
%     % gated lif params:
%     basegateDur = .05; 1; 50; % ms
%     tau = 1.2*5; % (msec)
%     gatedmembraneDecay = exp(-dt/tau);
%     gatedpostWeight = .5*0.0014;
% 
%     % non-gated res params:
%     weightMult = 7;
%     respostWeight= .5*1.5*0.0015; 0.002;
%     resbaseWeight = 0.003; 0.0024;
%     postc = -0.01;
%     postw = baselineFreq*2*pi/1000; % Hz
  elseif type==13
    simprefix = sprintf('filthaftingtraj_6vcoprec_n1_0noise_01res%s',datestr(now,'mmmmdd'));
    excitationClass = 2;
    useGapJunctions = 0;
    showAbstractFig = 1;
    plotType = 2;
    plotVelocity = 1;
    nVCOs = 6;
    ncells = 1;
    nPostIn = ncells;
    cheatHDRectification = 1;
    pcon = 0;
    useNoise = 0;
    uniqueNoiseSTD = 100*useNoise;
    useFilteredTrajectory = 1;
    load simple_model_RS2_FI_Jan09_n1.mat;
    g = 0;
    % shared grid variables: freq = baselineFreq+beta*speed*(cos(prefHD-curHD));
    baselineFreq = freqs(round(length(freqs)/2)); % Hz
    beta = 2; % Hz/(m/s)

    % non-gated lif params:
    tau = 10; 40*7.8989/baselineFreq; 35; % ms
    membraneDecay = exp(-dt/tau);
    postWeight = .5; 0.14;
    baseWeight = .1; .2; .8;

    % gated lif params:
    basegateDur = 40; % ms
    tau = 25; % (msec)
    gatedmembraneDecay = exp(-dt/tau);
    weightMult = 3; 2;
    gatedpostWeight = weightMult/ncells/nVCOs;

    % non-gated res params:
    weightMult = 1;
    respostWeight= 3; 3; 4; 1.4; weightMult*1/ncells/(nVCOs);
    resbaseWeight = 1.5; 2; 6; 5; 2*postWeight;
    postc = -0.01;
    postw = baselineFreq*2*pi/1000; % Hz (times 2pi for angular freq, /1000 to convert to /s
  elseif type==14
    simprefix = sprintf('filthaftingtraj_6vcoprec_post25_n250_1noise_02res_%s',datestr(now,'mmmmdd'));
    excitationClass = 2;
    useGapJunctions = 1;
    showAbstractFig = 0;
    plotType = 2;
    ncells = 250;
    nVCOs = 6;
    nPostIn = 1;
    pcon = 1;
    useNoise = 1;
    useFilteredTrajectory = 1;
    load simple_model_RS2gn_FI_Jan09_n250.mat;
    
    opponentVCOInhibition = 1;
    nInhib = nPostIn;
    opponentInhibAmp = -max(currents);
    uniqueNoiseSTD = 100*useNoise;
    g = 0.1;

    % shared grid variables: freq = baselineFreq+beta*speed*(cos(prefHD-curHD));
    baselineFreq = freqs(2+round(length(freqs)/2)); % Hz
    beta = 2; % Hz/(m/s)

    % non-gated lif params:
    tau = 12; % ms
    membraneDecay = exp(-dt/tau);
    postWeight = 0.15/10;
    baseWeight = .07/10;

    % gated lif params:
    basegateDur = 30; % ms
    tau = 10; % (msec)
    gatedmembraneDecay = exp(-dt/tau);
    %     weightMult = 1.4;
    gatedpostWeight = 0.05;% weightMult/nVCOs;

    % non-gated res params:
    weightMult = 1;
    respostWeight=.14;
    resbaseWeight = 0.01;
    postc = -0.01;
    postw = baselineFreq*2*pi/1000; % Hz (times 2pi for angular freq, /1000 to convert to /s
  elseif type==15
    simprefix = sprintf('filthaftingtraj_n1_0noise_01res%s',datestr(now,'mmmmdd'));
    excitationClass = 2;
    useGapJunctions = 0;
    showAbstractFig = 0;
    plotType = 1;
    plotVelocity = 1;
    ncells = 1;
    nPostIn = ncells;
    pcon = 0;
    useNoise = 0;
    uniqueNoiseSTD = 100*useNoise;
    useFilteredTrajectory = 1;
    load simple_model_RS2_FI_Jan09_n1.mat;
    g = 0;
    % shared grid variables: freq = baselineFreq+beta*speed*(cos(prefHD-curHD));
    baselineFreq = freqs(round(length(freqs)/2)); % Hz
    beta = 2; % Hz/(m/s)

    % non-gated lif params:
    tau = 40*7.8989/baselineFreq; 35; % ms
    membraneDecay = exp(-dt/tau);
    postWeight = 1.2*0.14;
    baseWeight = .8;

    % gated lif params:
    basegateDur = 40; % ms
    tau = 25; % (msec)
    gatedmembraneDecay = exp(-dt/tau);
    weightMult = 1.4;
    gatedpostWeight = weightMult/ncells/nVCOs;

    % non-gated res params:
    weightMult = 1;
    respostWeight=3*1.4; weightMult*1/ncells/(nVCOs);
    resbaseWeight = .5*6; 5; 2*postWeight;
    postc = -0.01;
    postw = baselineFreq*2*pi/1000; % Hz (times 2pi for angular freq, /1000 to convert to /s
  end

  typePrefix = sprintf('%d',excitationClass);
  if ncells>1
    if useGapJunctions
      typePrefix = [typePrefix 'g'];
    else
      typePrefix = [typePrefix 's'];
    end
  end
  if useNoise
    typePrefix = [typePrefix 'n'];
  end

  Cf=100; vr=-60; vt=-40; k=0.7; % parameters used for RS
  a=0.03; c=-50; d=100; % neocortical pyramidal neurons
  vpeak=35; % spike cutoff
  if abs(excitationClass)==1
    b=-2;
  elseif abs(excitationClass)==2
    b = 2;
  end

  if excitationClass<0
    d = 60;
  end

  if ~isempty(loadC)
    load(loadC)
  else
    % connectivity matrix
    if pcon==1
      C = 1;
    elseif pcon==0
      C = 0;
    else
      C = (rand(ncells)<pcon)>0;
      C = sparse(C.*(1-eye(ncells))); % remove autoconnections
    end
  end

  stabilities = zeros(1,nruns);
  stabilities2 = zeros(1,nruns);
  for run=1:nruns
    % activity of postsynaptic I&F
    npost = 3;
    post = zeros(npost, round(simdur/dt)+1);
    if stablePost==0
      postu = zeros(npost, round(simdur/dt)+1);
      postinhib = zeros(npost, round(simdur/dt)+1);
    end
    vpost = zeros(nVCOs,round(simdur/dt)+1);
    vpostb = zeros(1,round(simdur/dt)+1);

    % spike times of postsynaptic I&F
    postspikes = spalloc(1, round(simdur/dt)+1, 20*simdur); % expect firing at, say, 20 Hz

    t = 0; % current time in simulation
    v = vr*ones(ncells,nVCOs); % vr*rand(ncells,1);
    u = 0*v; % initial values
    vb = vr*ones(ncells,1); % vr*rand(ncells,1);
    ub = 0*vb; % initial values

    % save state of one cell over simulation:
    state = zeros(1+nVCOs,simdur/dt);

    % binary array indicating which cells fire on each time steps
    spikes = zeros(ncells, nVCOs);
    spikesb = zeros(ncells, 1);

    % not keeping full array anymore, but still want to know times any cell
    % spiked:
    clear VCOSpikeTimes;
    VCOSpikeTimes{nVCOs} = []; % implicitly create this struct/class/whatever thing
    BaseSpikeTimes = [];
    VCOinds = zeros(1,nVCOs);
    Baseind = 0;

    %% abstract grid model variables:
    % dbasePhi/dt = baselineFreq
    % dVCOPhi/dt = (baselineFreq+beta*speed*(cos(prefHD-curHD)))
    basePhi = 0; % baseline oscillator phase variable
    VCOPhi = zeros(1,nVCOs); % VCO phase variable(s)
    prefHDs = [0 2*pi/3 4*pi/3 pi/3 pi 5*pi/3]; % VCO preferred direction(s), (radians)
    if nVCOs>length(prefHDs)
      prefHDs = repmat(prefHDs,1,ceil(nVCOs/length(prefHDs)));
    end
    prefHDs = prefHDs(1:nVCOs);

    abstractGrid = zeros(1,round(simdur/dt));
    basehist = zeros(1,round(simdur/dt));
    vcohist = zeros(nVCOs,round(simdur/dt));

    x = 0; % m
    dx = 0; % m
    y = 0; % m
    dy = 0; % m

    commonNoise = 0;

    pos = zeros(2,round(simdur/dt));
    baseI = currents(find(freqs==baselineFreq));
    I = baseI;
    ISIs = zeros(1,nVCOs);

    speed = 0; % m/s

    basegateCount = 0;

    buffer = zeros(1,nVCOs+1);

    if ~isempty(loadVCOSpikesFileName)
      loadedSpikes = textread(loadVCOSpikesFileName,'%d');
      loadedSpikes = reshape(loadedSpikes,nVCOs+1,[])';
    end

    tic
    while t<simdur-2*dt
      t = t+dt; % advance to next time step
      tind = 1+round(t/dt);

      if mod(t,round(simdur)/100)<=dt
        disp(sprintf('t = %g, %g elapsed',t,toc));
      end

      % move virtual rat: pregenerated trajectory in variables dxs dys
      if useFilteredTrajectory
        %% filtered trajectory:
        dx = fdxs(tind);
        dy = fdys(tind);
      else
        %% unfiltered trajectory:
        dx = dxs(tind);
        dy = dys(tind);
      end
      x = x+dx;
      y = y+dy;


      speed = abs(sqrt((dx^2 + dy^2)/(dt/1000)^2)); % (m/s); abs because we include the direction via cosine later
      pos(:,tind) = [x; y];
      curHD = atan2(dy,dx);

      if runAbstract
        %% abstract grid model: (divide dt by 1000 because freqs are in Hz
        %% but dt is in ms)
        basePhi = basePhi + dt/1000*2*pi*baselineFreq;
        VCOAngularFreqs = 2*pi*(baselineFreq+beta*speed*(cos(prefHDs-curHD)));
        VCOPhi = VCOPhi + dt/1000*VCOAngularFreqs;
        abstractGrid(tind) = nVCOs*cos(basePhi)+sum(cos(VCOPhi));
        basehist(tind) = basePhi;
        vcohist(:,tind) = VCOPhi;
      end

      if isempty(loadVCOSpikesFileName)
        if runNetwork
          if commonNoiseSTD
            commonNoise = commonNoiseSTD*randn;
          end
          %% active VCO networks
          if runNetworkVCOs
            for vco=1:nVCOs
              oldv = v(:,vco);
              oldu = u(:,vco);
              % much faster if we interpolate the FI curve ourself, though we do lose
              % the ability to extrapolate, which will cause an error
              % here with one of lowind or highind being empty
              desiredFreq = baselineFreq + beta*speed*(cos(prefHDs(vco)-curHD));
              lowind = find(freqs<desiredFreq,1,'last');
              highind = find(freqs>desiredFreq,1,'first');
              proportion = (desiredFreq-freqs(lowind))/(freqs(highind)-freqs(lowind));
              I = currents(lowind) + proportion*(currents(highind)-currents(lowind));
              
              inhibI = zeros(size(v(:,vco)));
              if opponentVCOInhibition
                inhibI(1:nInhib) = opponentInhibAmp.*(cos(prefHDs(vco)-curHD)<0);
              end
              
              if useGapJunctions
                if pcon==1
                  % if all-to-all connected, the gap-junctions are equivalent
                  % to sum(v(:,vco))-ncells*v(:,vco)
                  v(:,vco) = oldv + dt*(k*(oldv-vr).*(oldv-vt) - oldu + I + inhibI + commonNoise + uniqueNoiseSTD*randn(ncells,1) + g*(sum(v(:,vco))-ncells*v(:,vco)))/Cf;
                else
                  % make matrix of differences in membrane potential
                  vdiffs = C.*(repmat(oldv',ncells,1) - repmat(oldv,1,ncells));
                  v(:,vco) = oldv + dt*(k*(oldv-vr).*(oldv-vt) - oldu + I + inhibI + commonNoise + uniqueNoiseSTD*randn(ncells,1) + g*sum(vdiffs,2))/Cf;
                  %% slowest way:
                  % v(:,vco) = oldv + dt*(k*(oldv-vr).*(oldv-vt) - oldu + I + inhibI + commonNoise + uniqueNoiseSTD*randn(ncells,1) + g*sum(C.*vdiffs,2))/Cf;
                end
              else % using synaptic connections:
                if pcon==1
                  v(:,vco) = oldv + dt*(k*(oldv-vr).*(oldv-vt) - oldu + I + inhibI + commonNoise + uniqueNoiseSTD*randn(ncells,1) + g*(sum(spikes(:,vco),1)-spikes(:,vco)))/Cf;
                else
                  v(:,vco) = oldv + dt*(k*(oldv-vr).*(oldv-vt) - oldu + I + inhibI + commonNoise + uniqueNoiseSTD*randn(ncells,1) + g*(C*double(spikes(:,vco))))/Cf;
                end
              end
              u(:,vco) = oldu + dt*a*(b*(oldv-vr)-oldu);

              % save and reset spikes when v>=vpeak
              spikes(:,vco) = v(:,vco)>=vpeak;
              if excitationClass>0
                u(v(:,vco)>=vpeak,vco) = u(v(:,vco)>=vpeak,vco)+d;
              else
                u(v(:,vco)>=vpeak,vco) = d;
              end
              oldv(v(:,vco)>=vpeak) = vpeak;
              v(v(:,vco)>=vpeak,vco) = c;
              if any(spikes(:,vco))
                if mod(VCOinds(vco),1000)==0
                  VCOSpikeTimes{vco}(length(VCOSpikeTimes{vco})+1000) = 0;
                end
                VCOinds(vco) = VCOinds(vco)+1;
                VCOSpikeTimes{vco}(VCOinds(vco)) = t;
              end
            end
          end
          spikesum = sum(spikes(1:nPostIn,:),1);

          if runNetworkBaseline
            %% baseline VCO network
            oldvb = vb;
            oldub = ub;
            if useGapJunctions
              if pcon==1
                % if all-to-all connected, the gap-junctions are equivalent
                % to sum(vb)-ncells*vb
                vb = oldvb + dt*(k*(oldvb-vr).*(oldvb-vt) - oldub + baseI + commonNoise + uniqueNoiseSTD*randn(ncells,1) + g*(sum(vb)-ncells*vb))/Cf;
              else
                % make matrix of differences in membrane potential
                vbdiffs = C.*(repmat(oldvb',ncells,1) - repmat(oldvb,1,ncells));
                vb = oldvb + dt*(k*(oldvb-vr).*(oldvb-vt) - oldub + baseI + commonNoise + uniqueNoiseSTD*randn(ncells,1) + g*sum(vbdiffs,2))/Cf;
                %% slowest way:
                % vb = oldvb + dt*(k*(oldvb-vr).*(oldvb-vt) - oldub + baseI + commonNoise + uniqueNoiseSTD*randn(ncells,1) + g*sum(C.*vbdiffs,2))/Cf;
              end
            else
              if pcon==1
                vb = oldvb + dt*(k*(oldvb-vr).*(oldvb-vt) - oldub + baseI + commonNoise + uniqueNoiseSTD*randn(ncells,1) + g*(sum(spikesb,1)-spikesb))/Cf;
              else
                vb = oldvb + dt*(k*(oldvb-vr).*(oldvb-vt) - oldub + baseI + commonNoise + uniqueNoiseSTD*randn(ncells,1) + g*(C*double(spikesb)))/Cf;
              end
            end
            ub = oldub + dt*a*(b*(oldvb-vr)-oldub);

            % save and reset spikes when v>=vpeak
            spikesb = vb>=vpeak;
            %       spikesb(:,tind) = vb>=vpeak;
            spikesumb = sum(spikesb(1:nPostIn));
            if excitationClass>0
              ub(vb>=vpeak) = ub(vb>=vpeak)+d;
            else
              ub(vb>=vpeak) = d;
            end
            oldvb(vb>=vpeak) = vpeak;
            vb(vb>=vpeak) = c;
%             if any(spikesb)
%               if mod(Baseind,1000)==0
%                 BaseSpikeTimes(length(BaseSpikeTimes)+1000) = 0;
%               end
%               Baseind = Baseind+1;
%               BaseSpikeTimes(Baseind) = t;
% 
%               if baselineLIFGating
%                 basegateCount = basegateDur/dt;
%               end
%             else
%               basegateCount = basegateCount-1;
%             end
% 
%             if basegateCount>0
%               basegate = 1;
%             else
%               basegate = 0;
%             end
          end
        end
      else
        spikesum = loadedSpikes(tind,1:end-1);
        spikesumb = loadedSpikes(tind,end);
      end

      if spikesumb
        if mod(Baseind,1000)==0
          BaseSpikeTimes(length(BaseSpikeTimes)+1000) = 0;
        end
        Baseind = Baseind+1;
        BaseSpikeTimes(Baseind) = t;

        if baselineLIFGating
          basegateCount = basegateDur/dt;
        end
      else
        basegateCount = basegateCount-1;
      end

      if basegateCount>0
        basegate = 1;
      else
        basegate = 0;
      end

      
      if ~isempty(saveVCOSpikesFileName)
        % each line: number of spikes from 1:nPostIn on this
        % time step for each VCO and the baseline
        % we'll buffer the output to speed things up
        buflin = [];
        if runNetworkVCOs
          for vcoi=1:nVCOs
            buflin = [buflin sum(spikes(1:nPostIn,vcoi))];
          end
        end
        if runNetworkBaseline
          buflin = [buflin sum(spikesb(1:nPostIn))];
        end
        buffer = [buffer; buflin];
        if size(buffer,1)==50
          fid = fopen(saveVCOSpikesFileName,'a');
          for br=1:size(buffer,1)
            for bc=1:size(buffer,2)
              fprintf(fid,'%d ',buffer(br,bc));
            end
            fprintf(fid,'\n');
          end
          fclose(fid);
          buffer = [];
        end
      end

      state(:,tind) = [v(1,:)'; vb(1)];

      if cheatHDRectification
        spikesum = spikesum.*(cos(prefHDs-curHD)>0);
      end
      
      % post cell 1 = non-gated lif
      if stablePost
        post(1,tind) = post(1,tind-1)*membraneDecay + postWeight*sum(spikesum)+baselineMult*baseWeight*spikesumb;
        if post(1,tind)>1
          post(1,tind) = -1e-10;
        end
      else
        GABAtau = 15; % ms
        % post cell 1 = non-gated simple model
        postinhib(1,tind) = postinhib(1,tind-1)*exp(-dt/GABAtau) + postWeight*sum(spikesum)+baselineMult*baseWeight*spikesumb;
        oldv = post(1,tind-1);
        oldu = postu(1,tind-1);
        post(1,tind) = oldv + dt*(k*(oldv-vr).*(oldv-vt) - oldu + baseI + postinhib(1,tind))/Cf;
        postu(1,tind) = oldu + dt*a*(b*(oldv-vr)-oldu);
        if post(1,tind)>vpeak
          postu(1,tind) = postu(1,tind) + d;
          post(1,tind-1) = vpeak;
          post(1,tind) = c;
        end
      end
      % post cell 2 = gated lif
      post(2,tind) = post(2,tind-1)*gatedmembraneDecay + basegate*gatedpostWeight*sum(spikesum);
      if post(2,tind)>1
        post(2,tind) = -1e-10;
      end
      % post cell 3 = non-gated resonant
      post(3,tind) = post(3,tind-1) + dt*((postc+postw*i)*post(3,tind-1) + respostWeight*sum(spikesum)+baselineMult*resbaseWeight*spikesumb);
      postspikes(tind) = real(post(3,tind))>1;
      if real(post(3,tind))>1
        post(3,tind) = i;
      end

      vpostb(tind) =  vpostb(tind-1)*membraneDecay + baseWeight*spikesumb;
      for vind=1:nVCOs
        vpost(vind,tind) = vpost(vind,tind-1)*membraneDecay + basegate*postWeight*spikesum(vind);
      end

      if runningPlot
        figure(runh);
        plot(pos(1,:),pos(2,:),'.',pos(1,find(postspikes>0)),pos(2,find(postspikes>0)),'R.'), title('network VCO spatial firing'); xlabel('Position (m)','fontsize',12), ylabel('Position (m)','fontsize',12);
        set(gca,'xlim',[-1 1],'ylim',[-1 1]);
        drawnow
      end
    end
    toc
    %% simulation is complete, now for some post-processing and plotting:

    if ~isempty(loadVCOSpikesFileName)
      for vcoi=1:nVCOs
        VCOSpikeTimes{vcoi} = dt*find(loadedSpikes(:,vcoi)>0);
      end
      BaseSpikeTimes = dt*find(loadedSpikes(:,end)>0);
    end
    
    % flush spike times buffer
    if ~isempty(saveVCOSpikesFileName)
      if size(buffer,1)>0
        fid = fopen(saveVCOSpikesFileName,'a');
        for br=1:size(buffer,1)
          for bc=1:size(buffer,2)
            fprintf(fid,'%d ',buffer(br,bc));
          end
        end
        fprintf(fid,'\n');
        fclose(fid);
        buffer = [];
      end
    end

    clear C

    if runNetwork
      for ce=1:npost
        spikeTimes = find(post(ce,:)==-1e-10);

        nspikes(ce) = length(spikeTimes);

        ISIs = dt*diff(spikeTimes)/1000;

        % count ISIs as occuring between initial
        % spikes of bursts (and additional spikes as occuring in a given
        % burst if they are < 50 ms apart);
        fISIs = [];
        csum = 0;
        for is=1:length(ISIs)
          csum = csum + ISIs(is);
          if ISIs(is)>.050
            fISIs = [fISIs csum];
            csum=0;
          end;
        end
        ISIs = fISIs;

        ISImeans(ce) = mean(ISIs);
        ISIstds(ce) = std(ISIs);
        ISIstabilities(ce) = 5*ISImeans(ce)^3/16/pi^2/(eps+ISIstds(ce))^2;
      end
    end
  end

  for vco=1:nVCOs
    VCOSpikeTimes{vco}(VCOSpikeTimes{vco}==0) = [];
    if length(VCOSpikeTimes{vco})>40000
      VCOSpikeTimes{vco} = VCOSpikeTimes{vco}(1:10:end);
    end
  end
  BaseSpikeTimes(BaseSpikeTimes==0) = [];
  if length(BaseSpikeTimes)>40000
    BaseSpikeTimes = BaseSpikeTimes(1:10:end);
  end

  % transform the simple model output so that our LIF spike-detection in
  % the figures works
  if stablePost==0
    post(post==vpeak) = -1e-10;
  end

  % tweak straight line trajectories to avoid error below
  if min(pos(1,:))==max(pos(1,:))
    pos(1,end) = pos(1,end)+eps;
  end
  if min(pos(2,:))==max(pos(2,:))
    pos(2,end) = pos(2,end)+eps;
  end
  

  %% FIGURES from here to the end of the file
  %% (and they are a bit of a mess!)
  % let's set up a few aspects of the appearance:
  set(0,'defaulttextfontsize',12);
  if journalChargesALotForColor
    % for images, white is largest value
    figure; set(0,'defaultfigurecolormap',flipud(gray)); close
    % we only ever plot <=3 colors at one time in the main figures
    set(0,'defaultaxescolororder',[0 0 1; 0 0.5 0; 1 0 0;]);
    %     set(0,'defaultaxescolororder',[0 0 0; .75 .75 .75]);
    %     set(0,'defaultaxeslinestyleorder',{'-','--'});
    %     % this means line 1 is solid/black, line 2 is solid/gray, line 3 dashed/black, line 4 dashed/gray
  end

  if plotVCOsAndPosts
    % Note: this won't plot if using loaded spikes, have to save this when
    % the VCO spikes are first generated!
    figure; plot(state(:,1:40000)'); set(gcf,'position',[82 404 3*257 125])
    xlabel('Time (ms)','fontsize',12); ylabel('Voltage (mV)','fontsize',12)
    set(gca,'ylim',[-60 50]);
    set(gca,'fontsize',12)
    set(gca,'box','off');
    if saveFigures && exist(['fig_RS' typePrefix '_' simprefix '_2D_vcos_trace_0a.eps'])~=2
      print_eps(['fig_RS' typePrefix '_' simprefix '_2D_vcos_trace_0a.eps'])
      saveas(gcf,['fig_RS' typePrefix '_' simprefix '_2D_vcos_trace_0a.fig'])
    end

    plotduration = 4000; % ms
    % plot non-gated lif post:
    figure; plot(dt:dt:plotduration,post(1,1:(plotduration/dt)));
    if any(post(1,1:(plotduration/dt))==-1e-10)
      hold on; plot(dt*find(post(1,1:(plotduration/dt))==-1e-10),1,'r.','MarkerSize',16);
    end
    set(gcf,'position',[82 404 3*257 125])
    xlabel('Time (ms)','fontsize',12); ylabel('Activity','fontsize',12)
    if stablePost
      set(gca,'ylim',[-0.1 1.1]);
    end
    set(gca,'fontsize',12)
    set(gca,'box','off');
    if saveFigures && exist(['fig_RS' typePrefix '_' simprefix '_2D_post_nogatelif_trace_0a.eps'])~=2
      print_eps(['fig_RS' typePrefix '_' simprefix '_2D_post_nogatelif_trace_0a.eps'])
      saveas(gcf,['fig_RS' typePrefix '_' simprefix '_2D_post_nogatelif_trace_0a.fig'])
    end

    % plot gated lif post:
    if nVCOs==1
      figure; plot(dt:dt:plotduration,post(2,1:(plotduration/dt))==-1e-10);
    else
%       figure; plot(dt:dt:plotduration,post(2,1:(plotduration/dt))==-1e-10);
      figure; plot(dt:dt:plotduration,post(2,1:(plotduration/dt)));      
    end
    if any(post(2,1:(plotduration/dt))==-1e-10)
      hold on; plot(dt*find(post(2,1:(plotduration/dt))==-1e-10),1,'r.','MarkerSize',16);
    end
    set(gcf,'position',[82 404 3*257 125])
    xlabel('Time (ms)','fontsize',12); ylabel('Activity','fontsize',12)
    set(gca,'ylim',[-0.1 1.1]);
    set(gca,'box','off');
    set(gca,'fontsize',12)
    if saveFigures && exist(['fig_RS' typePrefix '_' simprefix '_2D_post_gating_trace_0a.eps'])~=2
      print_eps(['fig_RS' typePrefix '_' simprefix '_2D_post_gating_trace_0a.eps'])
      saveas(gcf,['fig_RS' typePrefix '_' simprefix '_2D_post_gating_trace_0a.fig'])
    end

    % plot non-gated res post:
    figure; plot(dt:dt:plotduration,real(post(3,1:(plotduration/dt))));
    if any(postspikes)
      hold on; plot(dt*find(postspikes(1:(plotduration/dt))>0),1,'r.','MarkerSize',16);
    end
    set(gcf,'position',[82 404 3*257 125])
    xlabel('Time (ms)','fontsize',12); ylabel('Activity','fontsize',12)
    set(gca,'ylim',[-1.1 1.1]);
    set(gca,'box','off');
    set(gca,'fontsize',12)
    if saveFigures && exist(['fig_RS' typePrefix '_' simprefix '_2D_post_nogateres_trace_0a.eps'])~=2
      print_eps(['fig_RS' typePrefix '_' simprefix '_2D_post_nogateres_trace_0a.eps'])
      saveas(gcf,['fig_RS' typePrefix '_' simprefix '_2D_post_nogateres_trace_0a.fig'])
    end
  end

  if plotType==0
    % these in normalized coords:
    trajH = 3/4*1/2;
    trajW = 3/4*1/3;
    autoH = trajH;
    autoW = trajW;
    diffH = 2/3*1/3;
    diffW = trajW;

    trajL = 1/4*1/3;
    trajB = 1/6*1/2;
    autoL = 1/6*1/3;
    autoB = trajB;
    diffL = 1/5*1/3;
    diffB1 = .007+1/8*1/2;
    diffB2 = diffB1+0.005;
    diffB3 = trajB;

    figure('position',[82 404 647 484])
    if runNetwork
      posskip = 500;
      thr = 1.5;
      if stablePost
        subplot('position',[trajL+0 trajB+1/2 trajW trajH]); plot(pos(1,1:posskip:length(pos)),pos(2,1:posskip:length(pos)),'.','MarkerSize',8); hold on; plot(pos(1,find(post(2,:)==-1e-10)),pos(2,find(post(2,:)==-1e-10)),'R.','MarkerSize',15)
      else
        subplot('position',[trajL+0 trajB+1/2 trajW trajH]); plot(pos(1,1:posskip:length(pos)),pos(2,1:posskip:length(pos)),'.','MarkerSize',8); hold on; plot(pos(1,find(post(1,:)==-1e-10)),pos(2,find(post(1,:)==-1e-10)),'R.','MarkerSize',15)
      end
      if ~showAbstractFig
        title(['Spatial firing - 0 to ' num2str(simdur/1000) ' s '],'fontsize',12)
      else
        title('Network model spatial firing','fontsize',12)
      end
      xlabel('Position (m)','fontsize',12), ylabel('Position (m)','fontsize',12)
      axis equal
      set(gca,'xlim',[min(pos(1,:)) max(pos(1,:))],'ylim',[min(pos(2,:)) max(pos(2,:))])
      set(gca,'box','off');
      set(gca,'fontsize',12)

      edges{1} = linspace(min(pos(1,:)),max(pos(1,:)),50);
      edges{2} = linspace(min(pos(2,:)),max(pos(2,:)),50);
      if stablePost
        rate = hist3([pos(1,find(post(2,:)==-1e-10)); pos(2,find(post(2,:)==-1e-10))]','Edges',edges);
      else
        rate = hist3([pos(1,find(post(1,:)==-1e-10)); pos(2,find(post(1,:)==-1e-10))]','Edges',edges);
      end
      subplot('position',[autoL+1/3 autoB+1/2 autoW autoH]); imagesc(rot90(conv2(rate,rate,'same')))
      if showAbstractFig
        title('Network model autocorrelogram','fontsize',12)
      elseif stablePost
        title(['Spatial autocorrelogram - 0 to ' num2str(simdur/1000) ' s '],'fontsize',12)
      end
      set(gca,'xtick',[],'ytick',[]);
      set(gca,'fontsize',12)
      axis tight

      if ~showAbstractFig
        subplot('position',[trajL+0 trajB+0 trajW trajH]); plot(pos(1,1:posskip:length(pos)/4),pos(2,1:posskip:length(pos)/4),'.','MarkerSize',8);
        hold on;
        plot(pos(1,find(post(2,1:length(pos)/4)==-1e-10)),pos(2,find(post(2,1:length(pos)/4)==-1e-10)),'R.','MarkerSize',15)
        title('Spatial firing - 0 to 40 s ','fontsize',12)
        title(['Spatial firing - 0 to ' num2str(simdur/1000/4) ' s '],'fontsize',12)
        xlabel('Position (m)','fontsize',12), ylabel('Position (m)','fontsize',12)
        set(gca,'fontsize',12)
        axis equal
        set(gca,'xlim',[min(pos(1,:)) max(pos(1,:))],'ylim',[min(pos(2,:)) max(pos(2,:))])
        set(gca,'box','off');
        rate = hist3([pos(1,find(post(2,1:length(pos)/4)==-1e-10)); pos(2,find(post(2,1:length(pos)/4)==-1e-10))]','Edges',edges);
        subplot('position',[autoL+1/3 autoB+0 autoW autoH]); imagesc(rot90(conv2(rate,rate,'same')))
        title(['Spatial autocorrelogram - 0 to ' num2str(simdur/1000/4) ' s '],'fontsize',12)
        set(gca,'xtick',[],'ytick',[]);
        set(gca,'fontsize',12)
        axis tight
      end


      if ncells>1 && useGapJunctions==0 && useNoise>0
        % to keep filesizes lower when many noisy cells are synaptically
        % coupled:
        spikeskips = ncells/11;
      else
        spikeskips = 1;
      end
      if runNetworkVCOs
        if nVCOs>1
          subplot('position',[diffL+2/3 diffB1+2/3 diffW diffH]); plot(mod(vcohist(1,round(VCOSpikeTimes{1}(1:spikeskips:end)/dt)),2*pi),'.')
        else
          subplot('position',[diffL+2/3 trajB+1/2 diffW trajH]); plot(mod(vcohist(1,round(VCOSpikeTimes{1}(1:spikeskips:end)/dt)),2*pi),'.')
        end
        if plotVelocity
          hold on; plot(5-1e4*fdxs(round(VCOSpikeTimes{1}(1:spikeskips:end)/dt)),'k'); set(gca,'ylim',[0 6.2]);
        end
        title('VCO 1 phase drift','fontsize',12)
        ylabel({'Phase difference','(rad)'},'fontsize',12)
        set(gca,'xlim',[0 length(round(VCOSpikeTimes{1}(1:spikeskips:end)))],'ylim',[0 2*pi],'xtick',[])
        set(gca,'box','off');
        set(gca,'fontsize',12)
        if nVCOs>1
          subplot('position',[diffL+2/3 diffB2+1/3 diffW diffH]); plot(mod(vcohist(2,round(VCOSpikeTimes{2}(1:spikeskips:end)/dt)),2*pi),'.')
          title('VCO 2 phase drift','fontsize',12)
          ylabel({'Phase difference','(rad)'},'fontsize',12)
          set(gca,'xlim',[0 length(round(VCOSpikeTimes{2}(1:spikeskips:end)))],'ylim',[0 2*pi],'xtick',[])
          set(gca,'box','off');
          set(gca,'fontsize',12)
        end
      end

      if runNetworkBaseline
        if nVCOs>1
          subplot('position',[diffL+2/3 diffB3+0 diffW diffH]); plot(mod(basehist(round(BaseSpikeTimes(1:spikeskips:end)/dt)),2*pi),'.')
        else
          subplot('position',[diffL+2/3 diffB3+0 diffW trajH]); plot(mod(basehist(round(BaseSpikeTimes(1:spikeskips:end)/dt)),2*pi),'.')
        end
        title('Baseline phase drift','fontsize',12)
        xlabel('Time','fontsize',12)
        ylabel({'Phase difference','(rad)'},'fontsize',12)
        set(gca,'xlim',[0 length(round(BaseSpikeTimes(1:spikeskips:end)))],'ylim',[0 2*pi],'xtick',[])
        set(gca,'box','off');
        set(gca,'fontsize',12)
      end
    end
    if runAbstract && showAbstractFig
      spikeskip = 100;
      subplot('position',[trajL+0 trajB+0 trajW trajH]); plot(pos(1,1:posskip:length(pos)),pos(2,1:posskip:length(pos)),'.','MarkerSize',8);
      hold on; plot(pos(1,spikeskip*find(abstractGrid(1:spikeskip:end)>abstractThr)),pos(2,spikeskip*find(abstractGrid(1:spikeskip:end)>abstractThr)),'R.','MarkerSize',15)
      title('Abstract model spatial firing','fontsize',12)
      xlabel('Position (m)','fontsize',12), ylabel('Position (m)','fontsize',12)
      axis equal
      set(gca,'xlim',[min(pos(1,:)) max(pos(1,:))],'ylim',[min(pos(2,:)) max(pos(2,:))])
      set(gca,'box','off');
      set(gca,'fontsize',12)

      rate = hist3([pos(1,find(abstractGrid>abstractThr)); pos(2,find(abstractGrid>abstractThr))]','Edges',edges);
      %   rate = hist3([pos(1,find(abstractGrid>abstractThr)); pos(2,find(abstractGrid>abstractThr))]',[50 50]);
      subplot('position',[autoL+1/3 autoB+0 autoW autoH]); imagesc(rot90(conv2(rate,rate,'same')))
      title('Abstract model autocorrelogram','fontsize',12)
      set(gca,'xtick',[],'ytick',[]);
      set(gca,'fontsize',12)
      axis tight
    end
    if saveFigures && exist(['fig_RS' typePrefix '_' simprefix '_2D_all_0a.eps'])~=2
      save(['fig_RS' typePrefix '_' simprefix '_2D_variables_0a.mat'])
      % free some memory:
      %   clear pos abstractGrid vhist vcohist basehist post
      print_eps(['fig_RS' typePrefix '_' simprefix '_2D_all_0a.eps'])
      saveas(gcf,['fig_RS' typePrefix '_' simprefix '_2D_all_0a.fig'])
      close
    end
  elseif plotType==1
    % these in normalized coords:
    trajH = 3/4*1/2;
    trajW = 3/4*1/3;
    diffH = 2/3*1/3;
    diffW = trajW;

    col1L = 1/4*1/3;
    col2L = 1/3+1/5*1/3;
    col3L = 2/3+1/5*1/3;
    row2B = 1/6*1/2;
    diffL = 1/5*1/3;
    diffB1 = .007+1/8*1/2;
    diffB2 = diffB1+0.005;

    if ncells>1 && useNoise>0 % && useGapJunctions==0
      % to keep filesizes lower when many noisy cells are synaptically
      % coupled:
      spikeskips = ncells/23;
    else
      spikeskips = 1;
    end

    mainplot = 2;
    if mainplot==1
      if stablePost
        plottype = 'Integrator postsynaptic cell';
      else
        plottype = 'Inhibitory oscillators';
      end
    else
      plottype = 'Gated postsynaptic cell';
    end
    figure('position',[82 404 647 484])
    posskip = 500;
    if nVCOs==2
      % gated lif:
      subplot('position',[col1L+0 row2B+1/2 trajW trajH]);
      plot(pos(1,1:posskip:length(pos)),pos(2,1:posskip:length(pos)),'.','MarkerSize',8);
      hold on; plot(pos(1,find(post(mainplot,:)==-1e-10)),pos(2,find(post(mainplot,:)==-1e-10)),'R.','MarkerSize',15)
      title(plottype,'fontsize',12)
      xlabel('Position (m)','fontsize',12), ylabel('Position (m)','fontsize',12)
      set(gca,'xlim',[min(pos(1,:)) max(pos(1,:))],'ylim',[min(pos(2,:)) max(pos(2,:))])
      set(gca,'box','off');
      set(gca,'fontsize',12)

      % gated lif xcorr:
      edges{1} = linspace(min(pos(1,:)),max(pos(1,:)),50);
      edges{2} = linspace(min(pos(2,:)),max(pos(2,:)),50);
      rate = hist3([pos(1,find(post(mainplot,:)==-1e-10)); pos(2,find(post(mainplot,:)==-1e-10))]','Edges',edges);
      subplot('position',[col2L+0 row2B+1/2 trajW trajH]); imagesc(rot90(conv2(rate,rate,'same')))
      if showAbstractFig
        title('Network model autocorrelogram','fontsize',12)
      else
        title([plottype ' autocorrelogram'],'fontsize',12)
      end
      set(gca,'xtick',[],'ytick',[]);
      set(gca,'fontsize',12)
      axis tight

      if showAbstractFig
        % abstract:
        abstractThr = 3;
        subplot('position',[col1L+0 row2B+0 trajW trajH]);
        plot(pos(1,1:posskip:length(pos)),pos(2,1:posskip:length(pos)),'.','MarkerSize',8);
        hold on; plot(pos(1,find(abstractGrid>abstractThr)),pos(2,find(abstractGrid>abstractThr)),'R.','MarkerSize',15)
        title('Abstract model','fontsize',12)
        xlabel('Position (m)','fontsize',12), ylabel('Position (m)','fontsize',12)
        set(gca,'xlim',[min(pos(1,:)) max(pos(1,:))],'ylim',[min(pos(2,:)) max(pos(2,:))])
        set(gca,'box','off');
        set(gca,'fontsize',12)

        rate = hist3([pos(1,find(abstractGrid>abstractThr)); pos(2,find(abstractGrid>abstractThr))]','Edges',edges);
        subplot('position',[col2L+0 row2B trajW trajH]); imagesc(rot90(conv2(rate,rate,'same')))
        title(['Abstract model autocorrelogram'],'fontsize',12)
        set(gca,'xtick',[],'ytick',[]);
        set(gca,'fontsize',12)
        axis tight
      else
        % non-gated res:
        subplot('position',[col1L+0 row2B+0 trajW trajH]);
        plot(pos(1,1:posskip:length(pos)),pos(2,1:posskip:length(pos)),'.','MarkerSize',8);
        hold on; plot(pos(1,find(postspikes>0)),pos(2,find(postspikes>0)),'R.','MarkerSize',15)
        title('Resonator postsynaptic cell','fontsize',12)
        xlabel('Position (m)','fontsize',12), ylabel('Position (m)','fontsize',12)
        set(gca,'xlim',[min(pos(1,:)) max(pos(1,:))],'ylim',[min(pos(2,:)) max(pos(2,:))])
        set(gca,'box','off');
        set(gca,'fontsize',12)

        rate = hist3([pos(1,find(postspikes==1)); pos(2,find(postspikes==1))]','Edges',edges);
        subplot('position',[col2L+0 row2B trajW trajH]); imagesc(rot90(conv2(rate,rate,'same')))
        title(['Resonator autocorrelogram'],'fontsize',12)
        set(gca,'xtick',[],'ytick',[]);
        set(gca,'fontsize',12)
        axis tight
      end

      % vco 1 phase diff:
      subplot('position',[col3L diffB2+2/3 diffW diffH]);
      plot(mod(vcohist(1,round(VCOSpikeTimes{1}(1:spikeskips:end)/dt)),2*pi),'.')
      title('VCO 1 phase drift','fontsize',12)
      ylabel({'Phase difference','(rad)'},'fontsize',12)
      set(gca,'xlim',[0 length(round(VCOSpikeTimes{1}(1:spikeskips:end)))],'ylim',[0 2*pi],'xtick',[])
      set(gca,'box','off');
      set(gca,'fontsize',12)

      subplot('position',[diffL+2/3 diffB2+1/3 diffW diffH]);
      plot(mod(vcohist(2,round(VCOSpikeTimes{2}(1:spikeskips:end)/dt)),2*pi),'.')
      title('VCO 2 phase drift','fontsize',12)
      ylabel({'Phase difference','(rad)'},'fontsize',12)
      set(gca,'xlim',[0 length(round(VCOSpikeTimes{2}(1:spikeskips:end)))],'ylim',[0 2*pi],'xtick',[])
      set(gca,'box','off');
      set(gca,'fontsize',12)

      % base phase diff:
      subplot('position',[col3L diffB2 diffW diffH]);
      plot(mod(basehist(round(BaseSpikeTimes(1:spikeskips:end)/dt)),2*pi),'.')
      title('Baseline phase drift','fontsize',12)
      xlabel('Time','fontsize',12)
      ylabel({'Phase difference','(rad)'},'fontsize',12)
      set(gca,'xlim',[0 length(round(BaseSpikeTimes(1:spikeskips:end)))],'ylim',[0 2*pi],'xtick',[])
      set(gca,'box','off');
      set(gca,'fontsize',12)
    end
    if saveFigures && exist(['fig_RS' typePrefix '_' simprefix '_2D_all_0a.eps'])~=2
      save(['fig_RS' typePrefix '_' simprefix '_2D_variables_0a.mat'])
      % free some memory:
      %   clear pos abstractGrid vhist vcohist basehist post
      print_eps(['fig_RS' typePrefix '_' simprefix '_2D_all_0a.eps'])
      saveas(gcf,['fig_RS' typePrefix '_' simprefix '_2D_all_0a.fig'])
      close
    end
  elseif plotType==2
    % these in normalized coords:
    trajH = 3/4*1/2; % <1/2
    trajW = 3/4*1/3; % <1/3

    col1L = 1/4*1/3; % <1/3-trajW
    col2L = 1/3+1/5*1/3;
    col3L = 2/3+1/5*1/3;
    row2B = 1/6*1/2; % <1/3-trajH

    figure('position',[82 404 647 484])
    posskip = 500;
    spikeskip = 1;
    % non-gated lif:
    subplot('position',[col1L+0 row2B+1/2 trajW trajH]);
    plot(pos(1,1:posskip:length(pos)),pos(2,1:posskip:length(pos)),'.','MarkerSize',8);
    hold on; plot(pos(1,spikeskip*find(post(1,1:spikeskip:end)==-1e-10)),pos(2,spikeskip*find(post(1,1:spikeskip:end)==-1e-10)),'R.','MarkerSize',15)
    if stablePost
      title('Integrator postsynaptic cell','fontsize',12)
    else
      title('Inhibitory oscillators','fontsize',12)
    end
    xlabel('Position (m)','fontsize',12), ylabel('Position (m)','fontsize',12)
    set(gca,'xlim',[min(pos(1,:)) max(pos(1,:))],'ylim',[min(pos(2,:)) max(pos(2,:))])
    set(gca,'box','off');
    set(gca,'fontsize',12)

    % gated lif:
    subplot('position',[col2L+0 row2B+1/2 trajW trajH]);
    plot(pos(1,1:posskip:length(pos)),pos(2,1:posskip:length(pos)),'.','MarkerSize',8);
    hold on; plot(pos(1,spikeskip*find(post(2,1:spikeskip:end)==-1e-10)),pos(2,spikeskip*find(post(2,1:spikeskip:end)==-1e-10)),'R.','MarkerSize',15)
    title('Gated postsynaptic cell','fontsize',12)
    xlabel('Position (m)','fontsize',12), ylabel('Position (m)','fontsize',12)
    set(gca,'xlim',[min(pos(1,:)) max(pos(1,:))],'ylim',[min(pos(2,:)) max(pos(2,:))])
    set(gca,'box','off');
    set(gca,'fontsize',12)

    % non-gated res:
    subplot('position',[col3L+0 row2B+1/2 trajW trajH]);
    plot(pos(1,1:posskip:length(pos)),pos(2,1:posskip:length(pos)),'.','MarkerSize',8);
    hold on; plot(pos(1,spikeskip*find(postspikes(1:spikeskip:end)>0)),pos(2,spikeskip*find(postspikes(1:spikeskip:end)>0)),'R.','MarkerSize',15)
    title('Resonator postsynaptic cell','fontsize',12)
    xlabel('Position (m)','fontsize',12), ylabel('Position (m)','fontsize',12)
    set(gca,'xlim',[min(pos(1,:)) max(pos(1,:))],'ylim',[min(pos(2,:)) max(pos(2,:))])
    set(gca,'box','off');
    set(gca,'fontsize',12)

    % abstract:
    if nVCOs==1
      abstractThr = 1.8;
    else
      abstractThr = 3;
    end
    spikeskip = 100;
    subplot('position',[col1L+0 row2B+0 trajW trajH]);
    plot(pos(1,1:posskip:length(pos)),pos(2,1:posskip:length(pos)),'.','MarkerSize',8);
    hold on; plot(pos(1,spikeskip*find(abstractGrid(1:spikeskip:end)>abstractThr)),pos(2,spikeskip*find(abstractGrid(1:spikeskip:end)>abstractThr)),'R.','MarkerSize',15)
    title('Abstract model','fontsize',12)
    xlabel('Position (m)','fontsize',12), ylabel('Position (m)','fontsize',12)
    set(gca,'xlim',[min(pos(1,:)) max(pos(1,:))],'ylim',[min(pos(2,:)) max(pos(2,:))])
    set(gca,'box','off');
    set(gca,'fontsize',12)

    % vco 1 phase diff:
    spikeskips = 1;
    subplot('position',[col2L+0 row2B+0 trajW trajH]);
    plot(1:spikeskips:length(VCOSpikeTimes{1}),mod(vcohist(1,round(VCOSpikeTimes{1}(1:spikeskips:end)/dt)),2*pi),'.')
    title('VCO 1 phase drift','fontsize',12)
    xlabel('Time','fontsize',12)
    ylabel({'Phase difference','(rad)'},'fontsize',12)
    set(gca,'xlim',[0 length(round(VCOSpikeTimes{1}(1:spikeskips:end)))],'ylim',[0 2*pi],'xtick',[])
    set(gca,'box','off');
    set(gca,'fontsize',12)

    % base phase diff:
    subplot('position',[col3L row2B+0 trajW trajH]);
    plot(1:spikeskips:length(BaseSpikeTimes),mod(basehist(round(BaseSpikeTimes(1:spikeskips:end)/dt)),2*pi),'.')
    title('Baseline phase drift','fontsize',12)
    xlabel('Time','fontsize',12)
    ylabel({'Phase difference','(rad)'},'fontsize',12)
    set(gca,'xlim',[0 length(round(BaseSpikeTimes(1:spikeskips:end)))],'ylim',[0 2*pi],'xtick',[])
    set(gca,'box','off');
    set(gca,'fontsize',12)

    if saveFigures && exist(['fig_RS' typePrefix '_' simprefix '_2D_all_0a.eps'])~=2
      save(['fig_RS' typePrefix '_' simprefix '_2D_variables_0a.mat'])
      % if print_eps complains, you can free some memory: (or even clear all)
      %   clear pos abstractGrid vhist vcohist basehist post
      print_eps(['fig_RS' typePrefix '_' simprefix '_2D_all_0a.eps'])
      saveas(gcf,['fig_RS' typePrefix '_' simprefix '_2D_all_0a.fig'])
      close
    end
  elseif plotType==3
    % these in normalized coords:
    trajH = 3/4*1/3; % <1/2
    trajW = 3/4*1/3; % <1/3

    col1L = 1/4*1/3; % <1/3-trajW
    col2L = 1/3+1/5*1/3;
    col3L = 2/3+1/5*1/3;
    row2B = 1/6*1/3; % <1/3-trajH

    figure('position',[82 404 647 484])
    posskips = 500;
    %     if nVCOs==1
    % non-gated lif:
    subplot('position',[col1L+0 row2B+2/3 trajW trajH]);
    plot(pos(1,1:posskips:length(pos)),pos(2,1:posskips:length(pos)),'.','MarkerSize',8);
    hold on; plot(pos(1,find(post(1,:)==-1e-10)),pos(2,find(post(1,:)==-1e-10)),'R.','MarkerSize',15)
    if stablePost
      title('Integrator postsynaptic cell','fontsize',12)
    else
      title('Inhibitory oscillators','fontsize',12)
    end
    xlabel('Position (m)','fontsize',12), ylabel('Position (m)','fontsize',12)
    set(gca,'xlim',[min(pos(1,:)) max(pos(1,:))],'ylim',[min(pos(2,:)) max(pos(2,:))])
    set(gca,'box','off');
    set(gca,'fontsize',12)

    edges{1} = linspace(min(pos(1,:)),max(pos(1,:)),50);
    edges{2} = linspace(min(pos(2,:)),max(pos(2,:)),50);
    rate = hist3([pos(1,find(post(1,:)==-1e-10)); pos(2,find(post(1,:)==-1e-10))]','Edges',edges);
    subplot('position',[col1L+0 row2B+1/3 trajW trajH]); imagesc(rot90(conv2(rate,rate,'same')))
    title('Integrator autocorrelogram','fontsize',12)
    set(gca,'xtick',[],'ytick',[]);
    axis tight
    set(gca,'fontsize',12)

    % gated lif:
    subplot('position',[col2L+0 row2B+2/3 trajW trajH]);
    plot(pos(1,1:posskips:length(pos)),pos(2,1:posskips:length(pos)),'.','MarkerSize',8);
    hold on; plot(pos(1,find(post(2,:)==-1e-10)),pos(2,find(post(2,:)==-1e-10)),'R.','MarkerSize',15)
    title('Gated, integrator postsynaptic cell','fontsize',12)
    xlabel('Position (m)','fontsize',12), ylabel('Position (m)','fontsize',12)
    set(gca,'xlim',[min(pos(1,:)) max(pos(1,:))],'ylim',[min(pos(2,:)) max(pos(2,:))])
    set(gca,'box','off');
    set(gca,'fontsize',12)

    rate = hist3([pos(1,find(post(2,:)==-1e-10)); pos(2,find(post(2,:)==-1e-10))]','Edges',edges);
    subplot('position',[col2L+0 row2B+1/3 trajW trajH]); imagesc(rot90(conv2(rate,rate,'same')))
    title('Gated, integrator autocorrelogram','fontsize',12)
    set(gca,'xtick',[],'ytick',[]);
    axis tight
    set(gca,'fontsize',12)

    % non-gated res:
    subplot('position',[col3L+0 row2B+2/3 trajW trajH]);
    plot(pos(1,1:posskips:length(pos)),pos(2,1:posskips:length(pos)),'.','MarkerSize',8);
    hold on; plot(pos(1,find(postspikes>0)),pos(2,find(postspikes>0)),'R.','MarkerSize',15)
    title('Resonator postsynaptic cell','fontsize',12)
    xlabel('Position (m)','fontsize',12), ylabel('Position (m)','fontsize',12)
    set(gca,'xlim',[min(pos(1,:)) max(pos(1,:))],'ylim',[min(pos(2,:)) max(pos(2,:))])
    set(gca,'box','off');
    set(gca,'fontsize',12)

    rate = hist3([pos(1,find(postspikes>0)); pos(2,find(postspikes>0))]','Edges',edges);
    subplot('position',[col3L+0 row2B+1/3 trajW trajH]); imagesc(rot90(conv2(rate,rate,'same')))
    title('Resonator autocorrelogram','fontsize',12)
    set(gca,'xtick',[],'ytick',[]);
    axis tight
    set(gca,'fontsize',12)

    % abstract:
    if nVCOs==1
      abstractThr = 1.8;
    elseif nVCOs==2
      abstractThr = 3;
    else
      abstractThr = nVCOs*1.5;
    end
    subplot('position',[col1L+0 row2B+0 trajW trajH]);
    plot(pos(1,1:posskips:length(pos)),pos(2,1:posskips:length(pos)),'.','MarkerSize',8);
    hold on; plot(pos(1,find(abstractGrid>abstractThr)),pos(2,find(abstractGrid>abstractThr)),'R.','MarkerSize',15)
    title('Abstract model','fontsize',12)
    xlabel('Position (m)','fontsize',12), ylabel('Position (m)','fontsize',12)
    set(gca,'xlim',[min(pos(1,:)) max(pos(1,:))],'ylim',[min(pos(2,:)) max(pos(2,:))])
    set(gca,'box','off');
    set(gca,'fontsize',12)

    % vco 1 phase diff:
    if ncells<5000
      spikeskips = 1;
    else
      spikeskips = 25;
    end
    subplot('position',[col2L+0 row2B+0 trajW trajH]);
    plot(mod(vcohist(1,round(VCOSpikeTimes{1}(1:spikeskips:end)/dt)),2*pi),'.')
    title('VCO 1 phase drift','fontsize',12)
    xlabel('Time','fontsize',12)
    ylabel({'Phase difference','(rad)'},'fontsize',12)
    set(gca,'xlim',[0 length(round(VCOSpikeTimes{1}(1:spikeskips:end)))],'ylim',[0 2*pi],'xtick',[])
    set(gca,'box','off');
    set(gca,'fontsize',12)

    % base phase diff:
    subplot('position',[col3L row2B+0 trajW trajH]);
    plot(mod(basehist(round(BaseSpikeTimes(1:spikeskips:end)/dt)),2*pi),'.')
    title('Baseline phase drift','fontsize',12)
    xlabel('Time','fontsize',12)
    ylabel({'Phase difference','(rad)'},'fontsize',12)
    set(gca,'xlim',[0 length(round(BaseSpikeTimes(1:spikeskips:end)))],'ylim',[0 2*pi],'xtick',[])
    set(gca,'box','off');
    set(gca,'fontsize',12)

    if saveFigures && exist(['fig_RS' typePrefix '_' simprefix '_2D_all_0a.eps'])~=2
      save(['fig_RS' typePrefix '_' simprefix '_2D_variables_0a.mat'])
      % if print_eps complains, you can free some memory: (or even clear all)
      %   clear pos abstractGrid vhist vcohist basehist post
      print_eps(['fig_RS' typePrefix '_' simprefix '_2D_all_0a.eps'])
      saveas(gcf,['fig_RS' typePrefix '_' simprefix '_2D_all_0a.fig'])
      close
    end
  end
end

% return

% plot spike phase with respect to baseline vs time
% spikeindsa = find(post(1,:)==-1e-10); name = 'LIF grid cell';
% spikeindsa = find(post(2,:)==-1e-10); name = 'Gated LIF grid cell';
[pks, spikeindsa] = find(postspikes>0); name = 'RAF grid cell';
spiketimesa = dt*spikeindsa; % ms 
spiketimesb = BaseSpikeTimes(diff(BaseSpikeTimes)>20); % ms
xs = cumsum(fdxs); % m
ys = cumsum(fdys); % m
spikeCoords = [xs(spikeindsa)' ys(spikeindsa)']; % m
spikePhases = zeros(1,length(spikeCoords));
% for each spike of vco 1
for spi=1:length(spiketimesa)
lowerbound = find(spiketimesb<spiketimesa(spi),1,'last');
upperbound = lowerbound+1;
if isempty(lowerbound) || lowerbound==length(spiketimesb)
spikePhases(spi) = NaN;
else
spikePhases(spi) = 2*pi*(spiketimesa(spi)-spiketimesb(lowerbound))/(spiketimesb(upperbound)-spiketimesb(lowerbound)); % rad
end
end
spikePhases = spikePhases(:);
figure; plot(spiketimesa,spikePhases,'.'); title(name); set(gca,'ylim',[0 2*pi]);
figure; plot(spikeCoords(:,1),spikePhases,'.'); title(name); set(gca,'ylim',[0 2*pi]);

spikeindsa = find(post(1,:)==-1e-10); name = 'LIF grid cell';
spiketimesa = dt*spikeindsa; % ms 
spiketimesb = BaseSpikeTimes(diff(BaseSpikeTimes)>20); % ms
xs = cumsum(fdxs); % m
ys = cumsum(fdys); % m
spikeCoords = [xs(spikeindsa)' ys(spikeindsa)']; % m
spikePhases = zeros(1,length(spikeCoords));
% for each spike of vco 1
for spi=1:length(spiketimesa)
lowerbound = find(spiketimesb<spiketimesa(spi),1,'last');
upperbound = lowerbound+1;
if isempty(lowerbound) || lowerbound==length(spiketimesb)
spikePhases(spi) = NaN;
else
spikePhases(spi) = 2*pi*(spiketimesa(spi)-spiketimesb(lowerbound))/(spiketimesb(upperbound)-spiketimesb(lowerbound)); % rad
end
end
spikePhases = spikePhases(:);
figure; plot(spiketimesa,spikePhases,'.'); title(name); set(gca,'ylim',[0 2*pi]);
figure; plot(spikeCoords(:,1),spikePhases,'.'); title(name); set(gca,'ylim',[0 2*pi]);
