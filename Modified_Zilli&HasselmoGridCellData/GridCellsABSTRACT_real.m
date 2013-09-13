%% Code to Generate Grid Cell firing patterns
% Requires 3 head direction cells (VCOs) responding to 3 different phase angles

%% Canibalized from Zilli&Hasselmo Version 2010 June 1


%% Comes up later
pcon=1; % Tweak to speed things up as recommended

%% Length of simulation
simdur = 1*360e3; % ms
trajdt = 1; % sampling rate for returned trajectory ** must be same as dt
if simdur>4e3
  % Take the dy and dx from actual rat (Hafting et al. 2005 mouse study)
  [dxs,dys,fdxs,fdys] = hafting_trajectory(simdur,trajdt);
else
  [dxs,dys,fdxs,fdys] = hafting_trajectory(4e3,trajdt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Neurons fire every ms
dt = 1; % ms

  nruns = 1;
  nVCOs = 3; %% Three head direction cells responding to different orientations
  
  
  
  loadC = ''; %% Not sure what you need this var for

  % generally you will want to run both models so that you have the
  % respective VCO phase differences as a visual performance measure (in
  % addition to the spatial firing itself)
  runAbstract = 1;
  runNetwork = 0;

  % generally leave these = 1, unless debugging one or the other
  runNetworkVCOs = 1;
  runNetworkBaseline = 1;

  % if true, figures will be saved to disk if a figure of the same name does
  % not already exist
  saveFigures = 1;

  % plotType = 2 means all three post plus abstract spatial firing for
  % nVCOs=1 (overrides abstractThr) (Figure S5, S9)

  % if true, plots traces of cell #1 from each VCO and the postsynaptic
  % cells (each one gets its own figure)
  % this makes Figures 7, 12, 14)
  plotVCOsAndPosts = 0;

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
  commonNoiseSTD = 0;
  useFilteredTrajectory = 1;

  
% End crazy setup stuff

%% Setup type 7 scenario
    rand('seed',1);
    randn('seed',1);
    simprefix = sprintf('filthaftingtraj_n250_1noise_3vco_%s',datestr(now,'mmmmdd'));
    excitationClass = 2;
    useGapJunctions = 0;
    useNoise = 1;
    useFilteredTrajectory = 1;
    load simple_model_RS2sn_FI_Jan09_n250.mat;
    uniqueNoiseSTD = 100*useNoise;
    g = 256;
    baselineFreq = freqs(round(-1+length(freqs)/2)); % Hz
    beta = 2; % Hz/(m/s)
 %   beta = 3;

%     % non-gated lif params:
%     tau = 5; % ms
%     weightMult = 1; 0.8;
%     membraneDecay = exp(-dt/tau);
%     postWeight=weightMult*1/ncells/(2+nVCOs);
%     baseWeight = 1.5*postWeight; 2*postWeight;
% 
%     % gated lif params:
%     basegateDur = .05; 1; 50; % ms
%     tau = 5; % (msec)
%     gatedmembraneDecay = exp(-dt/tau);
%     gatedpostWeight = 0.0014;
% 
%     % non-gated res params:
%     weightMult = 7;
%     respostWeight=0.0015; 0.002;
%     resbaseWeight = 0.003; 0.0024;
%     postc = -0.01;
%     postw = baselineFreq*2*pi/1000; % Hz

%%Not sure what this does
typePrefix = sprintf('%d',excitationClass);

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
  
%% Setup matrix for saving firing data
%savedAbstract = zeros(simdur,10);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup preferred head directions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Maximum number of grid cells
maxcol = 100;

%% Setup phase angles of preferred firing
for h = 1:(maxcol/10)
    prefHDsMatrix(h,:) = [0 2*pi/3 4*pi/3] + (h*2*pi)/(3*10*2); 
    % Rotate half way through one rotation, ie from 0 to pi/3 (radians)
end
    prefHDsMatrix = repmat(prefHDsMatrix,10,1);

%% Setup basePhi preferred firing
for p = 1:(maxcol/10)
    basePhiMatrix(p) = pi+ p * 0.02 * pi;
end

basePhiMatrix = repmat(basePhiMatrix,10,1);

%% Setup Beta variation
for p = 1:(maxcol/10)
    betaMatrix(p) = 1.5 + (0.1*p);
end

betaMatrix = repmat(betaMatrix,10,1);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for col = 1:maxcol % increment grid cells
   
   % baselineFreq = freqs(0+30*col) % Hz % Pick freqs
    
    % activity of postsynaptic I&F
    npost = 3;
    post = zeros(npost, round(simdur/dt)+1);
    if stablePost==0
      postu = zeros(npost, round(simdur/dt)+1);
      postinhib = zeros(npost, round(simdur/dt)+1);
    end
    vpost = zeros(nVCOs,round(simdur/dt)+1);
    vpostb = zeros(1,round(simdur/dt)+1);

%     % spike times of postsynaptic I&F
%     postspikes = spalloc(1, round(simdur/dt)+1, 20*simdur); % expect firing at, say, 20 Hz

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
 %% basePhi is what matters for offsetting grids
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %basePhi = 0+1*col; % baseline oscillator phase variable
    basePhi = basePhiMatrix(col);
    
    beta = betaMatrix(col);
    
    % basePhi = 0;
    naming = basePhi;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    VCOPhi = zeros(1,nVCOs); % VCO phase variable(s)
    %VCOPhi = rand(1,nVCOs);
    prefHDs = prefHDsMatrix(col,:);
%     prefHDs = [0 2*pi/3 4*pi/3 pi/3 pi 5*pi/3]; %% VCO preferred direction(s), (radians)
%     if nVCOs>length(prefHDs)
%       prefHDs = repmat(prefHDs,1,ceil(nVCOs/length(prefHDs)));
%     end
%    prefHDs = prefHDs(1:nVCOs);

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
    %%%%%%%%%%%%%%%%
    %% Main Activity
    %%%%%%%%%%%%%%%%
    

    
    tic
    while t<simdur-2*dt
      t = t+dt; % advance to next time step
      tind = 1+round(t/dt);
      
      % Print to the user every so often
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

      %% Figure out rat's HD
      speed = abs(sqrt((dx^2 + dy^2)/(dt/1000)^2)); % (m/s); abs because we include the direction via cosine later
      pos(:,tind) = [x; y]; %% Add to rat's pos matrix
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
%       
%       if ~isempty(saveVCOSpikesFileName)
%         % each line: number of spikes from 1:nPostIn on this
%         % time step for each VCO and the baseline
%         % we'll buffer the output to speed things up
%         buflin = [];
%         if runNetworkVCOs
%           for vcoi=1:nVCOs
%             buflin = [buflin sum(spikes(1:nPostIn,vcoi))];
%           end
%         end
%         if runNetworkBaseline
%           buflin = [buflin sum(spikesb(1:nPostIn))];
%         end
%         buffer = [buffer; buflin];
%         if size(buffer,1)==50
%           fid = fopen(saveVCOSpikesFileName,'a');
%           for br=1:size(buffer,1)
%             for bc=1:size(buffer,2)
%               fprintf(fid,'%d ',buffer(br,bc));
%             end
%             fprintf(fid,'\n');
%           end
%           fclose(fid);
%           buffer = [];
%         end
%       end

      state(:,tind) = [v(1,:)'; vb(1)];

      if cheatHDRectification
        spikesum = spikesum.*(cos(prefHDs-curHD)>0);
      end
      
      
%       if runningPlot
%         figure(runh);
%         plot(pos(1,:),pos(2,:)), title('network VCO spatial firing'); xlabel('Position (m)','fontsize',12), ylabel('Position (m)','fontsize',12);
%         set(gca,'xlim',[-1 1],'ylim',[-1 1]);
%         drawnow
%       end
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
    
    %%%%%%%%% End the run
    
    
%% Extra stuff
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
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Saving Abstract Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


     savedAbstract(:,col) = (abstractGrid(1:end)-min(abstractGrid)) / (max(abstractGrid)-min(abstractGrid));
%    savedAbstractlocations = find(abstractGrid(1:end)>abstractThr);
%    savedAbstract(savedAbstractlocations,col) = 1;

display(col) % display current step

% End col increment (then repeat)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%

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
  

%% Type 2 plotting
    % these in normalized coords:
    trajH = 3/4*1/2; % <1/2
    trajW = 3/4*1/3; % <1/3

    col1L = 1/4*1/3; % <1/3-trajW
    col2L = 1/3+1/5*1/3;
    col3L = 2/3+1/5*1/3;
    row2B = 1/6*1/2; % <1/3-trajH

    posskip = 100;

    % abstract:
    abstractThr = 0.9; % Threshold for showing up on abstract picture
    spikeskip = 10;
for q = 0:9
    figure
    for t = 1:10
        subplot(3,4,t)% Plot each cell individually
        plot(pos(1,1:spikeskip:length(pos)),pos(2,1:spikeskip:length(pos)),'.','MarkerSize',8);
        hold on; plot(pos(1,spikeskip*find(savedAbstract(1:spikeskip:end, (t+(10*q)))>abstractThr)),pos(2,spikeskip*find(savedAbstract(1:spikeskip:end, (t+(10*q)))>abstractThr)),'R.','MarkerSize',15)
        title('Abstract model','fontsize',12)
        xlabel('Position (m)','fontsize',12), ylabel('Position (m)','fontsize',12)
        set(gca,'xlim',[min(pos(1,:)) max(pos(1,:))],'ylim',[min(pos(2,:)) max(pos(2,:))])
        set(gca,'box','off');
        set(gca,'fontsize',12)
    end
end
    

%     if saveFigures && exist(['fig_RS' typePrefix '_' simprefix '_2D_all_0a.eps'])~=2
%       save(['fig_RS' typePrefix '_' simprefix '_2D_variables_0a.mat'])
%       %   clear pos abstractGrid vhist vcohist basehist post
%       saveas(gcf,['fig_RS' typePrefix '_' simprefix '_2D_all_0a.fig'])
%     end

%%%%%%%%%%%%%%%%%%%%%%%
%% Save data at the end
%%%%%%%%%%%%%%%%%%%%%%%
    
savedAbstract(:,col+1) = pos(1,:);
savedAbstract(:,col+2) = pos(2,:);
    
    