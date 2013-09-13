% izhikevich simple model neuron network with random connectivity
% with time constant input at varying levels with LIF-measured stability
% eric zilli - june 28, 2009 using simple model MATLAB code from
% Izhikevich's book modified to be a network
%
% release version 1.0. check modeldb for updates.
%
% this source code is released into the public domain
% 
% This script calculates the FI curve of a network of simple model neurons.
%
% Note: this is a cleaned up version of the script used to generate the
% paper figures--it is possible errors were introduced during the cleaning!
% Feel free to contact me if there is any difficulty reproducing any
% results from the manuscript.
%
% Note 2: This can take a while, but it is highly parallelizable because
% each of the ~200 runs (4 Hz/0.02 (Hz/run)) are completely independent. This also makes it a
% quick matter to estimate how long running a complete FI will take: number
% of simulations times length of single simulation.
%
% Note 3: Careful of the n=1 noisy case where there is no coupling to
% cancel out the noise: the resolution of the FI curve you can get will be
% very limited unless you run extremely long simulations.
%
% Note 4: If you're trying to run simulations not shown in the paper (i.e.
% going off into your own territory), I advise paying particular attention
% to the estimated stability time over the FI curve. In particular, it's
% going to decrease at higher frequencies and you're going to want even the
% lowest value to be at least as long as the length of the
% spatial trajectory you're going to run (and if the stability time is,
% say, twice as high, all the better!). Try increasing ncells if you want
% higher values, but you may have to go back and find a new g value that
% works for it.
%
% Note 5: If generating FI curves for 2D grid simulations is going slow,
% try a coarser resolution then come back and run this again to fill in the
% blanks if it turns out the FI curve didn't have enough resolution (which
% will manifest as a large slope in the VCO phase differences vs. time: the
% overall slope should be essentially flat/zero if the resolution is high enough, I
% think). If you have sampled x:a:y, try (x+a/2):a:y then (x+a/4):a/2:y
% then (x+a/8):a/4:y, etc. (this doubles the resolution each time without
% overlapping a previous run).


clear all;
warning off


% type=1 - (Figure 1) class 2 excitable, n=1, noise-free
% type=2 - (Figure 1) class 2 excitable, n=1, noisy
% type=3 - (Figure 1) class 2 excitable, n=250, gap-junction, noisy
% type=4 - (Figure 1) class 2 excitable, n=250, synaptic, noisy
% type=5 - class 1 excitable, n=1, noise-free
% type=6 - class 1 excitable, n=1, noisy
% type=7 - class 1 excitable, n=250, gap-junction, noisy
% type=8 - class 1 excitable, n=250, synaptic, noisy
% type=9 - class 2 excitable, n=1, noise-free, extended past 11 Hz for history dependence script
% type=10 - class 2 excitable, n=5000, p=0.01, synaptic, noisy
% (type=11 is commented out below but might be of interest)

for type=[10];
  disp(sprintf('type = %d',type))
  
  % simulation time step
  dt = .1; % ms

  nruns = 1;

  % whether or not per-run statistics are saved to a file
  % (from this output file the FI curve can be found in the I and mean
  % period columns)
  fileOut = 1;

  % if fileOut==1, this uses the output saved in the output file to create
  % a .mat file of the FI curve that can be used in the 2D simulations
  generateFIfile = 1;
  
  % if fileOut==1 and generateFIfile==1, this will load a pre-existing FI
  % file (if it exists) and merge the new results in, e.g. to increase the
  % resolution of an existing file (new values will overwrite old values)
  mergeOld = 1;
  
  % set number of ISIs to skip from beginning of simulation, which are likely to
  % be corrupted by start-up transients
  skipISIs = 5;
  skipFirstHalf = 0;
  pcon = 1;
  useNoise = 1;
  commonNoiseSTD = 0;

  % if not-empty, this field will be loaded to provide the connectivity
  % matrix C
  loadC = '';
%   % to make a C:
   C = (rand(ncells)<pcon)>0;
   C = sparse(C-eye(size(C))>0); % remove autoconnections
   save C_simple_nX_pX.mat C
  
  if type==1
    excitationClass = 2;
    useGapJunctions = 1;
    g = 0;
    useNoise = 0;
    uniqueNoiseSTD = 100*useNoise;
    simdur = 2500; % ms
    ncells = 1;
    inputMags = 100:0.1:128;
  elseif type==2
    excitationClass = 2;
    useGapJunctions = 1;
    g = 0;
    uniqueNoiseSTD = 100*useNoise;
    simdur = 60000; % ms
    ncells = 1;
    inputMags = 100:0.1:128; % note: due to noise this resolution is probably too high
  elseif type==3
    excitationClass = 2;
    useGapJunctions = 1;
    g = 0.1;
    uniqueNoiseSTD = 100*useNoise;
    simdur = 15000; % ms
    ncells = 250;
    inputMags = 100:0.1:128;
  elseif type==4
    excitationClass = 2;
    useGapJunctions = 0;
    g = 256;
    uniqueNoiseSTD = 100*useNoise;
    simdur = 15000; % ms
    ncells = 250;
    inputMags = 2.^[6.5:.005:8];
  elseif type==5
    excitationClass = 1;
    useGapJunctions = 1;
    g = 0;
    useNoise = 0;
    uniqueNoiseSTD = 100*useNoise;
    simdur = 2500; % ms
    ncells = 1;
    inputMags = [60:.1:80];
  elseif type==6
    excitationClass = 1;
    useGapJunctions = 1;
    g = 0;
    uniqueNoiseSTD = 130*useNoise;
    simdur = 60000; % ms
    ncells = 1;
    inputMags = [60:.1:80]; % note: due to noise this resolution is probably too high
  elseif type==7
    excitationClass = 1;
    useGapJunctions = 1;
    g = 0.1;
    uniqueNoiseSTD = 130*useNoise;
    simdur = 15000; % ms
    ncells = 250;
    inputMags = 80; [56:.1:80];
  elseif type==8
    excitationClass = 1;
    useGapJunctions = 0;
    g = 128;
    uniqueNoiseSTD = 130*useNoise;
    simdur = 15000; % ms
    ncells = 250;
    inputMags = [70:.1:120];
  elseif type==9
    excitationClass = 2;
    useGapJunctions = 1;
    g = 0;
    useNoise = 0;
    uniqueNoiseSTD = 100*useNoise;
    simdur = 2500; % ms
    ncells = 1;
    inputMags = 100:0.1:133;
  elseif type==10
    excitationClass = 2;
    useGapJunctions = 0;
    g = 256;
    pcon = 0.01;
    loadC = 'C_simple_sn_n5000_p01.mat';
    useNoise = 1;
    uniqueNoiseSTD = 100*useNoise;
    simdur = 10000; % ms
    ncells = 5000;
    inputMags = [100.5:1:125];
  end
  
%% this case may be of interest to people hoping to eliminate
%% history-dependence on firing frequency:
%   elseif type==11
%     excitationClass = -1; % -1 means type 1 with u resetting instead of incrementing
%     useGapJunctions = 1;
%     g = 0.1;
%     uniqueNoiseSTD = 130*useNoise;
%     I = 62;
%     inputMags = [58:.05:86];


  
  outfilename = ['simple_model_RS' sprintf('%d',excitationClass)];
  if ncells>1
    if useGapJunctions
      outfilename = [outfilename 'g'];
    else
      outfilename = [outfilename 's'];
    end
  end
  if useNoise
    outfilename = [outfilename 'n'];
  end
  if type==9
    outfilename = [outfilename '_ext'];
  end
  outfilename = [outfilename '_FI_' datestr(now,'mmmdd') '_n' num2str(ncells) '.txt'];

  outfilename = 'simple_model_RS2sn_FI_Jul11_n5000.txt';
  disp('outfilename is being overridden!')
  
  disp(outfilename)
  
  % file to save statistics for every individual cell in network, don't
  % save if indivname = ''
  indivname = '';
  % indivname = [strtok(outfilename,'.') '_individualcells.txt'];

  Cf=100; vr=-60; vt=-40; k=0.7; % parameters used for RS
  a=0.03; c=-50; d=100; % neocortical pyramidal neurons
  vpeak=35; % spike cutoff
  if abs(excitationClass)==1
    b=-2;
  elseif abs(excitationClass)==2
    b = 2;
  end

  if excitationClass<0
    d = 60; % reset to 60 instead of 100
  end

  n=round(simdur/dt); % number of simulation steps

  if ~isempty(loadC)
    load(loadC)
  else
    % connectivity matrix
    C = (rand(ncells)<pcon)>0;
    C = sparse(C.*(1-eye(ncells))); % remove autoconnections
  end

  for inputi=1:length(inputMags)

    stabilities = zeros(1,nruns);
    stabilities2 = zeros(1,nruns);
    for run=1:nruns
      % activity of postsynaptic I&F
      npost = 1;
      post = zeros(npost, round(simdur/dt)+1);

      % spike times of postsynaptic I&F
      postspikes = spalloc(npost, round(simdur/dt)+1, 20*simdur); % expect firing at 20 Hz

      % I&F params
      tau = 4.5; % (ms)
      membraneDecay = exp(-dt/tau);
      postWeight=1.2/ncells;
      if type==10
        postWeight = 0.0006;
      end

      % magnitude of noise added to all cells on each step
      commonNoiseSTD = 0;
      commonNoise = 0;

      t = 0; % current time in simulation
      v = vr*ones(ncells,1);
      u = 0*v; % initial values

      % save state of one cell over simulation:
      state = zeros(2,simdur/dt);

      % binary array indicating which cells fire on each time steps
      spikes = zeros(ncells, 1);

      tic
      while t<simdur
        t = t+dt; % advance to next time step
        tind = 1+round(t/dt);

        if commonNoiseSTD
          commonNoise = commonNoiseSTD*randn;
        end

        I = inputMags(inputi); % externally applied input

        oldv = v;
        oldu = u;
        if useGapJunctions
          if pcon==1
            % if all-to-all connected, the gap-junctions are equivalent
            % to sum(v)-ncells*v
            v = oldv + dt*(k*(oldv-vr).*(oldv-vt) - oldu + I + commonNoise + uniqueNoiseSTD*randn(ncells,1) + g*(sum(v)-ncells*v))/Cf;
          else
            % make matrix of differences in membrane potential
            vdiffs = C.*(repmat(oldv',ncells,1) - repmat(oldv,1,ncells));
            v = oldv + dt*(k*(oldv-vr).*(oldv-vt) - oldu + I + commonNoise + uniqueNoiseSTD*randn(ncells,1) + g*sum(vdiffs,2))/Cf;
            %% slowest way:
            % v = oldv + dt*(k*(oldv-vr).*(oldv-vt) - oldu + I + commonNoise + uniqueNoiseSTD*randn(ncells,1) + g*sum(C.*vdiffs,2))/Cf;
          end
        else
          if pcon==1
            v = oldv + dt*(k*(oldv-vr).*(oldv-vt) - oldu + I + commonNoise + uniqueNoiseSTD*randn(ncells,1) + g*(sum(spikes,1)-spikes))/Cf;
          else
            v = oldv + dt*(k*(oldv-vr).*(oldv-vt) - oldu + I + commonNoise + uniqueNoiseSTD*randn(ncells,1) + g*(C*double(spikes)))/Cf;
          end
        end
        u = oldu + dt*a*(b*(oldv-vr)-oldu);

        % save and reset spikes when v>=vpeak
        spikes = v>=vpeak;
        if excitationClass>0
          u(v>=vpeak) = u(v>=vpeak)+d;
        else
          u(v>=vpeak) = d;
        end
        oldv(v>=vpeak) = vpeak;
        v(v>=vpeak) = c;

        state(:,tind) = [oldv(1); oldu(1)];

        post(:,tind) = post(:,tind-1)*membraneDecay + tau*postWeight*sum(spikes);
        % to save memory, instead of keeping a separate matrix for the
        % times the postsynaptic cell fired, we reset it to
        % ever-so-slightly below 0 to identify when it has spiked
        post(post(:,tind)>1,tind) = -1e-10;
      end
      toc
      % calculate mean and variance of ISIs for each postsynaptic cell
      npost = 1; %size(post,1);
      ISImeans = zeros(1,npost);
      ISIstds = zeros(1,npost);
      ISIstabilities = zeros(1,nruns);
      
      %% calculate ISI statistics for the postsynaptic LIF neurons
      for ce=1:npost % only run this once because we only save the last value in stabilities
        spikeTimes = find(post(ce,:)==-1e-10);
        ISIs = dt*diff(spikeTimes)/1000; % convert to seconds

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
          end
        end
        ISIs = fISIs;

        if skipISIs
          ISIs = ISIs(skipISIs:end);
        end
        if skipFirstHalf
          ISIs = ISIs(end/2:end);
        end


        ISImeans(ce) = mean(ISIs);
        ISIstds(ce) = std(ISIs);
        ISIstabilities(ce) = 5*ISImeans(ce)^3/16/pi^2/(eps+ISIstds(ce))^2;
      end

      % save to file:
      if fileOut
        fid = fopen(outfilename,'a');
        fprintf(fid,'%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\n', ...
          ncells,pcon,I,g,simdur/1000, useNoise*uniqueNoiseSTD,tau,postWeight,...
          Cf,vr,vt,k,a,b,c,d,vpeak,...
          mean(1./ISIs),ISImeans(1),ISIstds(1),ISIstabilities(1),length(ISIs),...
          length(spikeTimes));
        fclose(fid);
      end
      fprintf('%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\n', ...
        ncells,pcon,I,g,simdur/1000, useNoise*uniqueNoiseSTD,tau,postWeight,...
        Cf,vr,vt,k,a,b,c,d,vpeak,...
        mean(1./ISIs),ISImeans(1),ISIstds(1),ISIstabilities(1),length(ISIs), ...
        length(spikeTimes));

      if ~isempty(indivname)
        fid = fopen(indivname,'a');
        fprintf(fid,'%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\n', ...
          ncells,pcon,I,g,simdur/1000, useNoise*uniqueNoiseSTD,tau,postWeight,...
          mean(1./ISIs),ISImeans(1),ISIstds(1),ISIstabilities(1),length(ISIs),...
          length(spikeTimes));
        fclose(fid);
      end
    end

  end

  % now re-open the text file we just wrote and save the relevant bits as a
  % .mat file
  if fileOut && generateFIfile
    clear freqs currents;
    FIfilename = [strtok(outfilename,'.') '.mat'];
    data = textread(outfilename);

    if mergeOld && exist(FIfilename)==2
      load(FIfilename);
    else
      freqs = []; currents = [];
    end
    freqs = [freqs(:); data(:,18)];
    currents = [currents(:); data(:,3)];
    
    % in case of duplicates, use the newer one:
    [currents,m,n] = unique(currents);
    freqs = freqs(m);
    [freqs,m,n] = unique(freqs);
    currents = currents(m);

    if any(diff(freqs)<0)
      warning('badness happened! freqs is non-monotonic! generally noise causes this: increasing simdur may help; else use parameters to lower effects of noise. dropping bad samples')
      bads = find(diff(freqs)<0);
      freqs(bads+1) = [];
      currents(bads+1) = [];
    end
    
    save(FIfilename,'freqs','currents');
  end
end

warning on
return

%% To generate Figure 4 itself once the FI curves are generated
ms = 16;
figure; hold on;
load simple_model_RS2gn_FI_Jan05_n250.mat;
% load RS4gnFI_250b;
plot(currents,freqs,'b.-','MarkerSize',ms);
load simple_model_RS2n_FI_Jan05_n1.mat;
% load RS4nFI_1b;
plot(currents,freqs,'r.-','MarkerSize',ms);
load simple_model_RS2_FI_Jan05_n1.mat;
% load RS4FI_1b;
plot(currents,freqs,'k.-','MarkerSize',10); 
load simple_model_RS2sn_FI_Jan05_n250.mat;
% load RS4snFI_250b; % load last to use in other fig
plot(currents,freqs,'m.-','MarkerSize',14);
xlabel('Input current magnitude','fontsize',12); ylabel('Network frequency (Hz)','fontsize',12)
set(gcf,'position',[82 404 2*257 302])
set(gca,'xlim',[100 300])
set(gca,'ylim',[4 11])
set(gca,'fontsize',12)
if exist('fig4a_RS2all_freq_vs_I.eps')~=2
  print_eps('fig4a_RS2all_freq_vs_I.eps')
  saveas(gcf,'fig4a_RS2all_freq_vs_I.fig')
end

figure;
beta=2; omegab=freqs(1)/2+freqs(end)/2;
vels = (freqs-omegab)/beta;
plot(vels,currents);
xlabel('Velocity (m/s)','fontsize',12); ylabel('Input current magnitude','fontsize',12);
set(gcf,'position',[82 404 1*257 302])
set(gca,'xlim',[min(vels) max(vels)]);
set(gca,'ylim',[min(currents) max(currents)]);
set(gca,'fontsize',12)
if exist('fig4b_RS2all_vel_vs_I.eps')~=2
  print_eps('fig4b_RS2all_vel_vs_I.eps')
  saveas(gcf,'fig4b_RS2all_vel_vs_I.fig')
end

