% Acker et al 2003 model variant neuron network with random connectivity
% with time constant input at varying levels with LIF-measured stability
% eric zilli - oct 12, 2009
%
% release version 1.0. check modeldb for updates.
%
% this source code is released into the public domain
%
% This script calculates the FI curve of a network of biophysical neurons.
%
% The Acker et al model is a biophysical model of an entorhinal cortex
% layer II stellate. The model uses voltage-shifted Hodghkin-Huxley loligo
% Na and K channels for spiking, plus a non-inactivating persistent-sodium
% current and a two-component (fast and slow) H-current which are modified
% by Acker et al. from original current models by Erik Fransen. We modify the model by
% scaling down all the conductance densities by increasing the membrane
% capacitance. This slows the minimum firing frequency of the model down so
% that the noise level can be scaled to match our biological data (ISI mean
% 0.2 s, ISI std. dev. 0.036 s).
%
% This can take a while, but it is highly parallelizable because
% each of the ~200 runs (4 Hz/0.02 (Hz/run)) are completely independent. This also makes it a
% quick matter to estimate how long running a complete FI will take: number
% of simulations times length of single simulation.
%
% Interesting questions to continue this research: What changes to
% the biophysical parameters flatten the FI curve over the same range of
% currents that might also occur in the dorsoventral gradient of grid
% cells? The M-current may well play a role in this and might be added to
% the model to explore its effects on FI.
%
% Note: this is a cleaned up version of the script used to generate the
% paper figures--it is possible errors were introduced during the cleaning!
% Feel free to contact me if there is any difficulty reproducing any
% results from the manuscript.
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
% higher values, but you may have to go back to SI_acker_model_stability_vs_params.m
% and find a new g value that works for it.
%
% Note 5: If generating FI curves for 2D grid simulations is going slow,
% try a coarser resolution then come back and run this again to fill in the
% blanks if it turns out the FI curve didn't have enough resolution (which
% will manifest as a large slope in the VCO phase differences vs. time: the
% overall slope should be essentially flat/zero if the resolution is high enough, I
% think).

clear all;

simdur = 2000; % (ms)
dt = .01; % ms

nruns=1;

% if set to 1, ISI calculations will skip the first ISI, which is likely to
% be corrupted by start-up transients
skipFirstISI = 1;

pcon = 1;
useNoise = 1;
spikeThreshold = 20; % mV

slowslow = 0;

useGapJunctions = 0;
uniqueNoiseScale = 3.44*useNoise;

g = 80;
ncells = 250;

inputMags = [-2.5:0.025:-0.5];

% if fileOut==1, this uses the output saved in the output file to create
% a .mat file of the FI curve that can be used in the 2D simulations
generateFIfile = 1;

% if fileOut==1 and generateFIfile==1, this will load a pre-existing FI
% file (if it exists) and merge the new results in, e.g. to increase the
% resolution of an existing file (new values will overwrite old values)
mergeOld = 1;

fileOut = 1;
outfilename = 'acker_model_higherCM_';
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
outfilename = [outfilename sprintf('_vary_I_%gnoise_%s.txt',useNoise,datestr(now,'mmmmdd'))];

% file to save statistics for every individual cell in network, don't
% save if indivname = ''
indivname = '';
% indivname = [strtok(outfilename,'.') '_individualcells.txt'];

Cm = 5; 1.5; % uF/cm^2
gNa = 52; % mS/cm^2
gNap = 0.5; % mS/cm^2
gH = 1.5; % mS/cm^2
gK = 11; % mS/cm^2
gLeak = 0.5; % mS/cm^2
EH = -20; % (mv)
ELeak = -65; % (mv)
EK = -90; % (mv)
ENa = 55; % (mv)

% connectivity matrix
C = (rand(ncells)<pcon)>0;
C = sparse(C.*(1-eye(ncells))); % remove autoconnections

for inputi=1:length(inputMags)
  stabilities = zeros(1,nruns);
  stabilities2 = zeros(1,nruns);
  for run=1:nruns
    % activity of postsynaptic I&F
    npost = 1;
    post = zeros(npost, round(simdur/dt)+1);

    % I&F params
    tau = 4.5; % (ms)
    membraneDecay = exp(-dt/tau);
    postWeight=1.2/ncells;

    % magnitude of noise added to all cells on each step
    commonNoiseScale = 0;

    t = 0; % current time in simulation

    % initial values
    v = -55*ones(1,ncells);
    mNa = zeros(1,ncells);
    hNa = zeros(1,ncells);
    mNap = zeros(1,ncells);
    nK = zeros(1,ncells);
    mHfast = zeros(1,ncells);
    mHslow = zeros(1,ncells);

    % save voltage of one cell over simulation:
    vhist = zeros(1,round(simdur/dt)+1);
    slowhist = zeros(1,round(simdur/dt)+1);

    % binary array indicating which cells fire on each time steps
    spikes = spalloc(ncells, round(simdur/dt)+1,ncells*simdur); % sparse to save memory in big simulations

    tic
    while t<simdur
      t = t+dt; % advance to next time step
      tind = 1+round(t/dt);

      I = inputMags(inputi); % externally applied input

      oldv = v;
      alphamNa = -0.1*(oldv+23)./(exp(-0.1*(oldv+23))-1);
      betamNa = 4*exp(-(oldv+48)/18);
      mNa = mNa + dt*(alphamNa.*(1-mNa) - betamNa.*mNa);

      alphahNa = 0.07*exp(-(oldv+37)/20);
      betahNa = 1./(exp(-0.1*(oldv+7))+1);
      hNa = hNa + dt*(alphahNa.*(1-hNa) - betahNa.*hNa);

      alphanK = -0.01*(oldv+27)./(exp(-0.1*(oldv+27))-1);
      betanK = 0.125*exp(-(oldv+37)/80);
      nK = nK + dt*(alphanK.*(1-nK) - betanK.*nK);

      alphamNap = 1./(0.15*(1+exp(-(oldv+38)/6.5)));
      betamNap = exp(-(oldv+38)/6.5)./(0.15*(1+exp(-(oldv+38)/6.5)));
      mNap = mNap + dt*(alphamNap.*(1-mNap) - betamNap.*mNap);

      minfHfast = 1./(1+exp((oldv+79.2)/9.78));
      mtauHfast = 0.51./(exp((oldv-1.7)/10) + exp(-(oldv+340)/52)) + 1;
      mHfast = mHfast + dt*((minfHfast-mHfast)./mtauHfast);

      minfHslow = 1./(1+exp((oldv+71.3)/7.9));
      if ~slowslow
        mtauHslow = 5.6./(exp((oldv-1.7)/14) + exp(-(oldv+260)/43)) + 1;
      else
        mtauHslow = 65.5813 + 248.0469*exp(-(-79.2190-oldv).^2/33.5178^2);
      end
      mHslow = mHslow + dt*((minfHslow-mHslow)./mtauHslow);

      slowhist(tind) = mHslow(1);

      vnoise = uniqueNoiseScale*randn(1,ncells);
      if useGapJunctions
        if pcon==0
          v = v + dt*(vnoise - gH*(0.65*mHfast + 0.35*mHslow).*(v-EH) - (gNap*mNap+gNa*mNa.^3.*hNa).*(v-ENa) - gK*nK.^4.*(v-EK) - gLeak*(v-ELeak) + I)/Cm;
        elseif pcon==1
          v = v + dt*(vnoise - gH*(0.65*mHfast + 0.35*mHslow).*(v-EH) - (gNap*mNap+gNa*mNa.^3.*hNa).*(v-ENa) - gK*nK.^4.*(v-EK) - gLeak*(v-ELeak) + I + g*(sum(v)-ncells*v))/Cm;
        else
          vdiffs = repmat(v',ncells,1) - repmat(v,1,ncells);
          v = v + dt*(vnoise - gH*(0.65*mHfast + 0.35*mHslow).*(v-EH) - (gNap*mNap+gNa*mNa.^3.*hNa).*(v-ENa) - gK*nK.^4.*(v-EK) - gLeak*(v-ELeak) + I + g*sum(C*vdiffs,2))/Cm;
        end
      else
        if pcon==0
          v = v + dt*(vnoise - gH*(0.65*mHfast + 0.35*mHslow).*(v-EH) - (gNap*mNap+gNa*mNa.^3.*hNa).*(v-ENa) - gK*nK.^4.*(v-EK) - gLeak*(v-ELeak) + I)/Cm;
        elseif pcon==1
          v = v + dt*(vnoise - gH*(0.65*mHfast + 0.35*mHslow).*(v-EH) - (gNap*mNap+gNa*mNa.^3.*hNa).*(v-ENa) - gK*nK.^4.*(v-EK) - gLeak*(v-ELeak) + I + g*(sum(spikes(:,tind-1),1)-spikes(:,tind-1)'))/Cm;
        else
          v = v + dt*(vnoise - gH*(0.65*mHfast + 0.35*mHslow).*(v-EH) - (gNap*mNap+gNa*mNa.^3.*hNa).*(v-ENa) - gK*nK.^4.*(v-EK) - gLeak*(v-ELeak) + I + g*(C*spikes(:,tind-1))')/Cm;
        end
      end
      vhist(tind) = v(1);

      % spikes if upward crossing spikeThreshold mV
      spikes(:,tind) = (oldv<spikeThreshold).*(v>=spikeThreshold);

      post(:,tind) = post(:,tind-1)*membraneDecay + tau*postWeight*sum(spikes(:,tind));
      post(post(:,tind)>1,tind) = -1e-10;
    end
    toc

    % calculate mean and variance of ISIs for each postsynaptic cell
    npost = 1; %size(post,1);
    ISImeans = zeros(1,npost);
    ISIstds = zeros(1,npost);
    ISIstabilities = zeros(1,nruns);

    %% calculate ISI statistics for each individual unit in the
    %% network:
    indivmeans = zeros(1,ncells);
    indivstds = zeros(1,ncells);
    indivstabilities = zeros(1,ncells);

    for ce=1:ncells % only run this once because we only save the last value in stabilities
      spikeTimes = find(spikes(ce,:)>0);
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
        end;
      end
      ISIs = fISIs;

      if skipFirstISI
        ISIs = ISIs(2:end);
      end

      % if we have 3 or fewer ISIs, let's just trash this cell:
      if length(ISIs)<=3
        indivmeans(ce) = NaN;
        indivstds(ce) = NaN;
        indivstabilities(ce) = NaN;
      else
        indivmeans(ce) = mean(ISIs);
        indivstds(ce) = std(ISIs);
        indivstabilities(ce) = 5*indivmeans(ce)^3/16/pi/pi/(eps+indivstds(ce))/(eps+indivstds(ce));
      end
    end
    % drop NaNs and sort:
    indivmeans(isnan(indivmeans)) = [];
    indivstds(isnan(indivstds)) = [];
    indivstabilities(isnan(indivstabilities)) = [];
    [indivstabilities,ix] = sort(indivstabilities);
    indivmeans = indivmeans(ix);
    indivstds = indivstds(ix);

    % actually, put NaNs back if there is nothing left:
    if isempty(indivmeans)
      % this was added so that the min/median/max of these vectors
      % wouldn't be empty and thus throw off alignment when writing
      % results to a file
      indivmeans = NaN*ones(1,ncells);
      indivstds = NaN*ones(1,ncells);
      indivstabilities = NaN*ones(1,ncells);
    end

    %% calculate ISI statistics for the postsynaptic LIF neurons
    for ce=1:npost % only run this once because we only save the last value in stabilities
      %           spikeTimes = find(postspikes(ce,:)>0);
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
        end;
      end
      ISIs = fISIs;

      if skipFirstISI
        ISIs = ISIs(2:end);
      end

      ISImeans(ce) = mean(ISIs);
      ISIstds(ce) = std(ISIs);
      ISIstabilities(ce) = 5*ISImeans(ce)^3/16/pi^2/(eps+ISIstds(ce))^2;
    end

    if length(ISIs)
      lastISIFreq = mean(1./ISIs(end/2:end));
    else
      lastISIFreq = NaN;
    end

    % save to file:
    if fileOut
      fid = fopen(outfilename,'a');
      fprintf(fid,'%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\n', ...
        ncells,pcon,I,g,simdur/1000, useNoise*uniqueNoiseScale,tau,postWeight,...
        spikeThreshold,...
        lastISIFreq,1/ISImeans(1),ISImeans(1),ISIstds(1),ISIstabilities(1),length(ISIs),...
        min(indivmeans),median(indivmeans),max(indivmeans),...
        min(indivstds),median(indivstds),max(indivstds),...
        min(indivstabilities),median(indivstabilities),max(indivstabilities),...
        length(find(post==-1e-10)));
      fclose(fid);
    end
    fprintf('%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\n', ...
      ncells,pcon,I,g,simdur/1000, useNoise*uniqueNoiseScale,tau,postWeight,...
      spikeThreshold,...
      lastISIFreq,1/ISImeans(1),ISImeans(1),ISIstds(1),ISIstabilities(1),length(ISIs), ...
      min(indivmeans),median(indivmeans),max(indivmeans),...
      min(indivstds),median(indivstds),max(indivstds),...
      min(indivstabilities),median(indivstabilities),max(indivstabilities),...
      length(find(post==-1e-10)));

    if ~isempty(indivname)
      fid = fopen(indivname,'a');
      fprintf(fid,'%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d', ...
        ncells,pcon,I,g,simdur/1000, useNoise*uniqueNoiseScale,tau,postWeight,...
        lastISIFreq,1/ISImeans(1),ISImeans(1),ISIstds(1),ISIstabilities(1),length(ISIs),...
        length(find(post==-1e-10)));
      for ce=1:length(indivmeans)
        fprintf(fid,'\t%g',indivmeans(ce));
      end
      for ce=1:ncells-length(indivmeans) % for alignment in ouput, put back the NaNs we removed
        fprintf(fid,'\tNaN');
      end
      fprintf(fid,'\t');
      for ce=1:length(indivmeans)
        fprintf(fid,'\t%g',indivstds(ce));
      end
      for ce=1:ncells-length(indivmeans)
        fprintf(fid,'\tNaN');
      end
      fprintf(fid,'\t');
      for ce=1:length(indivmeans)
        fprintf(fid,'\t%g',indivstabilities(ce));
      end
      for ce=1:ncells-length(indivmeans)
        fprintf(fid,'\tNaN');
      end
      fprintf(fid,'\n');
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