% Acker et al 2003 model variant neuron network with random connectivity
% 2d grid simulation
% eric zilli - october 13, 2009
%
% release version 1.0. check modeldb for updates.
%
% this source code is released into the public domain
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
% This script uses experimentally recorded rat trajectory (given by the
% hafting_trajectory function) as input to multiple VCOs which can
% interfere to produce spatial firing of a postsynaptic cell which,
% hopefully, is arrayed in the classical hexagonal grid pattern.
%
% Running this script requires that an FI curve of the network has been
% generated (e.g. by SI_acker_model_FI_relation.m) for a cell or network
% with identical parameters to those used here.
%
% This is generally way slow: behavioral timescales are big and neural
% timescales are tiny. Hints: using pcon=1 lets the script bypass the
% connectivity matrix which really speeds things up. If the resonant postsynaptic
% cell starts spiking too much the script will slow way down because the
% vector of spiketimes of the postsynaptic cell is a pre-allocated sparse
% vector. When adjusting postsynaptic cell parameters, start with
% short simulation to act as a sort of sanity check.
% For 2 VCO, 160 s simulations, my desktop takes about
% 16 hours. It is easily possible there are simple speedups that would make
% this more wieldy and the interested reader might learn a lot about the
% model by trying to find them!
%
% Interesting questions to continue this line of research:
% Can we analytically find appropriate parameters for the
% postsynaptic models given requirements such as: the field diameter should
% be about half the field spacing? (Here's a start: that means on the edge
% of a field between two neighboring fields,
% all the VCOs are at phase zero except for one which is at pi/2,
% corresponding to a fixed time given a baseline frequency, which means the
% decayed sum of the inputs from all but one VCO (decayed over the time corresponding
% to pi/2 radians, or one quarter the duration of the oscillator period)
% plus the input from one VCO should equal the firing threshold [but a tau
% too high will carry over activity from the previous cycle!]). If noise is
% correlated between the VCO inputs, what is the effect on the model and on
% the stability of the oscillations between the VCOs? Zilli et al 2009 only
% examined this question for abstract oscillators, not for these more
% complex models where the effects may be worse.
%
% More possibilities are allowed than the single set of preset parameters here,
% but postsynaptic parameters might need to be adjusted.
% Note that parameters will need to change if the baseline frequency
% changes.
%
% Note: this is a cleaned up version of the script used to generate the
% paper figures--it is possible errors were introduced during the cleaning!
% Feel free to contact me if there is any difficulty reproducing any
% results from the manuscript.
%
% This generates Figure 10 and has just that one "preset".

clear all;
% grab trajectory: this is a bit slow so best to do it once then comment it
% out

simdur = 240*1e3; % ms
trajdt = 0.1; % sampling rate for trajectory; ms
[dxs,dys,fdxs,fdys] = hafting_trajectory(simdur,trajdt);
clear dxs dys

dt = .01; % ms

useFilteredTrajectory = 1;

% only store 1 out of every storeskip timesteps
storeskip = 20;

% to extend a run, enter the total previous run duration and increase
% simdur above to the new end-of-simulation-time (and comment out the "clear
% all" line!)
extendPreviousRun = 0;
previousDur = 12000; % (ms)
if extendPreviousRun
  %% this may not work anymore! try on a really short simulation first if
  %% you want to use this
  'extending previous run instead of overwriting!'
end

nruns = 1;
nVCOs = 2;

runAbstract = 1;
runNetwork = 1;

% if ==1, abstract output will be plotted, ==0 means resonator output
% plotted in its place, ==-1 means nothing in those spaces
showAbstractFig = -1;

saveFigures = 0;

% live plot of network spatial activity
runningPlot = 0;

if runningPlot
  runh = figure;
end

useGapJunctions = 0;
slowslow = 0;
ncells = 250;
pcon = 1;
useNoise = 1;
commonNoiseScale = 0;
uniqueNoiseScale = 3.44*useNoise;
g = 80; % syn
I = -2.35;
spikeThreshold = 20; % mV

load Acker_sn_FI_n250.mat;

% prefix used for saving figures
simprefix = sprintf('filttraj_n%d_%gnoise_%ds_%s',ncells,useNoise,round(simdur/1000),datestr(now,'mmmmdd'));
if ~useFilteredTrajectory
  simprefix = ['un' simprefix];
end

% shared grid variables: freq = baselineFreq+beta*speed*(cos(prefHD-curHD));
baselineFreq = freqs(round(length(freqs)/2));
beta = 2; % Hz/(m/s)

baselineMult = 1;

% non-gated lif params:
tau = 40; % ms
membraneDecay = exp(-dt/tau);
posthWeight=0.15/ncells;
baseWeight = 0.8/ncells;

% gated lif params:
% I&F params
basegateDur = 40; % ms
tau = 50; % (msec)
gatedmembraneDecay = exp(-dt/tau);
weightMult = 1.2;
gatedposthWeight = weightMult/ncells/nVCOs;

% non-gated res params:
weightMult = 100;
resposthWeight=weightMult*1/ncells/(1+nVCOs);
resbaseWeight = 2*posthWeight;
posthc = -0.005;
posthw = baselineFreq*2*pi/1000; % Hz (times 2pi for angular freq, /1000 to convert to Hz

%% can downsample FI curve to test the needed level of precision
%   freqs = freqs(1:50:end);
%   currents = currents(1:50:end);

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
if pcon~=1
  C = (rand(ncells)<pcon)>0;
  C = sparse(C.*(1-eye(ncells))); % remove autoconnections
else
  C = 0;
end
stabilities = zeros(1,nruns);
stabilities2 = zeros(1,nruns);
for run=1
  % activity of posthsynaptic I&F
  npost = 3;
  if ~extendPreviousRun
    posth = zeros(npost, round(simdur/dt/storeskip)+1);
    post = zeros(npost,1);

    % spike times of posthsynaptic I&F
    posthspikes = spalloc(npost, round(simdur/dt/storeskip)+1, 6*simdur); % expect firing at 6 Hz

    t = 0; % current time in simulation
    tind = 1;
    storeind = 1;

    % initial values
    v = -55*ones(ncells,nVCOs);
    mNa = zeros(ncells,nVCOs);
    hNa = zeros(ncells,nVCOs);
    mNap = zeros(ncells,nVCOs);
    nK = zeros(ncells,nVCOs);
    mHfast = zeros(ncells,nVCOs);
    mHslow = zeros(ncells,nVCOs);

    vb = -55*ones(ncells,1);
    mNab = zeros(ncells,1);
    hNab = zeros(ncells,1);
    mNapb = zeros(ncells,1);
    nKb = zeros(ncells,1);
    mHfastb = zeros(ncells,1);
    mHslowb = zeros(ncells,1);

    % save voltage of one cell over simulation:
    vhist = zeros(nVCOs, round(simdur/dt/storeskip)+1);
    vhistb = zeros(1, round(simdur/dt/storeskip)+1);
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
    prefHDs = [0 2*pi/3]; % VCO preferred direction(s), (radians)
    abstractGrid = zeros(1, round(simdur/dt/storeskip)+1);
    basehist = zeros(1,round(simdur/dt));
    vcohist = zeros(nVCOs,round(simdur/dt));

    x = 0; % m
    dx = 0; % m
    y = 0; % m
    dy = 0; % m

    pos = zeros(2, round(simdur/dt/storeskip)+1);
    baseI = interp1(freqs,currents,baselineFreq,'linear');
    I = baseI;
    gatedI = baseI*ones(1,nVCOs);
    ISIs = zeros(1,nVCOs);

    speed = 0; % m/s
    gatedspeeds = zeros(1,nVCOs);

    % external path integration to hold I constant:
    XatLastSpike = zeros(1,nVCOs);
    YatLastSpike = zeros(1,nVCOs);

    basegateCount = 0;
  else
    posth(npost,round(simdur/dt)+1) = 0; % implicit extension
    abstractGrid(1,round(simdur/dt/storeskip)+1) = 0; % implicit extension
    vhist(nVCOs,round(simdur/dt/storeskip)+1) = 0; % implicit extension
    vhistb(nVCOs,round(simdur/dt/storeskip)+1) = 0; % implicit extension
    pos(2,round(simdur/dt/storeskip)+1) = 0; % implicit extension
  end

  tic
  while t<simdur
    t = t+dt; % advance to next time step
    tind = 1+tind;

    if mod(t,round(simdur)/100)<=dt
      disp(sprintf('t = %d, %g elapsed',t,toc));
    end

    %%% filtered, interpolated trajectory (scale dx/dy by mismatch between
    %%% time steps to properly scale speed)
    if useFilteredTrajectory
      dx = fdxs(ceil(t/trajdt))*(dt/trajdt);
      dy = fdys(ceil(t/trajdt))*(dt/trajdt);
    else
      dx = dxs(ceil(t/trajdt))*(dt/trajdt);
      dy = dys(ceil(t/trajdt))*(dt/trajdt);
    end

    x = x+dx;
    y = y+dy;

    speed = abs(sqrt((dx^2 + dy^2)/(dt/1000)^2)); % (m/s)
    curHD = atan2(dy,dx);
    if mod(tind,storeskip)==1
      storeind=storeind+1;
      pos(:,storeind) = [x; y];
    end

    if runAbstract
      %% abstract grid model: (divide dt by 1000 because freqs are in Hz
      %% but dt is in ms)
      basePhi = basePhi + dt/1000*2*pi*baselineFreq;
      VCOAngularFreqs = 2*pi*(baselineFreq+beta*speed*(cos(prefHDs-curHD)));
      VCOPhi = VCOPhi + dt/1000*VCOAngularFreqs;
      basehist(tind) = basePhi;
      vcohist(:,tind) = VCOPhi;
    end
    if mod(tind,storeskip)==1
      abstractGrid(storeind) = nVCOs*cos(basePhi)+sum(cos(VCOPhi));
    end

    if runNetwork
      %% active VCOs
      for vco=1:nVCOs
        oldv = v(:,vco);
        alphamNa = -0.1*(oldv+23)./(exp(-0.1*(oldv+23))-1);
        betamNa = 4*exp(-(oldv+48)/18);
        mNa(:,vco) = mNa(:,vco) + dt*(alphamNa.*(1-mNa(:,vco)) - betamNa.*mNa(:,vco));

        alphahNa = 0.07*exp(-(oldv+37)/20);
        betahNa = 1./(exp(-0.1*(oldv+7))+1);
        hNa(:,vco) = hNa(:,vco) + dt*(alphahNa.*(1-hNa(:,vco)) - betahNa.*hNa(:,vco));

        alphanK = -0.01*(oldv+27)./(exp(-0.1*(oldv+27))-1);
        betanK = 0.125*exp(-(oldv+37)/80);
        nK(:,vco) = nK(:,vco) + dt*(alphanK.*(1-nK(:,vco)) - betanK.*nK(:,vco));

        alphamNap = 1./(0.15*(1+exp(-(oldv+38)/6.5)));
        betamNap = exp(-(oldv+38)/6.5)./(0.15*(1+exp(-(oldv+38)/6.5)));
        mNap(:,vco) = mNap(:,vco) + dt*(alphamNap.*(1-mNap(:,vco)) - betamNap.*mNap(:,vco));

        minfHfast = 1./(1+exp((oldv+79.2)/9.78));
        mtauHfast = 0.51./(exp((oldv-1.7)/10) + exp(-(oldv+340)/52)) + 1;
        mHfast(:,vco) = mHfast(:,vco) + dt*((minfHfast-mHfast(:,vco))./mtauHfast);

        minfHslow = 1./(1+exp((oldv+71.3)/7.9));
        if ~slowslow
          mtauHslow = 5.6./(exp((oldv-1.7)/14) + exp(-(oldv+260)/43)) + 1;
        else
          mtauHslow = 65.5813 + 248.0469*exp(-(-79.2190-oldv).^2/33.5178^2);
        end
        mHslow(:,vco) = mHslow(:,vco) + dt*((minfHslow-mHslow(:,vco))./mtauHslow);

        desiredFreq = baselineFreq+beta*speed*(cos(prefHDs(vco)-curHD));
        lowind = max(find(freqs<desiredFreq));
        highind = min(find(freqs>desiredFreq));
        proportion = (desiredFreq-freqs(lowind))/(freqs(highind)-freqs(lowind));
        I = currents(lowind) + proportion*(currents(highind)-currents(lowind));

        vnoise = uniqueNoiseScale*randn(ncells,1);
        if useGapJunctions
          if pcon==0
            v(:,vco) = v(:,vco) + dt*(vnoise - gH*(0.65*mHfast(:,vco) + 0.35*mHslow(:,vco)).*(v(:,vco)-EH) - (gNap*mNap(:,vco)+gNa*mNa(:,vco).^3.*hNa(:,vco)).*(v(:,vco)-ENa) - gK*nK(:,vco).^4.*(v(:,vco)-EK) - gLeak*(v(:,vco)-ELeak) + I)/Cm;
          elseif pcon==1
            v(:,vco) = v(:,vco) + dt*(vnoise - gH*(0.65*mHfast(:,vco) + 0.35*mHslow(:,vco)).*(v(:,vco)-EH) - (gNap*mNap(:,vco)+gNa*mNa(:,vco).^3.*hNa(:,vco)).*(v(:,vco)-ENa) - gK*nK(:,vco).^4.*(v(:,vco)-EK) - gLeak*(v(:,vco)-ELeak) + I + g*(sum(v(:,vco))-ncells*v(:,vco)))/Cm;
          else
            vdiffs = repmat(v(:,vco)',ncells,1) - repmat(v(:,vco),1,ncells);
            v(:,vco) = v(:,vco) + dt*(vnoise - gH*(0.65*mHfast(:,vco) + 0.35*mHslow(:,vco)).*(v(:,vco)-EH) - (gNap*mNap(:,vco)+gNa*mNa(:,vco).^3.*hNa(:,vco)).*(v(:,vco)-ENa) - gK*nK(:,vco).^4.*(v(:,vco)-EK) - gLeak*(v(:,vco)-ELeak) + I + g*sum(C*vdiffs,2))/Cm;
          end
        else
          if pcon==0
            v(:,vco) = v(:,vco) + dt*(vnoise - gH*(0.65*mHfast(:,vco) + 0.35*mHslow(:,vco)).*(v(:,vco)-EH) - (gNap*mNap(:,vco)+gNa*mNa(:,vco).^3.*hNa(:,vco)).*(v(:,vco)-ENa) - gK*nK(:,vco).^4.*(v(:,vco)-EK) - gLeak*(v(:,vco)-ELeak) + I)/Cm;
          elseif pcon==1
            v(:,vco) = v(:,vco) + dt*(vnoise - gH*(0.65*mHfast(:,vco) + 0.35*mHslow(:,vco)).*(v(:,vco)-EH) - (gNap*mNap(:,vco)+gNa*mNa(:,vco).^3.*hNa(:,vco)).*(v(:,vco)-ENa) - gK*nK(:,vco).^4.*(v(:,vco)-EK) - gLeak*(v(:,vco)-ELeak) + I + g*(sum(spikes(:,vco))-spikes(:,vco)')')/Cm;
          else
            v(:,vco) = v(:,vco) + dt*(vnoise - gH*(0.65*mHfast(:,vco) + 0.35*mHslow(:,vco)).*(v(:,vco)-EH) - (gNap*mNap(:,vco)+gNa*mNa(:,vco).^3.*hNa(:,vco)).*(v(:,vco)-ENa) - gK*nK(:,vco).^4.*(v(:,vco)-EK) - gLeak*(v(:,vco)-ELeak) + I + g*(C*spikes(:,vco))')/Cm;
          end
        end
        % spikes if upward crossing spikeThreshold mV
        spikes(:,vco) = (oldv<spikeThreshold).*(v(:,vco)>=spikeThreshold);

        if any(spikes(:,vco))
          % faster I hope:
          if mod(VCOinds(vco),1000)==0
            VCOSpikeTimes{vco}(length(VCOSpikeTimes{vco})+1000) = 0;
          end
          VCOinds(vco) = VCOinds(vco)+1;
          VCOSpikeTimes{vco}(VCOinds(vco)) = t;

          % now save current position for end of this cycle:
          XatLastSpike(vco) = x;
          YatLastSpike(vco) = y;
        end
      end

      if mod(tind,storeskip)==1
        vhist(:,storeind) = v(1,:);
        if any(isnan(v(1,:)))
          error('it is bad again')
        end
      end

      %% baseline oscillator
      oldvb = vb;
      alphamNa = -0.1*(oldvb+23)./(exp(-0.1*(oldvb+23))-1);
      betamNa = 4*exp(-(oldvb+48)/18);
      mNab = mNab + dt*(alphamNa.*(1-mNab) - betamNa.*mNab);

      alphahNa = 0.07*exp(-(oldvb+37)/20);
      betahNa = 1./(exp(-0.1*(oldvb+7))+1);
      hNab = hNab + dt*(alphahNa.*(1-hNab) - betahNa.*hNab);

      alphanK = -0.01*(oldvb+27)./(exp(-0.1*(oldvb+27))-1);
      betanK = 0.125*exp(-(oldvb+37)/80);
      nKb = nKb + dt*(alphanK.*(1-nKb) - betanK.*nKb);

      alphamNap = 1./(0.15*(1+exp(-(oldvb+38)/6.5)));
      betamNap = exp(-(oldvb+38)/6.5)./(0.15*(1+exp(-(oldvb+38)/6.5)));
      mNapb = mNapb + dt*(alphamNap.*(1-mNapb) - betamNap.*mNapb);

      minfHfast = 1./(1+exp((oldvb+79.2)/9.78));
      mtauHfast = 0.51./(exp((oldvb-1.7)/10) + exp(-(oldvb+340)/52)) + 1;
      mHfastb = mHfastb + dt*((minfHfast-mHfastb)./mtauHfast);

      minfHslow = 1./(1+exp((oldvb+71.3)/7.9));
      if ~slowslow
        mtauHslow = 5.6./(exp((oldvb-1.7)/14) + exp(-(oldvb+260)/43)) + 1;
      else
        mtauHslow = 65.5813 + 248.0469*exp(-(-79.2190-oldvb).^2/33.5178^2);
      end
      mHslowb = mHslowb + dt*((minfHslow-mHslowb)./mtauHslow);

      vnoise = uniqueNoiseScale*randn(ncells,1);
      if useGapJunctions
        if pcon==0
          vb = vb + dt*(vnoise - gH*(0.65*mHfastb + 0.35*mHslowb).*(vb-EH) - (gNap*mNapb+gNa*mNab.^3.*hNab).*(vb-ENa) - gK*nKb.^4.*(vb-EK) - gLeak*(vb-ELeak) + baseI)/Cm;
        elseif pcon==1
          vb = vb + dt*(vnoise - gH*(0.65*mHfastb + 0.35*mHslowb).*(vb-EH) - (gNap*mNapb+gNa*mNab.^3.*hNab).*(vb-ENa) - gK*nKb.^4.*(vb-EK) - gLeak*(vb-ELeak) + baseI + g*(sum(vb)-ncells*vb))/Cm;
        else
          vdiffs = repmat(vb',ncells,1) - repmat(vb,1,ncells);
          vb = vb + dt*(vnoise - gH*(0.65*mHfastb + 0.35*mHslowb).*(vb-EH) - (gNap*mNapb+gNa*mNab.^3.*hNab).*(vb-ENa) - gK*nKb.^4.*(vb-EK) - gLeak*(vb-ELeak) + baseI + g*sum(C*vdiffs,2))/Cm;
        end
      else
        if pcon==0
          vb = vb + dt*(vnoise - gH*(0.65*mHfastb + 0.35*mHslowb).*(vb-EH) - (gNap*mNapb+gNa*mNab.^3.*hNab).*(vb-ENa) - gK*nKb.^4.*(vb-EK) - gLeak*(vb-ELeak) + baseI)/Cm;
        elseif pcon==1
          vb = vb + dt*(vnoise - gH*(0.65*mHfastb + 0.35*mHslowb).*(vb-EH) - (gNap*mNapb+gNa*mNab.^3.*hNab).*(vb-ENa) - gK*nKb.^4.*(vb-EK) - gLeak*(vb-ELeak) + baseI + g*(sum(spikesb,1)-spikesb')')/Cm;
        else
          vb = vb + dt*(vnoise - gH*(0.65*mHfastb + 0.35*mHslowb).*(vb-EH) - (gNap*mNapb+gNa*mNab.^3.*hNab).*(vb-ENa) - gK*nKb.^4.*(vb-EK) - gLeak*(vb-ELeak) + baseI + g*(C*spikesb)')/Cm;
        end
      end

      % spikes if upward crossing spikeThreshold mV
      spikesb = (oldvb<spikeThreshold).*(vb>=spikeThreshold);

      if any(spikesb)
        if mod(Baseind,1000)==0
          BaseSpikeTimes(length(BaseSpikeTimes)+1000) = 0;
        end
        Baseind = Baseind+1;
        BaseSpikeTimes(Baseind) = t;

        % if baseline spikes, active VCOs can influence posth for
        % basegateCount steps
        basegateCount = basegateDur/dt;
      else
        basegateCount = basegateCount-1;
      end


      if mod(tind,storeskip)==1
        vhistb(:,storeind) = vb(1);
      end

      % while basegateCount>0, the gate is open
      if basegateCount>0
        basegate = 1;
      else
        basegate = 0;
      end

      spikesum = sum(sum(spikes));
      spikebsum = sum(spikesb);
      % posth cell 1 = non-gated lif
      post(1) = post(1)*membraneDecay + posthWeight*spikesum+baselineMult*baseWeight*spikebsum;
      if post(1)>1
        posthspikes(1,storeind) = posthspikes(1,storeind) + post(1)>1;
        post(1) = 0;
      end
      % posth cell 2 = gated lif
      post(2) = post(2)*gatedmembraneDecay + basegate*gatedposthWeight*spikesum;
      if post(2)>1
        posthspikes(2,storeind) = posthspikes(2,storeind) + post(2)>1;
        post(2) = 0;
      end
      % posth cell 3 = non-gated resonant
      post(3) = post(3) + dt*((posthc+posthw*i)*post(3) + resposthWeight*spikesum+baselineMult*resbaseWeight*spikebsum);
      %       if mod(tind,storeskip)==1
      %% since this doesn't spike a lot, let's add up the spikes instead
      %% of periodically sampling whether the cell is spiking
      posthspikes(3,storeind) = posthspikes(3,storeind) + real(post(3))>1;
      %       end
      if real(post(3))>1
        post(3) = i;
      end

      if mod(tind,storeskip)==1
        posth(:,storeind) = post;
      end

      if runningPlot
        figure(runh);
        plot(pos(1,:),pos(2,:),'.',pos(1,find(posthspikes>0)),pos(2,find(posthspikes>0)),'R.'), title('network VCO spatial firing'); xlabel('position (m)'), ylabel('position (m)');
        set(gca,'xlim',[-1 1],'ylim',[-1 1]);
        drawnow
      end
    end
  end
  toc

  for vco=1:nVCOs
    VCOSpikeTimes{vco}(VCOSpikeTimes{vco}==0) = [];
  end
  BaseSpikeTimes(BaseSpikeTimes==0) = [];

  clear C

  if runNetwork
    for ce=1:npost
      spikeTimes = find(posthspikes(ce,:)>0);

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

if simdur>4000
  plotduration = 4000; % ms
else
  plotduration  = simdur;
end
figure; plot(vhist(:,1:(plotduration/dt/storeskip))');
hold on; plot(vhistb(1:(plotduration/dt/storeskip))','r');
set(gcf,'position',[82 404 3*257 2/3*302])
set(gca,'fontsize',12,'box','off')
xlabel('Time (ms)','fontsize',12); ylabel('Voltage (mV)','fontsize',12)
if saveFigures && exist(['fig_Acker' typePrefix '_' simprefix '_2D_vcos_trace_0a.eps'])~=2
  print_eps(['fig_Acker' typePrefix '_' simprefix '_2D_vcos_trace_0a.eps'])
  saveas(gcf,['fig_Acker' typePrefix '_' simprefix '_2D_vcos_trace_0a.fig'])
  close
end

% plot non-gated lif posth:
figure; plot(dt:dt:plotduration/storeskip,posth(1,1:(plotduration/dt/storeskip)));
if any(posthspikes(1,1:(plotduration/dt/storeskip))>0)
  hold on; plot(dt*find(posthspikes(1,1:(plotduration/dt/storeskip))>0),1,'r.','MarkerSize',16);
end
set(gcf,'position',[82 404 3*257 2/3*302])
set(gca,'fontsize',12,'box','off')
xlabel('Time (ms)','fontsize',12); ylabel('Voltage (mV)','fontsize',12)
set(gca,'ylim',[-0.1 1.1]);
if saveFigures && exist(['fig_Acker' typePrefix '_' simprefix '_2D_posth_nogatelif_trace_0a.eps'])~=2
  print_eps(['fig_Acker' typePrefix '_' simprefix '_2D_posth_nogatelif_trace_0a.eps'])
  saveas(gcf,['fig_Acker' typePrefix '_' simprefix '_2D_posth_nogatelif_trace_0a.fig'])
end

% plot gated lif posth:
figure; plot(dt:dt:plotduration/storeskip,posth(2,1:(plotduration/dt/storeskip)));
if any(posthspikes(2,1:(plotduration/dt/storeskip))>0)
  hold on; plot(dt*find(posthspikes(2,1:(plotduration/dt/storeskip))>0),1,'r.','MarkerSize',16);
end
set(gcf,'position',[82 404 3*257 2/3*302])
set(gca,'fontsize',12,'box','off')
xlabel('Time (ms)','fontsize',12); ylabel('Voltage (mV)','fontsize',12)
set(gca,'ylim',[-0.1 1.1]);
if saveFigures && exist(['fig_Acker' typePrefix '_' simprefix '_2D_posth_gating_trace_0a.eps'])~=2
  print_eps(['fig_Acker' typePrefix '_' simprefix '_2D_posth_gating_trace_0a.eps'])
  saveas(gcf,['fig_Acker' typePrefix '_' simprefix '_2D_posth_gating_trace_0a.fig'])
end

% plot non-gated res posth:
figure; plot(dt:dt:plotduration/storeskip,real(posth(3,1:(plotduration/dt/storeskip))));
if any(posthspikes(3,:))
  hold on; plot(dt*find(posthspikes(3,1:(plotduration/dt/storeskip))>0),1,'r.','MarkerSize',16);
end
set(gcf,'position',[82 404 3*257 2/3*302])
set(gca,'fontsize',12,'box','off')
xlabel('Time (ms)','fontsize',12); ylabel('Voltage (mV)','fontsize',12)
set(gca,'ylim',[-1.1 1.1]);
if saveFigures && exist(['fig_Acker' typePrefix '_' simprefix '_2D_posth_nogateres_trace_0a.eps'])~=2
  print_eps(['fig_Acker' typePrefix '_' simprefix '_2D_posth_nogateres_trace_0a.eps'])
  saveas(gcf,['fig_Acker' typePrefix '_' simprefix '_2D_posth_nogateres_trace_0a.fig'])
end

% these in normalized coords:
trajH = 3/4*1/2; % <1/2
trajW = 3/4*1/3; % <1/3
diffH = 2/3*1/3; % <1/3
diffW = trajW; % <1/3

col1L = 1/4*1/3; % <1/3-trajW
col2L = 1/3+1/5*1/3;
col3L = 2/3+1/5*1/3;
row2B = 1/6*1/2; % <1/3-trajH
diffL = 1/5*1/3;
diffB1 = .007+1/8*1/2;
diffB2 = diffB1+0.005;
%   diffB3 = trajB;

figure('position',[82 404 3*257 2*277])
if nVCOs==2
  % gated lif:
  subplot('position',[col1L+0 row2B+1/2 trajW trajH]);
  plot(pos(1,1:500:length(pos)),pos(2,1:500:length(pos)),'.','MarkerSize',8);
  hold on; plot(pos(1,round(find(posthspikes(2,:)>0))),pos(2,round(find(posthspikes(2,:)>0))),'R.','MarkerSize',15)
  title('Gated posthsynaptic cell','fontsize',12)
  xlabel('Position (m)','fontsize',12), ylabel('Position (m)','fontsize',12)
  set(gca,'fontsize',12,'box','off')
  set(gca,'xlim',[min(pos(1,:)) max(pos(1,:))],'ylim',[min(pos(2,:)) max(pos(2,:))])

  % gated lif xcorr:
  edges{1} = linspace(min(pos(1,:)),max(pos(1,:)),50);
  edges{2} = linspace(min(pos(2,:)),max(pos(2,:)),50);
  rate = hist3([pos(1,round(find(posthspikes(2,:)>0))); pos(2,round(find(posthspikes(2,:)>0)))]','Edges',edges);
  subplot('position',[col2L+0 row2B+1/2 trajW trajH]); imagesc(rot90(conv2(rate,rate,'same')))
  if showAbstractFig
    title('Network model autocorrelogram','fontsize',12)
  else
    title(['Spatial autocorrelogram'],'fontsize',12)
  end
  set(gca,'xtick',[],'ytick',[]);
  set(gca,'fontsize',12)
  axis tight

  if showAbstractFig==1
    % abstract:
    abstractThr = 1.8;
    subplot('position',[col1L+0 row2B+0 trajW trajH]);
    plot(pos(1,1:500:length(pos)),pos(2,1:500:length(pos)),'.','MarkerSize',8);
    hold on; plot(pos(1,(find(abstractGrid>abstractThr))),pos(2,(find(abstractGrid>abstractThr))),'R.','MarkerSize',15)
    title('Abstract model','fontsize',12)
    xlabel('Position (m)','fontsize',12), ylabel('Position (m)','fontsize',12)
    set(gca,'xlim',[min(pos(1,:)) max(pos(1,:))],'ylim',[min(pos(2,:)) max(pos(2,:))])
    set(gca,'fontsize',12,'box','off')
    
    rate = hist3([pos(1,find(abstractGrid>abstractThr)); pos(2,find(abstractGrid>abstractThr))]','Edges',edges);
    subplot('position',[col2L+0 row2B trajW trajH]); imagesc(rot90(conv2(rate,rate,'same')))
    title(['Abstract model autocorrelogram'],'fontsize',12)
    set(gca,'xtick',[],'ytick',[]);
    set(gca,'fontsize',12)
    axis tight
  elseif showAbstractFig==0
    % non-gated res:
    subplot('position',[col1L+0 row2B+0 trajW trajH]);
    plot(pos(1,1:500:length(pos)),pos(2,1:500:length(pos)),'.','MarkerSize',8);
    hold on; plot(pos(1,round(find(posthspikes(3,:)>0))),pos(2,round(find(posthspikes(3,:)>0))),'R.','MarkerSize',15)
    title('Resonator posthsynaptic cell','fontsize',12)
    xlabel('Position (m)','fontsize',12), ylabel('Position (m)','fontsize',12)
    set(gca,'xlim',[min(pos(1,:)) max(pos(1,:))],'ylim',[min(pos(2,:)) max(pos(2,:))])
    set(gca,'fontsize',12,'box','off')

    rate = hist3([pos(1,round(find(posthspikes(3,:)==1))); pos(2,round(find(posthspikes(3,:)==1)))]','Edges',edges);
    subplot('position',[col2L+0 row2B trajW trajH]); imagesc(rot90(conv2(rate,rate,'same')))
    title(['Resonator autocorrelogram'],'fontsize',12)
    set(gca,'xtick',[],'ytick',[]);
    set(gca,'fontsize',12)
    axis tight
  end

  % vco 1 phase diff:
  spikeskips = 25;
  subplot('position',[col3L diffB2+2/3 diffW diffH]);
  plot(mod(vcohist(1,round(VCOSpikeTimes{1}(1:spikeskips:end)/dt)),2*pi),'.')
  %   hold on; plot(2-1e4*dxs(round(VCOSpikeTimes{1}(1:spikeskips:end)/dt)),'k.');
  title('VCO 1 phase drift','fontsize',12)
  ylabel('Phase difference (rad)','fontsize',12)
  set(gca,'xlim',[0 length(round(VCOSpikeTimes{1}(1:spikeskips:end)))],'ylim',[0 2*pi],'xtick',[])
  set(gca,'fontsize',12,'box','off')
  
  subplot('position',[diffL+2/3 diffB2+1/3 diffW diffH]);
  plot(mod(vcohist(2,round(VCOSpikeTimes{2}(1:spikeskips:end)/dt)),2*pi),'.')
  title('VCO 2 phase drift','fontsize',12)
  ylabel('Phase difference (rad)','fontsize',12)
  set(gca,'xlim',[0 length(round(VCOSpikeTimes{2}(1:spikeskips:end)))],'ylim',[0 2*pi],'xtick',[])
  set(gca,'fontsize',12,'box','off')

  % base phase diff:
  subplot('position',[col3L diffB2 diffW diffH]);
  plot(mod(basehist(round(BaseSpikeTimes(1:spikeskips:end)/dt)),2*pi),'.')
  title('Baseline phase drift','fontsize',12)
  xlabel('Time','fontsize',12)
  ylabel('Phase difference (rad)','fontsize',12)
  set(gca,'xlim',[0 length(round(BaseSpikeTimes(1:spikeskips:end)))],'ylim',[0 2*pi],'xtick',[])
  set(gca,'fontsize',12,'box','off')
end
if saveFigures && exist(['fig_Acker' typePrefix '_' simprefix '_2D_all_0a.eps'])~=2
  save(['fig_Acker' typePrefix '_' simprefix '_2D_variables_0a.mat'])
  % free some memory:
  %   clear pos abstractGrid vhist vcohist basehist posth
  print_eps(['fig_Acker' typePrefix '_' simprefix '_2D_all_0a.eps'])
  saveas(gcf,['fig_Acker' typePrefix '_' simprefix '_2D_all_0a.fig'])
  close
end
