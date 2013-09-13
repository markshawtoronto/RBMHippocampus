% Acker et al 2003 model variant neuron network with random connectivity
% with integrate and fire postsynaptic cell
% eric zilli - oct 10, 2009
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
% Interesting questions to continue this line of research: Can we calculate
% these statistics analytically? Approximately? Can we predict the
% statistics of a single cell in a network based on numerical measurements
% of a single uncoupled cell?
%
% This script runs a network of optionally-noisy cells while the injected
% current is held constant. It allows the number of cells, the coupling
% type and strength, the connectivity probability, amount of noise, etc. to
% be changed.
%
% Note: this is a cleaned up version of the script used to generate the
% paper figures--it is possible errors were introduced during the cleaning!
% Feel free to contact me if there is any difficulty reproducing any
% results from the manuscript.

clear all;
warning off

simdur = 4000; % (ms)
dt = .01; % ms

nruns = 1;

% these allow initial transients and such to be ignored:
startAnalysisTime = dt; % (ms)
endAnalysisTime = simdur; % (ms)

useNoise = 1;
useFrozenNoise = 0; % requires more memory
spikeThreshold = 20; % mV

% if true, use a slower version of the H current
slowslow = 0;

% if true, analyze only the last half of the ISIs
lastHalf = 1;

useGapJunctions = 0;
if ~slowslow
  uniqueNoiseSTD = 3.44*useNoise;
else
  uniqueNoiseSTD = 10*useNoise;
end

% approximate probability any two cells are connected (autoconnections will
% be removed)
pcons = 1;

% type = 1 n=150 uncoupled, set I and noise to match biology
% type = 2 n=250, synaptic, varying g (to find coupling for Figure 10)
% type = 3 n=250, gap-junction, varying g
for type=[3]
  if type==1
    useGapJunctions = 0;
    allncells = 150;
    pcons = 0;
    gs = 0;
    I = -2.37;
    simdur = 4000; % ms
    useNoise = 1;
    uniqueNoiseSTD = 3.44*useNoise;
  elseif type==2
    useGapJunctions = 0;
    allncells = 250;
    pcons = 1;
    gs = 2.^[3:.5:8];
    I = -2.37;
    simdur = 4000; % ms
    useNoise = 1;
    uniqueNoiseSTD = 3.44*useNoise;
  elseif type==3
    useGapJunctions = 1;
    allncells = 250;
    pcons = 1;
    gs = 2.^[-6:.5:4];
    I = -2.37;
    simdur = 4000; % ms
    useNoise = 1;
    uniqueNoiseSTD = 3.44*useNoise;
  end
end

fileOut = 0;
if ~fileOut
  'no fileout'
end
outfilename = 'acker_model_';
if ncells>1
  if useGapJunctions
    outfilename = [outfilename 'g'];
  else
    outfilename = [outfilename 's'];
  end
  endif useNoise
  outfilename = [outfilename 'n'];
end
outfilename = [outfilename sprintf('_vary_gncellspcon_%gnoise_%s.txt',useNoise,datestr(now,'mmmmdd'))];

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

mins = zeros(length(allncells),length(pcons),length(gs));
medians = zeros(length(allncells),length(pcons),length(gs));
maxs = zeros(length(allncells),length(pcons),length(gs));

for ncellsi=1:length(allncells)
  ncells = allncells(ncellsi);
  for gi=1:length(gs)
    g = gs(gi);
    for pind=1:length(pcons)
      pcon = pcons(pind);

      stabilities = zeros(1,nruns);
      ISIstabilities = zeros(1,nruns);
      for run=1:nruns
        % activity of postsynaptic I&F
        post = zeros(1, round(simdur/dt));

        % spike times of postsynaptic I&F
        postspikes = spalloc(1, round(simdur/dt), 6*simdur); % expect firing at 6 Hz

        % I&F time constant
        tau = 4.5; % (ms)
        membraneDecay = exp(-dt/tau);
        postWeight=1.2/ncells;

        % connectivity matrix
        if pcon==1
          C = 1; % will be ignored
        elseif pcon>0
          C = (rand(ncells)<pcon)>0;
          C = sparse(C.*(1-eye(ncells))); % remove autoconnections
        else
          C = spalloc(ncells,ncells,1);
        end

        state = zeros(7,round(simdur/dt));

        if useFrozenNoise
          randn('seed',3);
          frozenNoise = randn(ncells,round(simdur/dt)+1);
        end

        % magnitude of noise added to all cells on each step
        commonNoiseSTD = 0;

        t = 0; % current time in simulation
        tind = 1;

        % initial values
        v = -55*ones(1,ncells);
        mNa = zeros(1,ncells);
        hNa = zeros(1,ncells);
        mNap = zeros(1,ncells);
        nK = zeros(1,ncells);
        mHfast = zeros(1,ncells);
        mHslow = zeros(1,ncells);
        state(:,1) = [v(1); mNa(1); hNa(1); mNap(1); nK(1); mHfast(1); mHslow(1)];

        % binary array indicating which cells fire on each time steps
        spikes = sparse(ncells, simdur/dt,ncells*simdur); % sparse to save memory in big simulations

        tic
        while t<simdur
          t = t+dt; % advance to next time step
          tind = tind+1; 1+round(t/dt);

          % if our LIF neuron is spiking constantly, let's just skip this
          % one (there is no apparent synchrony, which may either reflect the network
          % itself or the LIF parameters)
          if tind==500 && sum(postspikes)>(497)
            %             postspikes = NaN;
            break
          end

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
            mtauHslow = 65.5813 + 248.0469*exp(-(-79.2190-oldv).^2/33.5178^2); % orig
          end
          mHslow = mHslow + dt*((minfHslow-mHslow)./mtauHslow);

          if useFrozenNoise
            vnoise = frozenNoise(:,tind);
          else
            vnoise = uniqueNoiseSTD*randn(1,ncells);
          end
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
          state(:,tind) = [v(1); mNa(1); hNa(1); mNap(1); nK(1); mHfast(1); mHslow(1)];

          % spikes if upward crossing spikeThreshold mV
          spikes(:,tind) = (oldv<spikeThreshold).*(v>=spikeThreshold);

          post(:,tind) = post(:,tind-1)*membraneDecay + tau*postWeight*sum(spikes(:,tind));
          % faster(?) if we don't need postspikes:
          post(post(:,tind)>1,tind) = -1e-10;
          % slower way if postspikes is needed:
          %           postspikes(:,tind) = post(:,tind)>1;
          %           post(postspikes(:,tind)>1,tind) = 0; % reset to 0

        end
        toc

        clear C

        % calculate mean and variance of ISIs for each postsynaptic cell
        npost = 1;
        ISImeans = zeros(1,npost);
        ISIstds = zeros(1,npost);
        ISIstabilities = zeros(1,nruns);


        %% calculate ISI statistics for each individual unit in the
        %% network:
        indivmeans = zeros(1,ncells);
        indivstds = zeros(1,ncells);
        indivstabilities = zeros(1,ncells);
        for ce=1:ncells % only run this once because we only save the last value in stabilities
          spikeTimes = find(spikes(ce,(startAnalysisTime/dt):(endAnalysisTime/dt))>0);
          ISIs = dt*diff(spikeTimes)/1000; % convert to seconds

          % instead of dropping, count ISIs as occuring between initial
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
          if lastHalf && length(ISIs)
            ISIs = ISIs(end/2:end);
          end

          % if we have 3 or fewer ISIs, let's just trash this cell:
          if length(ISIs)<=3
            indivmeans(ce) = NaN;
            indivstds(ce) = NaN;
            indivstabilities(ce) = NaN;
          else
            indivmeans(ce) = mean(ISIs);
            indivstds(ce) = std(ISIs);
            indivstabilities(ce) = 5*indivmeans(ce)^3/16/pi^2/(eps+indivstds(ce))^2;
          end
        end
        % drop NaNs and sort:
        indivmeans(isnan(indivmeans)) = [];
        indivstds(isnan(indivstds)) = [];
        indivstabilities(isnan(indivstabilities)) = [];
        [indivstabilities,ix] = sort(indivstabilities);
        indivmeans = indivmeans(ix);
        indivstds = indivstds(ix);

        %% calculate ISI statistics for the postsynaptic LIF neurons
        for ce=1:npost % only run this once because we only save the last value in stabilities
          %           spikeTimes = find(postspikes(ce,(startAnalysisTime/dt):(endAnalysisTime/dt))>0);
          spikeTimes = find(post(ce,(startAnalysisTime/dt):(endAnalysisTime/dt))==-1e-10);
          ISIs = dt*diff(spikeTimes)/1000; % convert to seconds

          % instead of dropping, count ISIs as occuring between initial
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

        % save to file:
        if fileOut
          fid = fopen(outfilename,'a');
          fprintf(fid,'%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\n', ...
            ncells,pcon,I,gs(gi),simdur/1000, useNoise*uniqueNoiseSTD,tau,postWeight,...
            spikeThreshold,...
            mean(1./ISIs),ISImeans(1),ISIstds(1),ISIstabilities(1),length(ISIs),...
            min(indivmeans),median(indivmeans),max(indivmeans),...
            min(indivstds),median(indivstds),max(indivstds),...
            min(indivstabilities),median(indivstabilities),max(indivstabilities),...
            length(find(postspikes>0)));
          fclose(fid);
        end
        fprintf('%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\n', ...
          ncells,pcon,I,gs(gi),simdur/1000, useNoise*uniqueNoiseSTD,tau,postWeight,...
          spikeThreshold,...
          mean(1./ISIs),ISImeans(1),ISIstds(1),ISIstabilities(1),length(ISIs), ...
          min(indivmeans),median(indivmeans),max(indivmeans),...
          min(indivstds),median(indivstds),max(indivstds),...
          min(indivstabilities),median(indivstabilities),max(indivstabilities),...
          length(find(postspikes>0)));

        if ~isempty(indivname) && fileOut
          fid = fopen(indivname,'a');
          fprintf(fid,'%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d', ...
            ncells,pcon,I,gs(gi),simdur/1000, useNoise*uniqueNoiseSTD,tau,postWeight,...
            mean(1./ISIs),ISImeans(1),ISIstds(1),ISIstabilities(1),length(ISIs),...
            length(find(postspikes>0)));
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

      medians(ncellsi,pind,gi) = median(stabilities);
      maxs(ncellsi,pind,gi) = max(stabilities);
      mins(ncellsi,pind,gi) = min(stabilities);
    end
  end
end % end ncells loop

warning on