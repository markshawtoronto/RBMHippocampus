% Acker et al 2003 model variant neuron network with random connectivity
% with time varying input with LIF-measured stability
% eric zilli - oct ??, 2009
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
% Interesting questions to continue this research: Setting slowslow=1
% decreases the adaptation's initial magnitude but the effect
% has an effect lasting many ISIs long. Is
% the corresponding time constant model realistic? Experiments looking at H
% current time constants at spiking voltages do not appear to exist, but
% gating variable "friction" in the extended HH model accounts for a
% possible non-zero base for time constants which might play a major role
% here. In-vitro experiments are needed.
% Under what conditions does the ISI period or frequency overshoot and when does it undershoot?
% What is the best way to quantify the history dependence on firing rate
% effect?
% Is it possible for the intrinsic properties of one neuron A to be exactly
% the inverse of B in the sense of this adaptation such that an
% input signal I to A would be distorted by the adaptation into an
% output I' which would then be un-distorted by B to recover a constant function
% of I? (E.g. I goes from driving at 7 Hz to 11 Hz, A's firing rate adapts from
% 11.5 Hz down to 11 Hz. Can there be a B receiving A's spikes that fires at 11 Hz
% given input at 11.5 Hz and which adapts at just the right rate that B continues to
% fire at 11 Hz as the input from A slows down slightly?)
%
% This script starts a network with a specified current level (or
% frequency) and after a specified number of spikes allows the current
% level to be changed. This can be done with a single starting and ending
% current level (Figure S2A-D), a single starting and multiple ending levels
% (Figure S2E-F), multiple starting and single ending level, or a grid of injected
% current levels (Figure S3). These correspond to three types of figures that
% might be of interest.
%
% If the FI curve of the cell or network is known, it can be used to find
% the appropriate current levels. Otherwise they must be set by hand. The
% FI curve needn't be terribly high in resolution, so generating them as
% needed should be quick and easy.
%
% Note: this is a cleaned up version of the script used to generate the
% paper figures--it is possible errors were introduced during the cleaning!
% Feel free to contact me if there is any difficulty reproducing any
% results from the manuscript.
%

clear all;

% type=1 - (Figure S2B,D) class 2 excitable, n=1, noise-free, 1 start, 1 end
% type=2 - (Figure S2F) class 2 excitable, n=1, noise-free, 1 start, 5 end, low-to-high
% type=3 - class 2 excitable, n=1, noise-free, 7:0.5:11 Hz start, 7:0.5:11 Hz end
for type=[3]
  simdur = 10000; % (ms)
  dt = .01; % ms

  nruns =1;

  saveFigures = 0;

  useNoise = 0;
  spikeThreshold = 20; % mV
  slowslow = 0;

  numDropISIs = 10;
  numPreISIs = 10;
  numPostISIs = 10;

  if type==1
    excitationClass = 2;
    useGapJunctions = 1;
    g = 0;
    uniqueNoiseSTD = 0*useNoise;
    ncells = 1;
    pcon = 0;
    load Acker_FI_n1.mat;
    inputMags = interp1(freqs,currents,7,'cubic');
    inputMags2 = interp1(freqs,currents,11,'cubic');
  elseif type==2
    excitationClass = 2;
    useGapJunctions = 1;
    g = 0;
    uniqueNoiseSTD = 0*useNoise;
    ncells = 1;
    pcon = 0;
    load Acker_FI_n1.mat;
    desiredFreqs = 7:11;
    inputMags = interp1(freqs,currents,7,'cubic');
    inputMags2 = interp1(freqs,currents,desiredFreqs,'cubic');
  elseif type==3
    excitationClass = 2;
    useGapJunctions = 1;
    g = 0;
    uniqueNoiseSTD = 0*useNoise;
    ncells = 1;
    pcon = 0;
    load Acker_FI_n1.mat;
    desiredFreqs = [7 11];
    inputMags = interp1(freqs,currents,desiredFreqs,'cubic');
    inputMags2 = interp1(freqs,currents,desiredFreqs,'cubic');
  end

  fileOut = 0;
  typePrefix = [sprintf('%d',excitationClass)];
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
  outfilename = sprintf('acker_model_%s_history_dependence_%dn_%gnoise_%s.txt',typePrefix,ncells,useNoise,datestr(now,'mmmmdd'));

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

  slopes = zeros(length(inputMags));
  if length(inputMags)==1 && length(inputMags2)>1
    figISIs = zeros(length(inputMags2),numPreISIs+numPreISIs);
  elseif length(inputMags)>1 && length(inputMags2)==1
    figISIs = zeros(length(inputMags),numPreISIs+numPreISIs);
  elseif length(inputMags)>1 && length(inputMags2)>1
    figISIs = zeros(length(inputMags),length(inputMags2));
  end

  % connectivity matrix
  if pcon==1
    C = 1;
  else
    C = (rand(ncells)<pcon)>0;
    C = sparse(C.*(1-eye(ncells))); % remove autoconnections
  end

  for inputi=1:length(inputMags)
    for inputi2=1:length(inputMags2)

      stabilities = zeros(1,nruns);
      stabilities2 = zeros(1,nruns);
      for run=1:nruns
        % activity of postsynaptic I&F
        npost = 1;
        post = zeros(npost, round(simdur/dt)+1);

        % spike times of postsynaptic I&F
        postspikes = 0;

        % I&F params
        tau = 2.5; % (sec)
        membraneDecay = exp(-dt/tau);
        postWeight=4/ncells;

        % magnitude of noise added to all cells on each step
        commonNoiseSTD = 0;

        numISIs = -1;

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

        % save state of one cell over simulation:
        state = zeros(2,simdur/dt);

        % binary array indicating which cells fire on each time steps
        spikes = spalloc(ncells, simdur/dt,ncells*simdur); % sparse to save memory in big simulations

        I = inputMags(inputi);
        Is = zeros(1,round(simdur/dt));

        tic
        while t<simdur
          t = t+dt; % advance to next time step
          tind = 1+tind; %round(t/dt);

          %% done below now:
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

          vnoise = uniqueNoiseSTD*randn(1,ncells);
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
          state(:,tind) = [v(1) mHslow];

          % spikes if upward crossing spikeThreshold mV
          spikes(:,tind) = (oldv<spikeThreshold).*(v>=spikeThreshold);

          %% do numPreISIs then change I and do numPostISIs
          if any(spikes(:,tind))
            numISIs = numISIs + 1; % starts equal to -1 so second spike is first ISI
          end
          if numISIs<(numPreISIs+numDropISIs)
            I = inputMags(inputi);
          else
            I = inputMags2(inputi2);
          end
          if numISIs==(numPreISIs+numPostISIs+numDropISIs)
            break
          end

          post(:,tind) = post(:,tind-1)*membraneDecay + tau*postWeight*sum(spikes(:,tind));
          % faster(?) if we don't need postspikes:
          post(post(:,tind)>1,tind) = -1e-10;
          % slower way if postspikes is needed:
          %       postspikes(:,tind) = post(:,tind)>1;
          %       post(postspikes(:,tind)>1,tind) = 0; % reset to 0

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
        ISIs = dt*diff(find(spikes(1,:)>0))/1000;

        if length(inputMags)==1 && length(inputMags2)>1 % add so we can take the average over runs
          figISIs(inputi2,:) = figISIs(inputi2,:) + ISIs((numDropISIs+1):end);
        elseif length(inputMags)>1 && length(inputMags2)==1
          figISIs(inputi,:) = figISIs(inputi,:) + ISIs((numDropISIs+1):end);
        elseif length(inputMags)>1 && length(inputMags2)>1
          % using a complex number to hold two values in one matrix:
          figISIs(inputi2,inputi) = figISIs(inputi2,inputi) + ISIs(numDropISIs+numPreISIs+1) + ISIs(end-1)*i;
        end

        % save to file:
        if fileOut
          fid = fopen(outfilename,'a');
          fprintf(fid,'%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\n', ...
            ncells,pcon,inputMags(inputi),inputMags2(inputi2),g,simdur/1000, useNoise*uniqueNoiseSTD,tau,postWeight,...
            mean(ISIs(2:end-numPostISIs-2)), ISIs(end-numPostISIs-2), ISIs(end-numPostISIs+1), ISIs(end-numPostISIs+2),ISIs(end),length(ISIs),...
            (inputMags2(inputi2)-inputMags(inputi)), ISIs(end-numPostISIs+1)-ISIs(end-numPostISIs+2),ISIs(end-numPostISIs+1)-ISIs(end));
          fclose(fid);
        end
        fprintf('%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\n', ...
          ncells,pcon,inputMags(inputi),inputMags2(inputi2),g,simdur/1000, useNoise*uniqueNoiseSTD,tau,postWeight,...
          mean(ISIs(2:end-numPostISIs-2)), ISIs(end-numPostISIs-2), ISIs(end-numPostISIs+1), ISIs(end-numPostISIs+2),ISIs(end),length(ISIs),...
          (inputMags2(inputi2)-inputMags(inputi)), ISIs(end-numPostISIs+1)-ISIs(end-numPostISIs+2),ISIs(end-numPostISIs+1)-ISIs(end));

        slopes(inputi,inputi2) = (ISIs(end-numPostISIs-1)-ISIs(end-numPostISIs))/(inputMags2(inputi2)-inputMags(inputi));
      end
    end

    % average ISIs over multiple runs if noisy network is in use:
    if length(inputMags)>1 || length(inputMags2)>1
      figISIs = figISIs./nruns;
    end
  end
end

state = state(:,1:tind);

if length(inputMags)==1 && length(inputMags2)==1
  start=260000;stop=360000;
  figure;
  plot((start:stop)*dt,state(1,start:stop)');
  set(gca,'xlim',[start*dt stop*dt]);
  set(gcf,'position',[82 404 257 185]);
  xlabel('Time (ms)','fontsize',12); ylabel('Voltage (mV)','fontsize',12);
  set(gca,'fontsize',12,'box','off')
  if slowslow
    if saveFigures && exist(['fig_Ackerslow_history_dependence_vtrace_0a.eps'])~=2
      print_eps(['fig_Ackerslow_history_dependence_vtrace_0a.eps'])
      saveas(gcf,['fig_Ackerslow_history_dependence_vtrace_0a.fig'])
    end
  else
    if saveFigures && exist(['fig_Acker_history_dependence_vtrace_0a.eps'])~=2
      print_eps(['fig_Acker_history_dependence_vtrace_0a.eps'])
      saveas(gcf,['fig_Acker_history_dependence_vtrace_0a.fig'])
    end
  end
  figure;
  plot((start:stop)*dt,state(2,start:stop)');
  set(gcf,'position',[82 404 257 227]);
  set(gca,'ygrid','on','xlim',[start*dt stop*dt]);
  xlabel('Time (ms)','fontsize',12); ylabel('Slow H','fontsize',12);
  set(gca,'fontsize',12,'box','off')
  if slowslow
    if saveFigures && exist(['fig_Ackerslow_history_dependence_htrace_0a.eps'])~=2
      print_eps(['fig_Ackerslow_history_dependence_htrace_0a.eps'])
      saveas(gcf,['fig_Ackerslow_history_dependence_htrace_0a.fig'])
    end
  else
    if saveFigures && exist(['fig_Acker_history_dependence_htrace_0a.eps'])~=2
      print_eps(['fig_Acker_history_dependence_htrace_0a.eps'])
      saveas(gcf,['fig_Acker_history_dependence_htrace_0a.fig'])
    end
  end
elseif (length(inputMags)==1 && length(inputMags2)>1) || (length(inputMags)>1 && length(inputMags2)==1)
  figure; plot(1./figISIs','b.-','Markersize',16);
  xlabel('Interspike interval #','fontsize',12); ylabel('1/ISI duration (Hz)','fontsize',12)
  set(gcf,'position',[82 404 257 302])
  set(gca,'ylim',[min(desiredFreqs)-0.1 max(desiredFreqs)+0.5]);
  set(gca,'fontsize',12,'box','off')
  if slowslow
    if saveFigures && exist(['fig_Ackerslow_history_dependence_multi_0a.eps'])~=2
      print_eps(['fig_Ackerslow_history_dependence_multi_0a.eps'])
      saveas(gcf,['fig_Ackerslow_history_dependence_multi_0a.fig'])
    end
  else
    if saveFigures && exist(['fig_Acker_history_dependence_multi_0a.eps'])~=2
      print_eps(['fig_Acker_history_dependence_multi_0a.eps'])
      saveas(gcf,['fig_Acker_history_dependence_multi_0a.fig'])
    end
  end
elseif length(inputMags)>1 && length(inputMags2)>1
  figure; imagesc(desiredFreqs,desiredFreqs,abs(real(figISIs)-imag(figISIs)));
  title('ISI overshoot magnitude (ms)','fontsize',12)
  xlabel('Starting frequency (Hz)','fontsize',12); ylabel('Ending frequency (Hz)','fontsize',12)
  set(gcf,'position',[82 404 322 302])
  set(gca,'fontsize',12,'box','off')
  if saveFigures && exist(['fig_RS' typePrefix '_histdep_ISI_colorplot.eps'])~=2
    print_eps(['fig_RS' typePrefix '_histdep_ISI_colorplot.eps'])
    saveas(gcf,['fig_RS' typePrefix '_histdep_ISI_colorplot.fig'])
  end

  figure; imagesc(desiredFreqs,desiredFreqs,abs(1./real(figISIs)-1./imag(figISIs)));
  title('Frequency overshoot magnitude (Hz)','fontsize',12)
  xlabel('Starting frequency (Hz)','fontsize',12); ylabel('Ending frequency (Hz)','fontsize',12)
  set(gcf,'position',[82 404 322 302])
  set(gca,'fontsize',12,'box','off')
  if saveFigures && exist(['fig_RS' typePrefix '_histdep_freq_colorplot.eps'])~=2
    print_eps(['fig_RS' typePrefix '_histdep_freq_colorplot.eps'])
    saveas(gcf,['fig_RS' typePrefix '_histdep_freq_colorplot.fig'])
  end
end
