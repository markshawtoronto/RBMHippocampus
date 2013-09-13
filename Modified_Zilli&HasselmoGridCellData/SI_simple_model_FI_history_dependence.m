% izhikevich simple model neuron network with random connectivity
% with time varying input with LIF-measured stability
% eric zilli - june 28, 2009 using simple model MATLAB code from
% Izhikevich's book modified to be a network
%
% release version 1.0. check modeldb for updates.
%
% this source code is released into the public domain
%
% This script starts a network with a specified current level (or
% frequency) and after a specified number of spikes allows the current
% level to be changed. This can be done with a single starting and ending
% current level (Figure S3A-D), a single starting and multiple ending levels
% (Figure S3E-F), multiple starting and single ending level, or a grid of injected
% current levels (Figure 4). These correspond to three types of figures that
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

% type=1 - (Figure S2A,C) class 2 excitable, n=1, noise-free, 7 Hz start, 11 Hz end
% type=2 - (Figure S2E) class 2 excitable, n=1, noise-free, 7 Hz start, 7:11 Hz end
% type=3 - class 2 excitable, n=1, noise-free, 11 Hz start, 7:11 Hz end
% type=4 - (Figure S3) class 2 excitable, n=1, noise-free, 7:0.5:11 Hz start, 7:0.5:11 Hz end
for type=[4];
  simdur = 8000; % (ms)
  dt = .1; % ms

  fileOut = 0;
  saveFigures = 0;
  
  % if journalChargesALotForColor, gray-scale plots will be produced
  journalChargesALotForColor=1;

  nruns =1;

  % number of ISIS to ignore (to ignore start-up transients)
  numDropISIs = 10;
  % number of ISIs (after numDropISIs) to use first current injection:
  numPreISIs = 10;
  % number of ISIs (after numDropISIs+numPreISIs) for second current injection:
  numPostISIs = 10;

  allncells = 1;
  pcon = 1;
  useNoise = 0;

  if type==1
    excitationClass = 2;
    useGapJunctions = 1;
    g = 0;
    uniqueNoiseSTD = 0*useNoise;
    I = 125;
    allncells = 1;
    pcon = 1;
    load simple_model_RS2_ext_FI_Jan09_n1.mat;
    inputMags = interp1(freqs,currents,7,'cubic');
    inputMags2 = interp1(freqs,currents,11,'cubic');
  elseif type==2
    excitationClass = 2;
    useGapJunctions = 1;
    g = 0;
    uniqueNoiseSTD = 0*useNoise;
    I = 125;
    allncells = 1;
    pcon = 1;
    load simple_model_RS2_ext_FI_Jan09_n1.mat;
    desiredFreqs = 7:11;
    inputMags = interp1(freqs,currents,7,'cubic');
    inputMags2 = interp1(freqs,currents,desiredFreqs,'cubic');
  elseif type==3
    excitationClass = 2;
    useGapJunctions = 1;
    g = 0;
    uniqueNoiseSTD = 0*useNoise;
    I = 125;
    allncells = 1;
    pcon = 1;
    load simple_model_RS2_ext_FI_Jan09_n1.mat;
    desiredFreqs = 7:11;
    inputMags = interp1(freqs,currents,11,'cubic');
    inputMags2 = interp1(freqs,currents,desiredFreqs,'cubic');
  elseif type==4
    excitationClass = 2;
    useGapJunctions = 1;
    g = 0;
    uniqueNoiseSTD = 0*useNoise;
    I = 125;
    allncells = 1;
    pcon = 1;
    load simple_model_RS2_ext_FI_Jan09_n1.mat;
    desiredFreqs = 7:0.5:11;
    inputMags = interp1(freqs,currents,desiredFreqs,'cubic');
    inputMags2 = interp1(freqs,currents,desiredFreqs,'cubic');
  end
   
  typePrefix = [sprintf('%d',excitationClass)];
	if allncells>1
    if useGapJunctions
      typePrefix = [typePrefix 'g'];
    else
      typePrefix = [typePrefix 's'];
    end
  end
  if useNoise
    typePrefix = [typePrefix 'n'];
  end
  outfilename = ['simple_model_RS' typePrefix '_LIF_constant_input_dec20_fi_histdep_draft1.txt'];

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

  n=round(simdur/dt); % number of simulation steps

  ncells = allncells;

  slopes = zeros(length(inputMags));
  if length(inputMags)==1 && length(inputMags2)>1
    figISIs = zeros(length(inputMags2),numPreISIs+numPreISIs);
  elseif length(inputMags)>1 && length(inputMags2)==1
    figISIs = zeros(length(inputMags),numPreISIs+numPreISIs);
  elseif length(inputMags)>1 && length(inputMags2)>1
    figISIs = zeros(length(inputMags),length(inputMags2));
  end

  % connectivity matrix
  C = (rand(ncells)<pcon)>0;
  C = sparse(C.*(1-eye(ncells))); % remove autoconnections

  for inputi=1:length(inputMags)
    for inputi2=1:length(inputMags2)

      stabilities = zeros(1,nruns);
      stabilities2 = zeros(1,nruns);
      for run=1:nruns
        % activity of postsynaptic I&F
        npost = 1;
        post = zeros(npost, round(simdur/dt)+1);

        % I&F params
        tau = 4.5; % (sec)
        membraneDecay = exp(-dt/tau);
        postWeight=4/ncells;

        % magnitude of noise added to all cells on each step
        commonNoiseSTD = 0;
        commonNoise = 0;

        numISIs = -1;

        t = 0; % current time in simulation
        tind = 1;

        v = vr*ones(ncells,1);
        u = 0*v; % initial values

        % save state of one cell over simulation:
        state = zeros(2,simdur/dt);

        % binary array indicating which cells fire on each time steps
        spikes = spalloc(ncells, simdur/dt,ncells*simdur); % sparse to save memory in big simulations

        I = inputMags(inputi);

        tic
        while t<simdur
          t = t+dt; % advance to next time step
          tind = 1+tind;

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
              v = oldv + dt*(k*(oldv-vr).*(oldv-vt) - oldu + I + commonNoise + uniqueNoiseSTD*randn(ncells,1) + g*(C*spikes))/Cf;
            end
          end
          u = oldu + dt*a*(b*(oldv-vr)-oldu);

          % save and reset spikes when v>=vpeak
          spikes(:,tind) = v>=vpeak;
          if excitationClass>0
            u(v>=vpeak) = u(v>=vpeak)+d;
          else
            u(v>=vpeak) = d;
          end
          oldv(v>=vpeak) = vpeak;
          v(v>=vpeak) = c;

          %% do numPreISIs then change I and do numPostISIs
          if any(spikes(1,tind))
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

          state(:,tind) = [oldv(1); oldu(1)];

          post(:,tind) = post(:,tind-1)*membraneDecay + tau*postWeight*sum(spikes(:,tind));
          % faster if we don't need postspikes:
          post(post(:,tind)>1,tind) = -1e-10;
        end
        toc

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
          fprintf(fid,'%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\n', ...
            ncells,pcon,inputMags(inputi),inputMags2(inputi2),g,simdur/1000, useNoise*uniqueNoiseSTD,tau,postWeight,...
            Cf,vr,vt,k,a,b,c,d,vpeak,...
            mean(ISIs(2:end-2)), ISIs(numDropISIs+numPreISIs+1), ISIs(end-1), ISIs(end),length(ISIs));
          fclose(fid);
        end
        fprintf('%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\n', ...
          ncells,pcon,inputMags(inputi),inputMags2(inputi2),g,simdur/1000, useNoise*uniqueNoiseSTD,tau,postWeight,...
          Cf,vr,vt,k,a,b,c,d,vpeak,...
          mean(ISIs(2:end-2)), ISIs(numDropISIs+numPreISIs+1), ISIs(end-1), ISIs(end),length(ISIs));

        if inputMags2(inputi2)~=inputMags(inputi)
          slopes(inputi,inputi2) = (ISIs(end-1)-ISIs(end))/(inputMags2(inputi2)-inputMags(inputi));
        else
          slopes(inputi,inputi2) = 0;
        end
      end
    end

    % average ISIs over multiple runs if noisy network is in use:
    if length(inputMags)>1 || length(inputMags2)>1
      figISIs = figISIs./nruns;
    end
  end

end

if length(inputMags)==1 && length(inputMags2)==1
  start=25000;stop=35000; % set these as desired
  figure;
  plot((start:stop)*dt,state(1,start:stop)');
  set(gca,'xlim',[start*dt stop*dt]);
  set(gca,'fontsize',12)
  set(gcf,'position',[82 404 257 185]);
  xlabel('Time (ms)','fontsize',12); ylabel('Voltage (mV)','fontsize',12);
  if saveFigures && exist(['fig_RS' typePrefix '_histdep_vtrace.eps'])~=2
    print_eps(['fig_RS' typePrefix '_histdep_vtrace.eps'])
    saveas(gcf,['fig_RS' typePrefix '_histdep_vtrace.fig'])
  end
  figure;
  plot((start:stop)*dt,state(2,start:stop)');
  set(gcf,'position',[82 404 257 227]);
  set(gca,'ygrid','on','xlim',[start*dt stop*dt]);
  set(gca,'fontsize',12)
  xlabel('Time (ms)','fontsize',12); ylabel('Recovery variable','fontsize',12);
  if saveFigures && exist(['fig_RS' typePrefix '_histdep_utrace.eps'])~=2
    print_eps(['fig_RS' typePrefix '_histdep_utrace.eps'])
    saveas(gcf,['fig_RS' typePrefix '_histdep_utrace.fig'])
  end
elseif (length(inputMags)==1 && length(inputMags2)>1) || (length(inputMags)>1 && length(inputMags2)==1)
  figure; plot(1./figISIs','b.-','Markersize',16);
  xlabel('Interspike interval #','fontsize',12); ylabel('1/ISI duration (Hz)','fontsize',12)
  set(gcf,'position',[82 404 257 302])
  set(gca,'ylim',[min(desiredFreqs)-0.1 max(desiredFreqs)+0.5]);
  set(gca,'fontsize',12)
  if saveFigures && exist(['fig_RS' typePrefix '_histdep_multi.eps'])~=2
    print_eps(['fig_RS' typePrefix '_histdep_multi.eps'])
    saveas(gcf,['fig_RS' typePrefix '_histdep_multi.fig'])
  end
elseif length(inputMags)>1 && length(inputMags2)>1
  figure; imagesc(desiredFreqs,desiredFreqs,1000*abs(real(figISIs)-imag(figISIs)));
  title('ISI overshoot magnitude (ms)','fontsize',12)
  xlabel('Starting frequency (Hz)','fontsize',12); ylabel('Ending frequency (Hz)','fontsize',12)
  set(gcf,'position',[82 404 322 302])
  set(gca,'fontsize',12)
  if journalChargesALotForColor
    colormap(flipud(gray))
  end
  h = colorbar; set(h,'fontsize',12);
  if saveFigures && exist(['fig_RS' typePrefix '_histdep_ISI_colorplot.eps'])~=2
    print_eps(['fig_RS' typePrefix '_histdep_ISI_colorplot.eps'])
    saveas(gcf,['fig_RS' typePrefix '_histdep_ISI_colorplot.fig'])
  end

  figure; imagesc(desiredFreqs,desiredFreqs,abs(1./real(figISIs)-1./imag(figISIs)));
  title('Frequency overshoot magnitude (Hz)','fontsize',12)
  xlabel('Starting frequency (Hz)','fontsize',12); ylabel('Ending frequency (Hz)','fontsize',12)
  set(gcf,'position',[82 404 322 302])
  h = colorbar; set(h,'fontsize',12);
  set(gca,'fontsize',12)
  if journalChargesALotForColor
    colormap(flipud(gray))
  end
  if saveFigures && exist(['fig_RS' typePrefix '_histdep_freq_colorplot.eps'])~=2
    print_eps(['fig_RS' typePrefix '_histdep_freq_colorplot.eps'])
    saveas(gcf,['fig_RS' typePrefix '_histdep_freq_colorplot.fig'])
  end
end
