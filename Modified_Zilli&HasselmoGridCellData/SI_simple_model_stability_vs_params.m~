% izhikevich simple model neuron network with random connectivity
% with integrate and fire postsynaptic cell
% eric zilli - june 27, 2009 using simple model MATLAB code from
% Izhikevich's book modified to be a network
%
% release version 1.0. check modeldb for updates.
%
% this source code is released into the public domain
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

clear all
warning off

% type=1 - class 1 excitable, n=100, synaptic, varying g
% type=2 - class 1 excitable, n=100, gap-junction, varying g
% type=3 - (Fig S4A-C) class 2 excitable, n=250, synaptic, varying g
% type=4 - (Fig S4D-F) class 2 excitable, n=250, gap-junction, varying g
% type=5 - (Fig 5A,D) class 2 excitable, n = 1, ISI histograms
% type=6 - (Fig 5B,E) class 2 excitable, n = 250, gap-junction, ISI histograms
% type=7 - (Fig 5C,F) class 2 excitable, n = 250, synaptic, ISI histograms
% type=8 - (Fig S8) class 2 excitable, varying n and noise, gap-junction g=0.1
% type=9 - (Discussion) class 2 excitable, n=250, gap-junctions, 100% indep. noise, 0 corr. noise
% type=10 - (Discussion) class 2 excitable, n=250, gap-junctions, 95% indep. noise var, 5% corr. noise var
% type=11 - (Discussion) class 2 excitable, n=1, 5% realistic noise var
% type=12 - (Discussion) class 2 excitable, n=2500, gap-junctions, 5% realistic corr. noise var
% type=13 - class 2 excitable, n=250, inhibitory synaptic, varying g
% type=14 - (Fig 8D,H) class 2 excitable, n = 5000, p=0.01, synaptic, ISI histograms
for type=[1]
  % Note: some parameters may be overwritten in the 'if type' block:
  simdur = 8000; % (ms)
  dt = .1; % ms

  nruns = 1;

  % for statistics about each run
  fileOut = 1;
    
  % for ISI histograms:
  saveFigures = 0;
  numISIHists = 0;

  % these allow initial transients and such to be ignored:
  startAnalysisTime = dt; % (ms)
  endAnalysisTime = simdur; % (ms)

  % sets whether or not noise is used
  useNoise = 0;

  % approximate probability any two cells are connected (autoconnections will
  % be removed)
  % set to 0 (or set gs=0) to uncouple cells to, e.g., find noise level and
  % input level such that the median matches the biological data
  pcons = 1;

  weightmult = 1;
  
  % std. dev. of per-step noise term correlated among all network cells:
  commonNoiseSTD = 0*useNoise;
  commonNoise = 0;

  if type==1
    excitationClass = 1;
    useGapJunctions = 0;
    uniqueNoiseSTD = 130*useNoise;
    I = 62;
    % rough range of good gs = 150:1180
    gs = 150:50:1300;
    allncells = 1;
    nruns = 1;
  elseif type==2
    excitationClass = 1;
    useGapJunctions = 1;
    uniqueNoiseSTD = 130*useNoise;
    I = 62;
    % rough range of good gs = 0.022:16
    gs = 2.^[-5.5:0.5:4];
    allncells = 100;
    nruns = 10;
  elseif type==3
    excitationClass = 2;
    useGapJunctions = 0;
    uniqueNoiseSTD = 100*useNoise;
    I = 95.8;
    gs = 0:25:450;
    allncells = 250;
    nruns = 10;
  elseif type==4
    excitationClass = 2;
    useGapJunctions = 1;
    uniqueNoiseSTD = 100*useNoise;
    I = 95.8;
    gs = 2.^[-6:0.5:4];
    allncells = 250;
    nruns = 10;
  elseif type==5
    excitationClass = 2;
    useGapJunctions = 0;
    uniqueNoiseSTD = 100*useNoise;
    I = 95.8;
    gs = 0;
    allncells = 1;
    nruns = 1;
    simdur = 60000; % ms
    numISIHists = 1;
  elseif type==6
    excitationClass = 2;
    useGapJunctions = 1;
    uniqueNoiseSTD = 100*useNoise;
    I = 95.8;
    gs = 0.1;
    allncells = 250;
    nruns = 1;
    simdur = 60000; % ms
    numISIHists = 1;
  elseif type==7
    excitationClass = 2;
    useGapJunctions = 0;
    uniqueNoiseSTD = 100*useNoise;
    I = 95.8;
    gs = 256;
    allncells = 250;
    nruns = 1;
    simdur = 60000; % ms
    numISIHists = 1;
  elseif type==8
    excitationClass = 2;
    useGapJunctions = 1;
    %% run one at a time:
    uniqueNoiseSTD = 400*useNoise;
%     uniqueNoiseSTD = 200*useNoise;
%     uniqueNoiseSTD = 100*useNoise; % original level
%     uniqueNoiseSTD = 50*useNoise;
%     uniqueNoiseSTD = 25*useNoise;
%     uniqueNoiseSTD = 12.5*useNoise;
    I = 125;
    allncells = [16 32 64 128 256 512];
    gs = 0.1;
  elseif type==9
    excitationClass = 2;
    useGapJunctions = 1;
    uniqueNoiseSTD = 100*useNoise;
    commonNoiseSTD = 0*useNoise;
    I = 95.8;
    gs = 0.1;
    allncells = 250;
  elseif type==10
    excitationClass = 2;
    useGapJunctions = 1;
    uniqueNoiseSTD = sqrt(0.95)*100*useNoise;
    commonNoiseSTD = sqrt(0.05)*100*useNoise;
    I = 95.8;
    gs = 0.1;
    allncells = 250;
  elseif type==11
    excitationClass = 2;
    useGapJunctions = 1;
    uniqueNoiseSTD = 0*useNoise;
    commonNoiseSTD = sqrt(0.05)*100*useNoise;
    I = 95.8;
    gs = 0;
    allncells = 1;
  elseif type==12
    excitationClass = 2;
    useGapJunctions = 1;
    uniqueNoiseSTD = 0*useNoise;
    commonNoiseSTD = sqrt(0.05)*100*useNoise;
    I = 95.8;
    gs = 0.1;
    allncells = 2500;
  elseif type==13
    excitationClass = 2;
    useGapJunctions = 0;
    useNoise = 1;
    uniqueNoiseSTD = 100*useNoise;
    I = 110;
    gs = -(2.^[9:15]); % inhibitory synapses
    allncells = 1000;
    nruns = 1;
    simdur = 15000;
    weightmult = 8;
    endAnalysisTime=simdur;
elseif type==14
    excitationClass = 2;
    useGapJunctions = 0;
    uniqueNoiseSTD = 100*useNoise;
    I = 95.8;
    gs = 256;
    allncells = 5000;
    nruns = 1;
    simdur = 6000;0; % ms
    pcon = 0.01;
    numISIHists = 1;
end

  % file to save statistics for every individual cell in network, don't
  % save if indivname = ''
  indivname = '';
  % indivname = [strtok(outfilename,'.') '_individualcells.txt'];

  Cf=100; vr=-60; vt=-40; k=0.7;  % parameters used for RS
  a=0.03; c=-50; d=100;           % neocortical pyramidal neurons
  if excitationClass==1
    b=-2; % integrator
  elseif excitationClass==2
    b = 2; % resonator
  end
  vpeak=35; % spike cutoff

  % store min/med/max stability time of individual cells in each simulation
  mins = zeros(length(allncells),length(pcons),length(gs));
  medians = zeros(length(allncells),length(pcons),length(gs));
  maxs = zeros(length(allncells),length(pcons),length(gs));

  for ncellsi=1:length(allncells)
    ncells = allncells(ncellsi);
    
    typeString = '';
    if ncells>1
      if useGapJunctions
        typeString = [typeString 'g'];
      else
        typeString = [typeString 's'];
      end
    end
    if useNoise
      typeString = [typeString 'n'];
    end
    simprefix = sprintf('simple_model_RS%s_n%d_%s',typeString,ncells,datestr(now,'mmmdd'))

    if ~fileOut
      disp('no fileout')
    end
    outfilename = ['simple_model_RS' sprintf('%d',excitationClass) typeString];
    outfilename = [outfilename sprintf('_LIF_%s_vary_gandnoise_%gnoise.txt',datestr(now,'mmmdd'),useNoise)];

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

          % LIF time constant
          tau = 4.5; % (ms)
          membraneDecay = exp(-dt/tau);
          postWeight=weightmult*3/ncells;

          % connectivity matrix
          if pcon>0
            C = (rand(ncells)<pcon)>0;
            C = sparse(C.*(1-eye(ncells))); % remove autoconnections
          else
            C = spalloc(ncells,ncells,1);
          end

          vhist = zeros(1,round(simdur/dt));

          t = 0; % current time in simulation

          v = vr*ones(ncells,1); % vr*rand(ncells,1);
          u = 0*v; % initial values

          % save state of one cell over simulation:
          state = zeros(2,simdur/dt);

          % binary array indicating which cells fire on each time steps
          spikes = sparse(ncells, simdur/dt,ncells*simdur); % sparse to save memory in big simulations

          tic
          while t<simdur
            t = t+dt; % advance to next time step
            tind = 1+round(t/dt);

            if commonNoiseSTD
              commonNoise = commonNoiseSTD*randn;
            end

            % if our LIF neuron is spiking constantly, let's just skip this
            % one (there is no apparent synchrony, which may either reflect the network
            % itself or the LIF parameters)
            if tind==500 && sum(postspikes)>(497)
              %             postspikes = NaN;
              break
            end

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
                v = oldv + dt*(k*(oldv-vr).*(oldv-vt) - oldu + I + commonNoise + uniqueNoiseSTD*randn(ncells,1) + g*(sum(spikes(:,tind-1),1)-spikes(:,tind-1)))/Cf; 
              else
                v = oldv + dt*(k*(oldv-vr).*(oldv-vt) - oldu + I + commonNoise + uniqueNoiseSTD*randn(ncells,1) + g*(C*spikes(:,tind-1)))/Cf;
              end
            end
            u = oldu + dt*a*(b*(oldv-vr)-oldu);

            vhist(tind) = v(1);

            % save and reset spikes when v>=vpeak
            spikes(:,tind) = v>=vpeak;
            u(v>=vpeak) = u(v>=vpeak)+d;
            oldv(v>=vpeak) = vpeak;
            v(v>=vpeak) = c;

            state(:,tind) = [oldv(1); oldu(1)];

            post(:,tind) = post(:,tind-1)*membraneDecay + tau*postWeight*sum(spikes(:,tind));
            postspikes(:,tind) = post(:,tind)>1;
            post(postspikes(:,tind)>1,tind) = 0; % reset to 0
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

            if ce<=numISIHists
              figure; hist(ISIs,linspace(.1,.3,200))
              title(['VCO Cell #'  num2str(ce) ' ISIs'],'fontsize',12)
              xlabel('ISI duration (s)','fontsize',12)
              set(gca,'fontsize',12);
              set(gca,'box','off');
              set(gcf,'position',[82 404 257 185]);
              fname = [simprefix '_vco_isihist_' num2str(ce)];
              if saveFigures && exist([fname '.eps'])~=2
                export_fig([fname '.eps'])
                saveas(gcf,[fname '.fig'])
              end
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


          %% calculate ISI statistics for the postsynaptic LIF neurons
          for ce=1:npost % only run this once because we only save the last value in stabilities
            spikeTimes = find(postspikes(ce,(startAnalysisTime/dt):(endAnalysisTime/dt))>0);
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

            if ce==1 && numISIHists
              figure; hist(ISIs,linspace(.1,.3,200))
              title('Postsynaptic LIF ISIs','fontsize',12)
              xlabel('ISI duration (s)','fontsize',12)
              set(gca,'fontsize',12)
              set(gca,'box','off')
              set(gcf,'position',[82 404 257 185]);
              fname = [simprefix '_post_isihist'];
              if saveFigures && exist([fname '.eps'])~=2
                export_fig([fname '.eps'])
                saveas(gcf,[fname '.fig'])
              end
            end

            ISImeans(ce) = mean(ISIs);
            ISIstds(ce) = std(ISIs);
            ISIstabilities(ce) = 5*ISImeans(ce)^3/16/pi^2/(eps+ISIstds(ce))^2;
          end

          % save to file:
          if fileOut
            fid = fopen(outfilename,'a');
            fprintf(fid,'%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\n', ...
              ncells,pcon,I,gs(gi),simdur/1000, useNoise*uniqueNoiseSTD,tau,postWeight,...
              Cf,vr,vt,k,a,b,c,d,vpeak,...
              1/ISImeans(1),ISImeans(1),ISIstds(1),ISIstabilities(1),length(ISIs),...
              min(indivmeans),median(indivmeans),max(indivmeans),...
              min(indivstds),median(indivstds),max(indivstds),...
              min(indivstabilities),median(indivstabilities),max(indivstabilities),...
              length(find(postspikes>0)));
            fclose(fid);
          end
          % print to screen:
          fprintf('%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\n', ...
            ncells,pcon,I,gs(gi),simdur/1000, useNoise*uniqueNoiseSTD,tau,postWeight,...
            Cf,vr,vt,k,a,b,c,d,vpeak,...
            1/ISImeans(1),ISImeans(1),ISIstds(1),ISIstabilities(1),length(ISIs), ...
            min(indivmeans),median(indivmeans),max(indivmeans),...
            min(indivstds),median(indivstds),max(indivstds),...
            min(indivstabilities),median(indivstabilities),max(indivstabilities),...
            length(find(postspikes>0)));

          % optinally save individual cell statistics to a file:
          if ~isempty(indivname) && fileOut
            fid = fopen(indivname,'a');
            fprintf(fid,'%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d', ...
              ncells,pcon,I,gs(gi),simdur/1000, useNoise*uniqueNoiseSTD,tau,postWeight,...
              1/ISImeans(1),ISImeans(1),ISIstds(1),ISIstabilities(1),length(ISIs),...
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
end % end run loop

% warning on
% return

%% Figure 1
data = textread('simple_model_RS1_LIF_Nov15_vary_gandnoise_0noise.txt');
figure; plot(data(:,4),data(:,18),'.');
set(gca,'xlim',[0 500])
set(gca,'ylim',[0 13])
set(gca,'fontsize',12,'box','off')
xlabel('Synaptic conductance','fontsize',12); ylabel('Network frequency (Hz)','fontsize',12)
set(gcf,'position',[82 404 257 302])
if exist('fig_freq_varyg.eps')~=2
  export_fig('fig_freq_varyg.eps')
  saveas(gcf,'fig_freq_varyg.fig')
end
figure; plot(data(:,4),data(:,21),'.'); xlabel('Synaptic conductance','fontsize',12); ylabel('Estimated stability time (s)','fontsize',12)
set(gcf,'position',[82 404 257 302])
set(gca,'xlim',[0 500])
set(gca,'ylim',[0 350])
set(gca,'fontsize',12,'box','off')
if exist('fig_stability_varyg.eps')~=2
  export_fig('fig_stability_varyg.eps')
  saveas(gcf,'fig_stability_varyg.fig')
end
figure; plot(data(:,4),data(:,20),'.'); xlabel('Synaptic conductance','fontsize',12); ylabel('Network period std. dev. (s)','fontsize',12)
set(gcf,'position',[82 404 257 302])
set(gca,'xlim',[0 500])
set(gca,'yscale','log')
set(gca,'fontsize',12,'box','off')
if exist('fig_stdev_varyg.eps')~=2
  export_fig('fig_stdev_varyg.eps')
  saveas(gcf,'fig_stdev_varyg.fig')
end
data = textread('simple_model_RS1_LIF_Nov15_vary_gandnoise_0noise.txt');
figure; plot(data(:,4),data(:,18),'.');
set(gca,'xlim',[1e-2 10])
set(gca,'xtick',[.01 .1 1 10])
set(gca,'xscale','log')
set(gca,'fontsize',12,'box','off')
% set(gca,'ylim',[0 9.1])
xlabel('Gap-junction conductance','fontsize',12); ylabel('Network frequency (Hz)','fontsize',12)
set(gcf,'position',[82 404 257 302])
if exist('fig_freq_varyg.eps')~=2
  export_fig('fig_freq_varyg.eps')
  saveas(gcf,'fig_freq_varyg.fig')
end
figure; plot(data(:,4),data(:,21),'.'); xlabel('Gap-junction conductance','fontsize',12); ylabel('Estimated stability time (s)','fontsize',12)
set(gcf,'position',[82 404 257 302])
set(gca,'xlim',[1e-2 20])
set(gca,'xtick',[.01 .1 1 10])
set(gca,'xscale','log')
set(gca,'fontsize',12,'box','off')
if exist('fig_stability_varyg.eps')~=2
  export_fig('fig_stability_varyg.eps')
  saveas(gcf,'fig_stability_varyg.fig')
end
figure; plot(data(:,4),data(:,20),'.'); xlabel('Gap-junction conductance','fontsize',12); ylabel('Network period std. dev. (s)','fontsize',12)
set(gcf,'position',[82 404 257 302])
set(gca,'xlim',[1e-2 20])
set(gca,'xtick',[.01 .1 1 10])
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'fontsize',12,'box','off')
if exist('fig_stdev_varyg.eps')~=2
  export_fig('fig_stdev_varyg.eps')
  saveas(gcf,'fig_stdev_varyg.fig')
end


%% Figure 2
% histogram plots are buried in the script above, search for numISIHists


%% Figure 3
% first run type=8 for each standard deviation of the noise, saving output
% we expect all files have the same number of rows
nfiles = 6; % i.e. number of noise levels
startrow = 0;
endrow = 6;

ISIfreqs = zeros(endrow-startrow,nfiles);
ISImeans = zeros(endrow-startrow,nfiles);
ISIstds = zeros(endrow-startrow,nfiles);
ISIstability = zeros(endrow-startrow,nfiles);

data = textread('simple_model_RS1_LIF_Nov15_vary_gandnoise_0noise.txt');
ISIfreqs(:,1) = data((startrow+1):endrow,18);
ISImeans(:,1) = data((startrow+1):endrow,19);
ISIstds(:,1) = data((startrow+1):endrow,20);
ISIstability(:,1) = data((startrow+1):endrow,21);
data = textread('simple_model_RS1_LIF_Nov15_vary_gandnoise_0noise.txt');
ISIfreqs(:,2) = data((startrow+1):endrow,18);
ISImeans(:,2) = data((startrow+1):endrow,19);
ISIstds(:,2) = data((startrow+1):endrow,20);
ISIstability(:,2) = data((startrow+1):endrow,21);
data = textread('simple_model_RS1_LIF_Nov15_vary_gandnoise_0noise.txt');
ISIfreqs(:,3) = data((startrow+1):endrow,18);
ISImeans(:,3) = data((startrow+1):endrow,19);
ISIstds(:,3) = data((startrow+1):endrow,20);
ISIstability(:,3) = data((startrow+1):endrow,21);
data = textread('simple_model_RS2gn_LIF_vary_gandnoise_0.5noise.txt');
ISIfreqs(:,4) = data((startrow+1):endrow,18);
ISImeans(:,4) = data((startrow+1):endrow,19);
ISIstds(:,4) = data((startrow+1):endrow,20);
ISIstability(:,4) = data((startrow+1):endrow,21);
data = textread('simple_model_RS2gn_LIF_vary_gandnoise_0.25noise.txt');
ISIfreqs(:,5) = data((startrow+1):endrow,18);
ISImeans(:,5) = data((startrow+1):endrow,19);
ISIstds(:,5) = data((startrow+1):endrow,20);
ISIstability(:,5) = data((startrow+1):endrow,21);
data = textread('simple_model_RS2gn_LIF_vary_gandnoise_0.125noise.txt');
ISIfreqs(:,6) = data((startrow+1):endrow,18);
ISImeans(:,6) = data((startrow+1):endrow,19);
ISIstds(:,6) = data((startrow+1):endrow,20);
ISIstability(:,6) = data((startrow+1):endrow,21);

ncellsvect = data((startrow+1):endrow,1);

figure; imagesc(log10(ISIstability(1:end,:))');
set(gca,'yticklabel',[400 200 100 50 25 12.5],'xtick',[1:6],'xticklabel',[16 32 64 128 256 512]);
xlabel('Number of cells','fontsize',12); ylabel('Noise standard deviation','fontsize',12);
title('Estimated stability time (log s)','fontsize',12)
set(gca,'fontsize',12)
h=colorbar; set(h,'fontsize',12);
colormap(flipud(gray));
set(gcf,'position',[82 404 322 275])
if exist('fig_RS2gn_varynnoise.eps')~=2
  export_fig('fig_RS2gn_varyn.eps')
  saveas(gcf,'fig_RS2gn_varyn.fig')
end
