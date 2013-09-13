% clean up and optionally upsample Hafting et al trajectory
% generate x,y trajectory at time points from 0 to simdur in increments of
% trajdt, returns dx and dy at each time point: cumsum(hafting_trajectory(simdur,trajdt))
% gives actual x and y positions
% eric zilli - oct 20, 2009
function [dxs,dys,fdxs,fdys] = hafting_trajectory(simdur,trajdt)

load rat_10925.mat

% grab this amount of trajectory form the Hafting data:
% can't be >600 s because Hafting data is 3 concatenated trajectories of 600
% s each

% the trajectory is scaled in time by this rate s.t. posskip=2 means the
% rat moves along the new trajectory at twice the rate of the original
posskip = 1;
dscale = 1;

trajts = 0:trajdt:simdur; % ms

endind = find(pos_timeStamps>simdur/1e3,1,'first');

dxs = diff(pos_x(1:posskip:posskip*endind))/100; % m
dys = diff(pos_y(1:posskip:posskip*endind))/100; % m
ts = pos_timeStamps(1:endind-1)*1e3; % ms

xinit = pos_x(1);
yinit = pos_y(1);

% resample:
if trajdt~=(ts(2)-ts(1))
  dxs = dscale*interp1(ts,dxs,trajts)*(trajdt/(ts(2)-ts(1))); % m
  dys = dscale*interp1(ts,dys,trajts)*(trajdt/(ts(2)-ts(1))); % m
end

% optionally filter trajectory
[fb,fa] = butter(3,.2/((1000/trajdt)/4),'low');
fdxs = filtfilt(fb,fa,dxs);
fdys = filtfilt(fb,fa,dys);

% report some information
disp(sprintf('mean unfiltered speed: %g m/s',mean(sqrt(dxs.^2+dys.^2)/trajdt*1000)))
disp(sprintf('mean filtered speed: %g m/s',mean(sqrt(fdxs.^2+fdys.^2)/trajdt*1000)))
% return
% figure; plot(xinit+cumsum(dxs),yinit+cumsum(dys),'.'); %set(gca,'xlim',[0 xmax],'ylim',[0 ymax]); xlabel('position (m)'); ylabel('position (m)')
% figure; plot(xinit+cumsum(fdxs),yinit+cumsum(fdys),'.'); %set(gca,'xlim',[0 xmax],'ylim',[0 ymax]); xlabel('position (m)'); ylabel('position (m)')
% figure; plot(xinit+cumsum(fdxs(1:(end/4))),yinit+cumsum(fdys(1:(end/4))),'.'); %set(gca,'xlim',[0 xmax],'ylim',[0 ymax]); xlabel('position (m)'); ylabel('position (m)')
% % figure; plot(xs,ys,'.')