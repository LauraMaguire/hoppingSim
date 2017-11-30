function [msd,dtime]=computeMSDFixedTimeOrigin(x, startPositions, quadFlag)
%takes array x (mxdxn array), t (1xn),
% startId = 1 use start
% startId = 2 use rand
% startId = 3 use end
% calculates mean-squared and quartic displacements vs dt
% m is number of particles
% d is dimension, will work on any dimension vector
% n is number time points
%  msd is nx5 mean squared displacement array vs dt from 0 to n-1
%  msd(:,1)=msd
%  msd(:,2)=std(squared displacement)
%  msd(:,3)=n intervals(dt)
%  msd(:,4)=mean quartic displacement
%  msd(:,5)=std(quartic displacement)
runs = size(x,1);
number_timepnts = size(x,3);
newx = [];
for run=1:runs
    startIndex = startPositions{run};
    for sp = 1:length(startIndex)
        datablock = nan(1,number_timepnts);
        truncateddata = x(run,1,startIndex(sp):end);
        datablock(1:length(truncateddata)) = truncateddata;
        newx = vertcat(newx, datablock);
    end
end
xx=zeros(size(newx,1), 1,number_timepnts);
for i=1:size(newx,1)
    xx(i,1,:) = reshape( newx(i,:,1), [1 1 number_timepnts] ) ;
end

number_delta_t  = number_timepnts - 1;
dtime= ( 1 : number_delta_t)' ;

if quadFlag
  msd=zeros(number_delta_t,5); %Store [mean, std, n]
else
  msd=zeros(number_delta_t,3); %Store [mean, std, n]
end

for dt = 1:number_delta_t
  try
    index_start = 1;
    index_end = index_start + dt;
   
    delta_coords = xx(:,:, index_end) - xx(:,:,index_start);
    % calculate displacement ^ 2
    squared_dis = sum(delta_coords.^2,2); % dx^2+dy^2+...
    % calculate displacement ^ 4 if flag
    if quadFlag
      quartic_dis = sum(delta_coords.^4,2); % dx^4+dy^4+..
      msd(dt,:) = [mean(squared_dis(:)); ... % average
        std(squared_dis(:)); ...; % std
        length(squared_dis(:)); ... % n (how many points used to compute mean)
        mean(quartic_dis(:)); ... %average
        std(quartic_dis(:))]'; %std
    else
      msd(dt,:) = [nanmean(squared_dis(:)); ... % average
        nanstd(squared_dis(:)); ...; % std
        length(squared_dis(:)) ]';
    end
  catch
  end
end
