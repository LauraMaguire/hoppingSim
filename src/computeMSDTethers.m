function [msd,dtime]=computeMSD(x, tetherDist, tetherThreshold, ...
  maxpts_msd)
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

number_timepnts = size(x,3);
number_delta_t  = number_timepnts - 1;
dtime= ( 1 : number_delta_t)' ;
msd=zeros(number_delta_t,3); %Store [mean, std, n]

% find all possible zeros
numParticles = size(x,1);
try
  for dt = 1:number_delta_t
    squared_dis_ind_start = 1;
    squared_dis = zeros( numParticles * number_delta_t,1);
    for ii = 1:numParticles
      particleICurrentInds = find( tetherDist(ii,:) < tetherThreshold );
      nonOverlappingStarts = particleICurrentInds(1);
      for jj = 2:length(particleICurrentInds)
        earliestStart = particleICurrentInds(jj) - dt;
        latestStart = number_delta_t - dt;
        if particleICurrentInds(jj) <= latestStart
          for kk = 1:jj-1
            if particleICurrentInds(kk) <= earliestStart
              nonOverlappingStarts = ...
                [nonOverlappingStarts particleICurrentInds(jj)];
              break
            end
          end % kk loop
        end % if jj < small enough
      end % jj loop
      nStartPoss = unique(nonOverlappingStarts);
      nwMax = length( nStartPoss );
      numberWindow = min(nwMax,maxpts_msd);
      randInd = randperm( nwMax, numberWindow );
      index_start = nStartPoss( randInd );
      index_end = index_start + dt;
      delta_coords = x(ii,:, index_end) - x(ii,:,index_start);
      squared_dis_ind_end = squared_dis_ind_start + numberWindow-1;
      squared_dis(squared_dis_ind_start:squared_dis_ind_end) = ...
        reshape( sum(delta_coords.^2,2), [numberWindow, 1] ); % dx^2+dy^2+...
      squared_dis_ind_start = squared_dis_ind_end +1;
    end % ii, loop over particles
    % calculate displacement ^ 4 if flag
    squared_dis_accept = squared_dis( squared_dis~=0 );
    msd(dt,:) = [mean(squared_dis_accept(:)); ... % average
      std(squared_dis_accept(:)); ...; % std
      length(squared_dis_accept(:)) ]';
  end % dt
catch err
  fprintf('%s\n',err.getReport('extended') );
end
