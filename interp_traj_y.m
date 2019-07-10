function [ traj_interp_y_iTime ] = interp_traj_y( ftraj_tseq_iTime, Edges_y )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

        numScatterEvent = size(ftraj_tseq_iTime, 1);  
        traj_interp_y_iTime = zeros(1e+5, 3);
            
        y0 = ftraj_tseq_iTime(:,2);                 % extract y location
        bins = discretize(y0, Edges_y);             % for each particle, compute which bin it belongs to
        difBins = diff(bins);                           % compute how many points we need to interpolate    
        curLen = 1;
        for i = 1 : numScatterEvent - 1           % interpolate trajectories between successive scattering events
            numCrossing = abs(difBins(i));      % # number of points we need to interpolate
            if (numCrossing ~= 0)                        
                v = ftraj_tseq_iTime(i+1, :) - ftraj_tseq_iTime(i, :);        % v: direction vector of two successive scattering events - [1,3] matrix
                if(sign(difBins(i)) < 0)                % compute z position for interpolating points along trajectories            
                    boundary_y = Edges_y(bins(i): -1: bins(i+1) + 1);       % boundary_z : [1, #numCrossing]
                else
                    boundary_y = Edges_y(bins(i) + 1: bins(i+1));
                end
                k = (boundary_y - ftraj_tseq_iTime(i, 2)) / v(2);    % compute amplitude of interpolation vector by z position - [1, #numCrossing]
                newDots = ftraj_tseq_iTime(i, :) + k'.*v;                % add vector to the origin particle - [#numCrossing, 3]
                traj_interp_y_iTime(curLen : curLen + numCrossing, :) = [ftraj_tseq_iTime(i, :); newDots];   % store interpolated points in traj_interp
                curLen = curLen + numCrossing + 1;                         
            else                                              % only add initial particle
                traj_interp_y_iTime(curLen, :) = ftraj_tseq_iTime(i, :); 
                curLen = curLen + 1;
            end
        end
        traj_interp_y_iTime(curLen, :) = ftraj_tseq_iTime(numScatterEvent, :);    % add the last scattering event to traj_interp
        traj_interp_y_iTime(curLen+1:end,:) = [];

end

