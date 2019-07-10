function [traj_interp_yx_iTime] = interp_traj_x( traj_interp_y_iTime, Edges_tran_x)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

        numScatterEvent = size(traj_interp_y_iTime, 1);  
        traj_interp_yx_iTime = zeros(1e+5, 3);
            
        x0 = traj_interp_y_iTime(:,1);                 % extract y location
        bins = discretize(x0, Edges_tran_x);             % for each particle, compute which bin it belongs to
        difBins = diff(bins);                           % compute how many points we need to interpolate    
        curLen = 1;
        for i = 1 : numScatterEvent - 1           % interpolate trajectories between successive scattering events
            numCrossing = abs(difBins(i));      % # number of points we need to interpolate
            if (numCrossing ~= 0)                        
                v = traj_interp_y_iTime(i+1, :) - traj_interp_y_iTime(i, :);        % v: direction vector of two successive scattering events - [1,3] matrix
                if(sign(difBins(i)) < 0)                % compute z position for interpolating points along trajectories            
                    boundary_x = Edges_tran_x(bins(i): -1: bins(i+1) + 1);       % boundary_z : [1, #numCrossing]
                else
                    boundary_x = Edges_tran_x(bins(i) + 1: bins(i+1));
                end
                k = (boundary_x - traj_interp_y_iTime(i, 1)) / v(1);    % compute amplitude of interpolation vector - [1, #numCrossing]
                newDots = traj_interp_y_iTime(i, :) + k'.*v;                % add vector to the origin particle - [#numCrossing, 3]
                traj_interp_yx_iTime(curLen : curLen + numCrossing, :) = [traj_interp_y_iTime(i, :); newDots];   % store interpolated points in traj_interp
                curLen = curLen + numCrossing + 1;                         
            else                                              % only add initial particle
                traj_interp_yx_iTime(curLen, :) = traj_interp_y_iTime(i, :); 
                curLen = curLen + 1;
            end
        end
        traj_interp_yx_iTime(curLen, :) = traj_interp_y_iTime(numScatterEvent, :);    % add the last scattering event to traj_interp        
        traj_interp_yx_iTime(curLen+1:end,:) = [];

end

