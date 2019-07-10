function [traj_interp_yxz_iTime] = interp_traj_z( traj_interp_yx_iTime, Edges_tran_z)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

        numScatterEvent = size(traj_interp_yx_iTime, 1);  
        traj_interp_yxz_iTime = zeros(1e+5, 3);
            
        z0 = traj_interp_yx_iTime(:,3);                 % extract y location
        bins = discretize(z0, Edges_tran_z);             % for each particle, compute which bin it belongs to
        difBins = diff(bins);                           % compute how many points we need to interpolate    
        curLen = 1;
        for i = 1 : numScatterEvent - 1           % interpolate trajectories between successive scattering events
            numCrossing = abs(difBins(i));      % # number of points we need to interpolate
            if (numCrossing ~= 0)                        
                v = traj_interp_yx_iTime(i+1, :) - traj_interp_yx_iTime(i, :);        % v: direction vector of two successive scattering events - [1,3] matrix
                if(sign(difBins(i)) < 0)                % compute z position for interpolating points along trajectories            
                    boundary_z = Edges_tran_z(bins(i): -1: bins(i+1) + 1);       % boundary_z : [1, #numCrossing]
                else
                    boundary_z = Edges_tran_z(bins(i) + 1: bins(i+1));
                end
                k = (boundary_z - traj_interp_yx_iTime(i, 3)) / v(3);    % compute amplitude of interpolation vector - [1, #numCrossing]
                newDots = traj_interp_yx_iTime(i, :) + k'.*v;                % add vector to the origin particle - [#numCrossing, 3]
                traj_interp_yxz_iTime(curLen : curLen + numCrossing, :) = [traj_interp_yx_iTime(i, :); newDots];   % store interpolated points in traj_interp
                curLen = curLen + numCrossing + 1;                         
            else                                              % only add initial particle
                traj_interp_yxz_iTime(curLen, :) = traj_interp_yx_iTime(i, :); 
                curLen = curLen + 1;
            end
        end
        traj_interp_yxz_iTime(curLen, :) = traj_interp_yx_iTime(numScatterEvent, :);    % add the last scattering event to traj_interp        
        traj_interp_yxz_iTime(curLen+1:end,:) = [];
end

