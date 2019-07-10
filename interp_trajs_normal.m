function traj_interp = interp_trajs_normal(Edges_y, Edges_tran_x, Edges_tran_z, ftraj_tseq)
%INTERP_TRAJS Interpolate moving trajectories in order to compute integration
%   INPUT:
%       Inc_y & Exit_y: incident plane and exit plane -> volume that needs to be interpolated
%       A_dis: amplitude of particle movement
%       numLayers: number of layers that we use to discretize the volume
%       ftraj_tseq: moving trajectories [m,n,k] mat
%                       m - # of scattering event
%                       n - 3, [x,y,z]
%                       k - # of time stamp
%   OUTPUT:
%       traj_interp: interpolated trajs at each time stamp [M, n, k] mat
%                       M - # of interpolations

            traj_interp = zeros(1e+5, 3, size(ftraj_tseq,3));
            numInterp = zeros(size(ftraj_tseq, 3) ,1);
            
            for iTime  = 1:size(ftraj_tseq, 3)             % for each time
                ftraj_tseq_iTime = ftraj_tseq(:,:,iTime);
                traj_interp_y_iTime = interp_traj_y( ftraj_tseq_iTime, Edges_y );
                traj_interp_yx_iTime = interp_traj_x( traj_interp_y_iTime, Edges_tran_x);
                traj_interp_yxz_iTime = interp_traj_z( traj_interp_yx_iTime, Edges_tran_z);
                numInterp(iTime) = size(traj_interp_yxz_iTime, 1);
                traj_interp(1:numInterp(iTime), : ,iTime) = traj_interp_yxz_iTime;
            end
            
            maxNumInterp = max(numInterp);                   % maximum number of interpolations
            indPad = find(numInterp < maxNumInterp);     % find those time stamp that contains fewer interpolation points
            for iPad = 1:length(indPad)                            % pad traj_interp to assure that for each time, number of points are the same for further processing
                iTPad = indPad(iPad);
                traj_interp(numInterp(iTPad)+1 : maxNumInterp, :, iTPad) = repmat(traj_interp(numInterp(iTPad), :, iTPad), maxNumInterp - numInterp(iTPad),1);
            end
            
            traj_interp(maxNumInterp+1 : end, :, :) = [];

end

