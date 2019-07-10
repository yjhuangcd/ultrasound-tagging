function OPL = compute_opl_normal(traj_interp, Edges_y, Edges_tran_x, Edges_tran_z, USpressure_transverse, USpressure_prop, nTime, RI_Coeff, n0)
%COMPUTE_OPL compute OPL given interpolated moving trajectories
%   INPUT:
%       traj_interp: interpolated trajs
%       Edges: boundaries of discretized layers
%       USpressure_prop: ultrasound pressure map @poisition and time
%       nTime: number of time stamps
%       RI_Coeff: delta_n = RI_Coeff * ultrasound pressure
%   OUTPUT:
%       OPL: [1, #tseq] mat

        segmat = traj_interp(2:end, :, :) - traj_interp(1:end-1, :, :);      % vector matrix between two successive interpolated points
        lenSegs = squeeze(sqrt(sum(segmat.^ 2, 2)));                       % [n,1,m] -> [n,m] mat, n - #segments, m - time
        midys = (traj_interp(2:end, 2, :) + traj_interp(1:end-1, 2, :)) / 2;    % compute midpoint of each segment, [n, m] mat
        midzs = (traj_interp(2:end, 3, :) + traj_interp(1:end-1, 3, :)) / 2; 
        midxs = (traj_interp(2:end, 1, :) + traj_interp(1:end-1, 1, :)) / 2; 
        bins_y = discretize(midys, Edges_y);                                           % compute bins for each segment, [n, m] mat
        bins_z = discretize(midzs, Edges_tran_z);
        bins_x = discretize(midxs, Edges_tran_x);
        numSegs = size(bins_y,1);
        
        % compute index mat to extract corresponding USpressure_prop @ each time and position
        binsExtract_y = reshape(bins_y, numel(bins_y), 1);                          % index along n direction (position)
        indTimeTemp = repmat(1:round(nTime), numSegs,1);            
        indTime = reshape(indTimeTemp, numel(indTimeTemp), 1);
        indExtract = sub2ind(size(USpressure_prop), binsExtract_y, indTime);
        Pressure_seg = reshape(USpressure_prop(indExtract), size(lenSegs));  % [n, m] mat
        
        % Transverse Gaussian field has a Gaussian profile
        binsExtract_z = reshape(bins_z, numel(bins_z), 1);
        binsExtract_x = reshape(bins_x, numel(bins_x), 1);
        indExtract_transverse = sub2ind(size(USpressure_transverse), binsExtract_x, binsExtract_z);
        Pressure_trans_seg = reshape(USpressure_transverse(indExtract_transverse), size(lenSegs));
        
        OPL = sum(lenSegs.* (1 + RI_Coeff * Pressure_seg.*Pressure_trans_seg), 1) *n0;
        
end

