%    2018/07/30

%     Return stats of used photons
%     When calculating trajectory changes due to particle motion, fix the starting point and ending point of the trajectory
%     Ultrasound and light propagates normal to each other (The result of running this should equals to flip the y and z dimension of parallel case and change inc_z and exit_z to the y range)
%     Only calculate OPL change in a period and rep mat to get all the OPL trace
%     Calculate bin for trajs at each time t.

%   Input:
%       para: optical, sample, ultrasound properties
%       processed_trajs: trajs mat with size of [N*4]. 
%              N: # of scattering event * # of photons
%              1: photon ID  2-4: x,y,z position
%   Output:
%       AvgPowerSpec: Averaged power spectrum
%       stats: record OPL and tagging efficiency for each photon

function [ AvgPowerSpec, stats, err] = cal_powerspec_normal_period_waitbar_errtrend( para, processed_trajs, addr, f)

% Parameters
    numLayers_y    = para.numLayers_y;             
    Inc_y            = para.Inc_y;                         % Ultrasound propagation direction
    Exit_y            = para.Exit_y;
    Inc_z            = para.Inc_z;                      % Tranverse direction (Light propagation: z)
    Exit_z           = para.Exit_z;
    Inc_x            = para.Inc_x;                       % Transverse-direction (x)
    Exit_x           = para.Exit_x;
    sigma          = para.FWHM/2;                    % Transverse ultrasound field FWHM
    fsample         = para.fsample;                   % Sample frequency: 10MHz
    tlength           = para.tlength;                       % Length of time sequence, (s) 100us = tPeriod
    lambda          = para.lambda;                    % optical wavelength: mm (1e-3 mm = 1 um)
    fa          = para.fa;             % 1MHz
    va          = para.va;            % mm/s 
    rou        = para.rou;           % kg/m^3
    A_dis     = para.A_dis;       % mm, 100nm
    n0          = para.n0;           % Need to be the same as trajectory generating code
%     interval   = para.tlength * para.reptime * para.fa;   %frequency-space resolution                           
    n           = para.tlength * para.reptime * para.fsample;    % total number of points in power spectrum  
    nPeriod  =  para.tlength * para.fsample;          % number of time stamps in computed periods
    wa         = 2*pi*fa;                
    ka          = wa/va;
    k0 = 2*pi / lambda;            % k vector of light
    eta         = 1.466*1e-10 * rou * (para.va/1000)^2;      % no unit, Lihong's paper, using water as background, need to find parameters for agar
    RI_Coeff = eta * ka * A_dis;       % US on
    %RI_Coeff = 0;          %US off, only consider displacement
    tseq        = 0 : (1/fsample) : (tlength - 1/fsample); 
    AvgPowerSpec = 0;    
    
%     NumofUsedPhotons = length(addr) - 1;
    NumofUsedPhotons = 500;
    err = zeros(NumofUsedPhotons,1);

     % Preparation
     Edges_y = linspace(Inc_y - A_dis, Exit_y + A_dis, numLayers_y + 1);   % Discretizing along ultrasound propagation direction y (normal to the light propgation direction)
     ys = (Edges_y(1:end-1) + Edges_y(2:end)) / 2;                                    % mid-z position of each layer
     [Time, Position] = meshgrid(tseq, ys);                                            % for each position in the volume and each time
     USpressure_prop = sin(wa*Time - ka*Position);                                      % compute ultrasound pressure as a sine wave
     % Transverse Ultrasound field has Gaussian profile
     Edges_tran_z = Inc_z : para.distancePerLayer_z : Exit_z;
     zs = (Edges_tran_z(1:end-1) + Edges_tran_z(2:end)) / 2;
     Edges_tran_x = Inc_x : para.distancePerLayer_x : Exit_x;
     xs = (Edges_tran_x(1:end-1) + Edges_tran_x(2:end)) / 2;
     [transversePos_z, transversePos_x] = meshgrid(zs, xs);
     USpressure_transverse = exp(-((transversePos_x - mean(Edges_tran_x)).^2 + (transversePos_z - mean(Edges_tran_z)).^2)./(2*sigma^2));

       for iPhoton = 1 : NumofUsedPhotons
           % wait bar
           if getappdata(f,'canceling')
               break
           end    
            % Update waitbar and message
           waitbar(iPhoton/NumofUsedPhotons,f,sprintf('%d',iPhoton));

           ftraj = processed_trajs(addr(iPhoton) : addr(iPhoton + 1) - 1, 2:4);        

            % 1 - Generate trajectories over time sequence based on baseline
            ftraj_tseq = move_trajs_normal(ftraj, tseq, wa, ka, A_dis);
%             ftraj_tseq = move_trajs_normal_unfixSE(ftraj, tseq, wa, ka, A_dis);

            % 2 - Interpolation given number of layers 
            traj_interp = interp_trajs_normal(Edges_y, Edges_tran_x, Edges_tran_z, ftraj_tseq);           

            % 3 - Compute OPL
            OPL = compute_opl_normal(traj_interp, Edges_y, Edges_tran_x, Edges_tran_z, USpressure_transverse, USpressure_prop, nPeriod, RI_Coeff, n0);
            
            % 3 - Compute power Spectrum
            powerSpec = compute_powerSpec(OPL, k0, para.reptime);
            
            % 4- Ensemble powerSpec            
            AvgPowerSpec = AvgPowerSpec + powerSpec;
            
            % See whether the power spectrum convergent or not
            err(iPhoton) = sum(abs((AvgPowerSpec-powerSpec)/(iPhoton-1) - AvgPowerSpec/iPhoton));
                   
%             peakPos_order1 = uint32([n/2+1- interval, n/2+1 + interval]);            
            
            stats.meanOPLs(iPhoton) = mean(OPL);
            stats.rangeOPLs(iPhoton) = range(OPL);
            stats.numScatter(iPhoton) = size(ftraj,1) - 2;
%             stats.I1s(iPhoton) = sum(powerSpec(peakPos_order1));
            stats.I0s(iPhoton) = powerSpec(uint32(n/2)+1);       
       end
       
       AvgPowerSpec = AvgPowerSpec./NumofUsedPhotons;
  
end

% valid = 0;   % how many power spectrums are symmetric
% if(abs(powerSpec(peakPos_order1(1)) - powerSpec(peakPos_order1(2))) < 0.01)
%     valid = valid + 1;
%     AvgPowerSpec = AvgPowerSpec + powerSpec;
% else
%     iPhoton
% end
