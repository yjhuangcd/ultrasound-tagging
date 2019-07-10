%% Simulate ultrasound tagging efficiency in scattering medium giving trajectories
% Ultrasound propagates normal to the light propagation direction
% by Yujia Huang, yjhuang@caltech.edu

%% Parameters
% Ultrasound covered region dimention
% Define the incidient and exit planes
para.Inc_y = 0;         
para.Exit_y = 5;
para.Inc_z = 0;          % Transverse US field (light propagation direction: z)
para.Exit_z = 10;
para.Inc_x = 0;          % Transverse US field (x)
para.Exit_x = 5;  
para.FWHM = 4.2;      % FWHM of ultrasound field, unit: mm
para.distancePerLayer_y = 0.02;   % Discretizing step size when integrating for OPL, unit: mm
para.numLayers_y = (para.Exit_y - para.Inc_y)/para.distancePerLayer_y;
para.distancePerLayer_x = 1;   % Discretizing step size for the Gaussian lateral profile, unit: mm
para.distancePerLayer_z = 1;   % Discretizing step size for the Gaussian lateral profile, unit: mm
para.fsample = 50.25*1e+6;      % Sample frequency: 25.25MHz
para.tlength = 4*1e-6;               % Length of time sequence, (s) 4us, CW
para.lambda = 0.532*1e-3;        % optical wavelength: mm (1e-3 mm = 1 um)
para.fa = 1*1e+6;          % 1 MHz, 2.25 MHz, 3.5 MHz 
para.va = 1.480e+6;      % ultrasound velocity, unit: mm/s
para.rou = 1000;           % medium density, unit: kg/m^3
para.n0 = 1.33;             % medium refractive index
para.reptime = 50;

interval   = para.tlength * para.reptime * para.fa;            %frequency-space resolution                           
n           = para.tlength * para.reptime * para.fsample;    % total number of points in power spectrum

f = waitbar(0,'1','Name','Computing Power Spectrum...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f,'canceling',0);

%% Calculting power spectrum
% load the trajectories
% file format: trajectory id | x | y | z
load('.\data\processed_data_mus_4.mat');
% load the ultrasound pressure
load('.\Calibrated-Ultrasound-Field\Pressure_MPa_1MHz');
% calculate the displacement
dis_list = Pressure_MPa*100/(para.fa/1e+6);
for dis_i = 1:length(dis_list)    
    para.A_dis = dis_list(dis_i)*1e-6;    % unit: mm
    
    [~, StartID] = unique(processed_trajs(:,1));
    addr = [StartID ; size(processed_trajs,1) + 1];

    [AvgPowerSpec, stats, err] = cal_powerspec_normal_period_waitbar_errtrend(para, processed_trajs, addr, f);

    filename = strcat('.\Results_dis_', num2str(dis_list(dis_i)), '.mat');
    save(filename, 'AvgPowerSpec', 'para', 'stats', 'addr', 'err');
end
delete(f);

%%
tagEff = zeros(length(dis_list),1);
for dis_i = 1:length(dis_list)
    filename = strcat('.\Results_dis_', num2str(dis_list(dis_i)), '.mat');
    load(filename);
    zeroOrder = n/2+1;
    tagEff(dis_i) = 1-AvgPowerSpec(zeroOrder);
end
figure; plot(Pressure_MPa, tagEff, 'm-o'); ylim([0,1]);
xlabel('Pressure (MPa)'); ylabel('Tagging Efficiency');

