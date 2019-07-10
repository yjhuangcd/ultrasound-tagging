function [ ftraj_tseq ] = move_trajs_normal(ftraj, tseq, wa, ka, A_dis)
%Generate trajectories over time sequence based on baseline
%   Input: 
%       ftraj: filtered trajs that corresponds to one photon
%       tseq: time sequence
%       wa, ka: ultrasound frequency
%       A_dis: particle displacement

        ftraj = double(ftraj);                    % For displacement as small as 0.1 nm, need to use double precision
        y0move = ftraj(2:end-1, 2);         % Fix the starting point and the ending point
        ftraj_tseq = repmat(ftraj, [1,1, length(tseq)]);
        yt = y0move + A_dis * sin(wa*tseq - ka*y0move + pi/2*3);     % phase diffrence: 270 degrees
        ftraj_tseq(2:end-1,2,:) = yt;

end

