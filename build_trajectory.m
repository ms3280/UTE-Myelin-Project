%% matlab/build_trajectory.m
%
% Reconstructs the actual 3D k-space trajectory traversed during
% acquisition, by reading gradient waveform files exported from the
% sequence simulator (.dsv format) and integrating them to recover
% gradient-time-area (i.e. k-space position) along each axis.
%
% Timing constants below (spacing, first_tr, ro_len) were determined
% empirically for this specific sequence by inspecting the gradient
% waveform plots.
%
% Part of: UTE Myelin Project (MSc thesis, Imperial College London)

clear all
startdir = pwd;

%% Folder for simulation files
sim_name = fullfile(getenv('HOME'), 'Desktop', 'mri', 'scan data', 'phantom'); % <-- adjust path
cd(sim_name);

%% Read in X, Y, Z gradient coordinates
[XG] = read_gradient_axis('*_M1X.dsv');
[YG] = read_gradient_axis('*_M1Y.dsv');
[ZG] = read_gradient_axis('*_M1Z.dsv');

%% Build 3D array of k-space points: readout index x axis x TR
% Constants below were determined empirically from the gradient plots
spacing  = 10002; % number of points per TR (spacing between readouts)
first_tr = 933;   % sample index at which the first gradient readout begins
ro_len   = 516;   % length of a single readout, in samples

traj = [];
for ro_idx = 0:round(length(ZG)/spacing)-1
    traj(1:ro_len, :, ro_idx+1) = [ ...
        XG(first_tr+(ro_idx*spacing+(1:ro_len))).', ...
        YG(first_tr+(ro_idx*spacing+(1:ro_len))).', ...
        ZG(first_tr+(ro_idx*spacing+(1:ro_len))).'];
end

% Visualise the end-point of the trajectory for each readout (sanity check)
scatter3(squeeze(traj(end,1,:)), squeeze(traj(end,2,:)), squeeze(traj(end,3,:)), '.')
axis square equal
title('K-space trajectory endpoints (sanity check)')

%% Downsample trajectory to match acquired ADC points and save
adc_points = 66;
stride = (size(traj, 1) / adc_points);
traj_ds = traj(1:stride:end, :, :);
traj_single = single(traj_ds);
save('traj_nufft32.mat', 'traj_single');


%% ---- Local helper function ----
function G = read_gradient_axis(pattern)
    % Reads one gradient axis (.dsv file) and integrates it to recover
    % the cumulative gradient-time-area (i.e. k-space trajectory) for
    % that axis. The .dsv format stores run-length-encoded gradient
    % samples, which this function decodes before integrating.

    vf_txt = split(regexp(fileread(dir(pattern).name), "[^\n]*VERTFACTOR[^\r]*", "match"), "=");
    vert_scl = str2double(vf_txt{2});
    temp_mat = readmatrix(dir(pattern).name, "FileType", "text");

    G = [];
    G(1) = 0;
    midx = 1;
    gidx = 1;

    while midx < (size(temp_mat,1) - 1)
        midx = midx + 1;
        gidx = gidx + 1;
        G(gidx) = temp_mat(midx);
        if temp_mat(midx) == temp_mat(midx-1)
            G(gidx+1:gidx+temp_mat(midx+1)+1) = G(gidx);
            gidx = gidx + temp_mat(midx+1) + 1;
            G(gidx) = temp_mat(midx+2);
            midx = midx + 2;
        end
    end

    G = cumsum(G(1:end-1)) ./ vert_scl;
end
