%% matlab/load_kspace_data.m
%
% Reads raw Siemens .dat (Twix) k-space data using the third-party
% mapVBVD toolbox (Ehses, P. -- https://github.com/pehses/mapVBVD),
% then extracts and trims the multi-coil k-space array for reconstruction.
%
% Part of: UTE Myelin Project (MSc thesis, Imperial College London)

addpath('/Users/mahima/Desktop/mri/mapVBVD-main');  % <-- adjust path to your local mapVBVD install

%% Load raw data from Siemens .dat file (Twix file)
measFile = '/Users/mahima/Desktop/mri/scan data/phantom/meas_MID340_radial_fwsuppressed_tr6p05_te0p196_FID81532.dat'; % <-- adjust path as needed
twix = mapVBVD(measFile);

%% Extract full k-space data
kspace_data = twix.image.unsorted();

%% Remove first 2 ADC points along first dimension (readouts)
% These initial points occur before the gradient trajectory has
% stabilised and can contaminate the k-space data (see thesis Section 2.2.1).
kspace_data = kspace_data(3:end, :, :);
[adc_len, n_coils, n_spokes_short] = size(kspace_data);

%% Save trimmed k-space data
save('/Users/mahima/Desktop/mri/scan data/phantom/3/noise_vol3.mat', 'kspace_data');
fprintf('Trimmed k-space data saved as kspace_data.mat\n');

%% Explore header (check MeasYaps block for trajectory info)
disp('--- Trajectory / K-space Info ---');
try
    disp(twix.hdr.MeasYaps.sKSpace);
catch
    disp('sKSpace field not found.');
end
