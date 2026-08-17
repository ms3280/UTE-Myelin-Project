"""
reconstruction/nufft_reconstruction.py

Reconstructs 3D radial (non-Cartesian) MRI k-space data into images using
adjoint Non-Uniform FFT (NUFFT), run on HPC via torchkbnufft.

Base NUFFT usage pattern adapted from torchkbnufft's public example
notebooks (https://github.com/mmuckley/torchkbnufft); extended here for
multi-coil Siemens radial brain/phantom data, including the controlled
frequency-shift experiment used to test whether reconstruction ring
artefacts originate from off-resonance effects (see thesis Section 3.2).

Part of: UTE Myelin Project (MSc thesis, Imperial College London)
"""

# -- Set headless backend for HPC plotting --
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import torch
import torchkbnufft as tkbn
import scipy.io
from warnings import filterwarnings
filterwarnings("ignore")

# -------- Load .mat files (k-space data + trajectory) --------
kdata_dict = scipy.io.loadmat('meas_MID198_radial_TE2p197_TR25_FA16_300mmFOV_FID83463_for_hpc.mat')
traj_dict = scipy.io.loadmat('traj_interp_final_96.mat')

# -------- Extract and reshape data --------
kdata_np = kdata_dict['kspace_data']     # shape: (214, 8, 28952)
ktraj_np = traj_dict['traj_ds']          # shape: (214, 3, 28952)

print('Loaded kdata shape:', kdata_np.shape)
print('Loaded ktraj shape:', ktraj_np.shape)

# -------- Apply controlled frequency shift to k-space --------
# Used to test whether ring artefacts scale with / are caused by
# off-resonance phase accumulation during the radial readout.
dt = 1.0416667e-5  # seconds/sample (from 1000 Hz/Px x 96 base resolution)
df = 10            # frequency shift in Hz (swept over +/-15, +/-100, +/-200 Hz in the full study)

npts = kdata_np.shape[0]
t = np.arange(npts, dtype=np.float32).reshape(-1, 1, 1) * dt  # shape: (214, 1, 1)

# Apply phase modulation to simulate a frequency shift during readout
phase_shift = np.exp(1j * 2 * np.pi * df * t).astype(np.complex64)
kdata_np = kdata_np.astype(np.complex64) * phase_shift

# -------- Convert to torch tensors --------
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Reshape kdata: (npts, coils, spokes) -> (coils, spokes, npts)
kdata = np.transpose(kdata_np, (1, 2, 0))  # -> (8, 28952, 214)
kdata = torch.tensor(kdata, dtype=torch.complex64).to(device)

# Reshape ktraj: (npts, 3, spokes) -> (3, spokes, npts) -> (1, 3, N) for density compensation
ktraj_scaled = np.pi * ktraj_np / np.max(np.abs(ktraj_np))
ktraj_flat = np.transpose(ktraj_scaled, (1, 2, 0))  # -> (3, 28952, 214)
ktraj = torch.tensor(ktraj_flat, dtype=torch.float32).to(device)
ktraj_for_dcomp = ktraj.reshape(1, 3, -1)  # -> (1, 3, 28952*214)

# -------- Set image size --------
img_size = (96, 96, 96)

# -------- NUFFT operator --------
adjnufft_ob = tkbn.KbNufftAdjoint(im_size=img_size).to(device)

# -------- Density compensation --------
# Radial trajectories oversample the k-space centre, so density
# compensation weights are needed before gridding to avoid a blurred,
# centre-weighted image.
dcomp = tkbn.calc_density_compensation_function(ktraj=ktraj_for_dcomp, im_size=img_size)
dcomp = dcomp.to(torch.complex64).to(device)  # shape: (1, 1, N)
dcomp = dcomp.reshape(1, 1, ktraj.shape[1] * ktraj.shape[2])  # match (1, 1, 28952*214)

# -------- Reconstruct per coil --------
coil_imgs = []
for c in range(kdata.shape[0]):
    data_c = kdata[c].reshape(1, 1, -1)  # (1, 1, N)
    recon = adjnufft_ob(data_c * dcomp, ktraj_for_dcomp)  # (1, 1, 96, 96, 96)
    coil_imgs.append(recon)

# Stack coil images: (coils, 1, 96, 96, 96)
coil_imgs = torch.stack(coil_imgs, dim=0).squeeze(1)

# -------- Combine coils (root-sum-of-squares) --------
rss_img = torch.sqrt(torch.sum(torch.abs(coil_imgs) ** 2, dim=0))  # (96, 96, 96)

# -------- Save central slice --------
central_slice = rss_img[:, :, rss_img.shape[2] // 2].cpu().numpy()
plt.imshow(np.squeeze(central_slice) / np.max(central_slice), cmap='gray')
plt.title('Central Slice (RSS over coils)')
plt.axis('off')
plt.savefig('MID198recon_slice_rss_shift10Hz.png', bbox_inches='tight')

# -------- Save full 3D volume --------
np.save('MID198recon_volume_rss_shift10Hz.npy', rss_img.cpu().numpy())

print("Saved recon_slice_rss_shift10Hz.png and recon_volume_rss_shift10Hz.npy")
