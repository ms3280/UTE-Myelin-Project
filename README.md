# UTE Myelin Project

MSc thesis project (Imperial College London, Biomedical Engineering — Medical Physics), investigating whether combining spectral suppression with Ultrashort Echo Time (UTE) MRI can isolate the ultrashort-T2* signal from myelin in the brain.

Supervised by Dr. Pete Lally. Full methodology and results are written up in the accompanying thesis report.

## Background

Myelin — the lipid-rich sheath insulating neuronal axons — has an extremely short T2* (~0.3 ms), making it invisible to conventional MRI. This project tests whether frequency-selective suppression of the dominant water and fat signals, combined with a 3D radial UTE acquisition, can isolate the myelin signal directly, without the complex multi-component fitting required by existing techniques.

**Key finding:** neither spectral suppression nor inversion-recovery UTE fully isolated myelin — both left residual long-T2* signal that appeared as concentric ring artefacts. Controlled frequency-offset experiments confirmed these rings came from off-resonance effects during the radial readout, not suppression mistuning — motivating specific future directions (B0 field-map correction, dual-echo subtraction) discussed in the full report.

## Repository structure
```UTE-Myelin-Project/
├── README.md
├── matlab/
│   ├── load_kspace_data.m       # Reads raw Siemens .dat files via mapVBVD, extracts k-space
│   └── build_trajectory.m       # Converts scanner gradient waveforms (.dsv) into k-space trajectory
└── reconstruction/
    └── nufft_reconstruction.py  # NUFFT-based image reconstruction (PyTorch / torchkbnufft)
```

## Pipeline overview

1. **Raw data loading** (`matlab/load_kspace_data.m`) — raw Siemens `.dat` (Twix) files are parsed using the [mapVBVD](https://github.com/pehses/mapVBVD) toolbox (Ehses, P.), extracting multi-coil k-space arrays.
2. **Trajectory construction** (`matlab/build_trajectory.m`) — gradient waveform files exported from the sequence simulator are parsed and integrated to reconstruct the actual 3D k-space trajectory traversed during acquisition.
3. **Image reconstruction** (`reconstruction/nufft_reconstruction.py`) — since radial k-space data doesn't lie on a Cartesian grid, images are reconstructed using adjoint Non-Uniform FFT (NUFFT) with a Kaiser-Bessel kernel and density compensation, implemented via the [torchkbnufft](https://github.com/mmuckley/torchkbnufft) library, run on HPC.
4. **Analysis** — reconstructed volumes are used for voxel-wise SNR quantification and qualitative artefact analysis (see thesis report for full results).

## Attribution

- Raw data parsing uses the third-party **mapVBVD** toolbox (Ehses, P., MIT-licensed).
- The NUFFT reconstruction pipeline is built on the public **torchkbnufft** library and adapted from its example usage patterns for this project's specific acquisition format (multi-coil Siemens radial data), including the frequency-shift diagnostic experiments used to identify off-resonance artefacts.
- Experimental design, sequence development, data acquisition, scientific interpretation, and all discussion/conclusions are the author's own work as part of an MSc thesis at Imperial College London.

## Note on data

Raw scan data is not included in this repository, consistent with the study's institutional ethics approval (data used solely for within-project analysis). This repository contains processing and analysis code only.
