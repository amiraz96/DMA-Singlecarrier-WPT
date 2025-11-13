# DMA-Singlecarrier-WPT

MATLAB reference implementation that accompanies the paper "Energy Beamforming for RF Wireless Power Transfer With Dynamic Metasurface Antennas" (https://ieeexplore.ieee.org/abstract/document/10361540). The scripts reproduce the alternating semidefinite-programming routine proposed in the paper to design the beamforming of a dynamic metasurface antenna (DMA) for multi-user wireless power transfer.

## Repository structure

| File | Description |
| --- | --- |
| `DMA_optimization.m` | Main entry point. Sets up the system geometry, synthesizes channel responses, and runs the two-stage semidefinite programs that optimize the DMA radio-frequency chain (RFC) weights and programmable element responses under user power constraints. |
| `DMA_deploy.m` | Generates the spatial coordinates of the DMA's RFC feeds and passive radiating elements, and evaluates the per-element DMA transfer function. |
| `Do_Channels.m` | Computes the free-space near-field channels between each programmable DMA element and each user location, including phase accumulation and amplitude tapering based on boresight gain. |
| `H_DMA.m` | Evaluates the complex Lorentzian response of the DMA waveguide network as a function of frequency, aperture length, and element offset. |
| `Q_DMA.m` | Implements the Lorentzian-constrained relationship between the DMA's tunable phase parameter and the complex element coefficient. |

## Prerequisites

* MATLAB R2021b or newer (the scripts rely on recent syntax such as implicit expansion and `physconst`).
* [CVX](http://cvxr.com/cvx/) convex optimization toolbox with an SDP-capable solver (e.g., SeDuMi or MOSEK). Install CVX and run `cvx_setup` in MATLAB before executing the scripts.
* Optional: Parallel Computing Toolbox if you plan to extend the code with parallel Monte Carlo channel evaluations.

## Quick start

1. Open `DMA_optimization.m` in MATLAB.
2. Adjust the **System Parameters** section to match your operating frequency, aperture length, transmit power budget, and deployment region.
3. Specify user locations (`user_loc`) and received power thresholds (`RF_Thr`) for each energy receiver. Coordinates are expressed in meters.
4. Run the script. CVX will solve two semidefinite programs:
   * the first jointly optimizes the RFC weights (`W`) subject to the radio-frequency thresholds;
   * the second refines the programmable element coefficients (`Q2`) under the Lorentzian constraint enforced through the lookup stage.
5. Inspect the resulting optimized complex coefficients (`w`, `Q`, `qfinal`) to configure your DMA hardware or to simulate performance.

The code currently instantiates a single-user scenario (`M = 1`) with one DMA centred in a 10 m × 10 m × 3 m volume. Extend the `user_loc` matrix and power thresholds to cover multi-user cases, and update the region bounds and DMA placement as needed.

## Extending the scripts

* **Alternative propagation models:** Replace `Do_Channels.m` with a model that incorporates blockage, misalignment, or measured channel data if you are targeting indoor deployments.
* **Hardware constraints:** Modify `Q_DMA.m` to incorporate additional amplitude constraints or quantization steps that match your metasurface tuning hardware.
* **Performance evaluation:** Log the received power `trace(W'*reshape(Bk(k,:,:), [N_d, N_d]))` across iterations to generate the figures of merit reported in the paper.

## Citing

If you use this code, please cite the associated IEEE Transactions on Wireless Power Transfer article. A BibTeX entry is provided below; update the metadata once you import it from IEEE Xplore.

```bibtex
@article{DMAWPT2023,
  title={Single-Carrier Wireless Power Transfer With Dynamic Metasurface Antennas},
  author={Authors et al.},
  journal={IEEE Transactions on Wireless Power Transfer},
  year={2023},
  doi={10.1109/TWPT.2023.10361540}
}
```

## License

The repository inherits the usage rights granted by the original authors. Refer to the paper or contact the authors for explicit licensing terms.
