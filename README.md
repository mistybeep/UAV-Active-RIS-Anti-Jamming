# Anti-Jamming Transmission with UAV-Mounted Active Reconfigurable Intelligent Surfaces

[![MATLAB](https://img.shields.io/badge/MATLAB-R2023b-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Status: Core Algorithms Released](https://img.shields.io/badge/Status-Core_Algorithms_Released-orange.svg)](#)

---

## 📖 Introduction
This repository contains the MATLAB implementation for the research on utilizing Unmanned Aerial Vehicles (UAVs) equipped with active Reconfigurable Intelligent Surfaces (RIS) to combat malicious jamming attacks in multi-user cellular networks. 

Due to the inherent broadcast nature of wireless channels, communications are highly vulnerable to malicious jamming. This project investigates the synergistic roles of UAVs and active RIS to fundamentally enhance both spectral efficiency and physical-layer security.

## 🎯 System Model & Problem Formulation
![System Model](https://picui.ogmua.cn/s1/2026/03/06/69aa742358c98.webp)
The core challenge is formulated as a **worst-case Signal-to-Interference-plus-Noise Ratio (SINR) maximization** problem. To achieve robust communications and user fairness, the proposed system jointly optimizes:
1.  **Transmit Beamformer** at the Base Station (BS).
2.  **Reflection Coefficient Matrix** of the active RIS.
3.  **Deployment Location** (3D coordinates) of the UAV.

## ⚙️ Proposed Methodology
Since the formulated problem is highly coupled and non-convex, we develop an efficient **iterative optimization framework**:
* **Joint Beamforming & Reflection Design:** We propose a Successive Convex Approximation (SCA)-based Second-Order Cone Programming (SOCP) method to jointly update the transmit beamformer and the active RIS coefficients.
* **Deployment Optimization:** We utilize an SCA-based refinement algorithm to optimize the spatial deployment position of the UAV-mounted active RIS.
* **Results:** Numerical evaluations demonstrate that this scheme significantly improves system robustness and user fairness under various severe jamming conditions compared to existing benchmarks.

## 🚀 Dependencies
* MATLAB (Tested on R2023b)
* [CVX](http://cvxr.com/cvx/): MATLAB Software for Disciplined Convex Programming (required for solving SCA/SOCP problems).

---

## 📝 Citation
If you use these algorithms or ideas in your research, please cite:

```bibtex
@ARTICLE{11479934,
  author={Wang, Hao and Guo, Yifan and Luo, Junshan and Wang, Shilian},
  journal={IEEE Transactions on Vehicular Technology}, 
  title={Anti-Jamming Transmission With UAV-Mounted Active Reconfigurable Intelligent Surfaces}, 
  year={2026},
  volume={},
  number={},
  pages={1-6},
  keywords={Jamming;Electronic warfare;Antennas;Antennas and propagation;Broadcasting;Broadcast technology;Thermal noise;Communication systems;Internet of Things;Wireless communication;Unmanned aerial vehicle;active reconfigurable intelligent surface;anti-jamming;max-min fairness},
  doi={10.1109/TVT.2026.3682788}}

