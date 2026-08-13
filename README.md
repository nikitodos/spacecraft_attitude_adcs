# Spacecraft Attitude Dynamics — 6U CubeSat ADCS

Politecnico di Milano · Spacecraft Attitude Dynamics · A.Y. 2025/26  
**Project Group: Slingshooters**

---

## Overview

This repository contains the MATLAB/Simulink implementation of an **Attitude Determination and Control System (ADCS)** for a 6U CubeSat dedicated to environmental monitoring.

The project covers spacecraft rotational dynamics and kinematics, environmental and disturbance modelling, realistic sensor models, attitude determination, eclipse propagation, nonlinear and linear control, actuator modelling, and statistical verification.

### Key features

- Rigid-body attitude dynamics and DCM kinematics
- Earth orbit, eclipse, and geomagnetic-field modelling
- Gravity-gradient, magnetic, SRP, and atmospheric-drag disturbances
- Sun sensor, magnetometer, and gyroscope models
- White, pink, and brown sensor-noise modelling
- Oustaloup fractional-order approximation
- q-Method attitude determination
- Gyroscope-based eclipse propagation
- Nonlinear de-tumbling and re-pointing control
- Linear nadir-pointing control
- Reaction-wheel actuator and saturation modelling
- Bode/Nyquist robustness analysis
- Statistical verification campaign

## Mission & Spacecraft

| Parameter | Value |
|---|---:|
| Configuration | 6U CubeSat |
| Dimensions | 30 × 20 × 10 cm |
| Mass | 12 kg |
| Orbit semi-major axis | 9371 km |
| Eccentricity | 0.2 |
| Inclination | 30° |
| Pointing requirement | < 1° nadir |

## ADCS Architecture

```text
Environment
    │
    ▼
Spacecraft Dynamics ──► Spacecraft Kinematics
    │                         │
    ├──────────────┐          │
    ▼              ▼          ▼
 Sensors       Actuators   Attitude Determination
    │              ▲          │
    └──────────────┴──────────┘
                   │
                   ▼
             Control System
        De-tumbling / Re-pointing
             / Nadir Control
```

## Verification

The final verification campaign used:

- 30 independent simulations
- 2 orbital periods per simulation
- Randomised initial conditions
- Stochastic sensor effects

The reported results show that more than 95% of simulation timestamps remain below the 1° pointing requirement, with a maximum observed error of approximately 1.25°.

## Software

- MATLAB
- Simulink
- MATLAB Control System Toolbox

## Getting Started

Clone the repository, open MATLAB, and add the source tree to the MATLAB path:

```matlab
addpath(genpath(pwd))
```

The main simulation and Simulink models are contained in `src/`.

## Repository Structure

```text
.
├── src/
├── README.md
└── SpacecraftAttitudeDynamics_FinalReport_Redatto.pdf
```

## Academic Context

**Course:** Spacecraft Attitude Dynamics  
**Institution:** Politecnico di Milano  
**Academic Year:** 2025/2026

### Team

- Stefano Diambri
- Ludovico Drioli
- Michele Pellizzer
- Giovanni Nicola D'Aloisio

## Report

The complete technical report is included in the repository.

## Disclaimer

This repository contains academic spacecraft simulation and control work. The models and results are intended for educational and engineering-analysis purposes and should not be considered flight-qualified spacecraft software.
