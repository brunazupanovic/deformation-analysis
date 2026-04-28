# Geodetic Deformation Analysis (Hannover & IWST)

This repository contains a Python implementation of advanced geodetic methods for detecting and analyzing displacements in free networks.

## Overview
The project focuses on the statistical analysis of 1D (height) geodetic networks to distinguish between measurement noise and actual physical movement of object points.

### Key Features:
* **Network Adjustment:** Implements free network adjustment using the Moore-Penrose pseudoinverse.
* **Hannover Method:** Statistical approach for outlier and displacement detection.
* **IWST (Iterative Weighted Similarity Transformation):** An iterative algorithm for stable datum identification and localization of deformation.
* **Visualization:** Graphical representation of vertical displacements using `matplotlib`.

## Technical Details
The algorithm is built using:
* **NumPy:** For matrix algebra and linear systems.
* **SciPy:** For statistical F-distributions and T-tests.
* **Matplotlib:** For plotting the deformation results.

## Future Improvements
- [ ] Refactoring the code into an Object-Oriented Programming (OOP) structure.
- [ ] Integration with QGIS for spatial visualization of displacement vectors.
- [ ] Extending support to 2D/3D networks.

## Project Context
This tool was verified against the **PANDA (PAN)** software report, ensuring the mathematical correctness of the Python implementation.
