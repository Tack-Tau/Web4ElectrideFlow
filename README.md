# Inorganic Electride Database

A web application for browsing and analyzing inorganic electride structures discovered through the [ElectrideFlow](https://github.com/MaterSim/ElectrideFlow) high-throughput workflow.

**Live Demo**: [inorganic-electride.onrender.com](https://inorganic-electride.onrender.com)

![Database Statistics](https://img.shields.io/badge/Structures-290-purple)
![Python](https://img.shields.io/badge/Python-3.10+-blue)
![License](https://img.shields.io/badge/License-MIT-yellow)

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [How to Use](#how-to-use)
- [Data Description](#data-description)
- [Local Installation](#local-installation)
- [Deployment](#deployment)
- [Citation](#citation)

---

## Overview

This database contains **290 unique inorganic electride structures** (251 binary + 39 ternary) identified through the [ElectrideFlow](https://github.com/MaterSim/ElectrideFlow) workflow ([arXiv:2601.21077](https://arxiv.org/abs/2601.21077)).

### What are Electrides?

Electrides are ionic compounds where electrons occupy anionic sites in the crystal lattice, leading to unique properties:
- Low work functions
- High electron mobility  
- Catalytic activity
- Superconductivity potential

### Workflow

The database structures were discovered using a 4-stage hierarchical workflow:

1. **Structure Generation**: 50,000+ candidate structures from substitution-based prototypes
2. **MatterSim Pre-screening**: Rapid ML-based energy evaluation and stability filtering
3. **VASP DFT Calculations**: Accurate electronic structure and energy calculations
4. **Electride Analysis**: Automated detection of interstitial electron density

---

## Features

- **Browse & Filter**: Search by chemical system, composition, energy, space group
- **3D Visualization**: Interactive JSmol viewer with unit cell controls
- **Convex Hull Plots**: Phase stability diagrams with MP reference data
- **Export Data**: Download structures in CIF, XYZ, JSON formats
- **Database Downloads**: Export full database as DB, CSV, or JSON

---

## How to Use

### Browsing the Database

1. Visit [inorganic-electride.onrender.com](https://inorganic-electride.onrender.com)
2. Select **Binary** (251 structures) or **Ternary** (39 structures)
3. Use filters or search bar to find specific structures
4. Click on any structure ID to view details

### Viewing Structure Details

Each structure page shows:
- **Crystal Structure**: Interactive 3D viewer (rotate with mouse, zoom with wheel)
- **Properties**: Space group, band gap, energy above hull, electride characteristics
- **Convex Hull Plot**: Phase stability diagram for the chemical system
- **Export Options**: Download structure in various formats

### Searching by Chemical System

Enter a chemical system in the search bar:
- Binary: `Ca-P`, `Li-N`, etc.
- Ternary: `K-B-O`, `Al-Ca-P`, etc.

The database automatically switches between binary/ternary based on your query.

### Understanding the Data

Key properties displayed:
- **E_hull (DFT)**: Energy above convex hull from VASP calculations (eV/atom)
  - < 0.001: Thermodynamically stable
  - < 0.05: Potentially synthesizable
- **N_excess**: Number of excess electrons in interstitial sites
- **Band Gap**: Electronic band gap (eV)
- **Pearson Symbol**: Structure type classification

---

## Data Description

Each structure includes:

**Basic Properties**
- Composition, space group, Pearson symbol
- Crystal structure (atomic positions, lattice parameters)

**Electronic Properties**  
- Band gap (eV)
- Electron density at interstitial sites (e0.025, e0.5, e1.0)
- Bader charge analysis (band0, band1)

**Thermodynamic Properties**
- MatterSim energy and stability
- VASP DFT energy and stability  
- Energy above convex hull (full phase diagram)

**Electride Characteristics**
- Number of excess electrons (N_excess)
- Interstitial site locations and volumes
- Electride confirmation flag

### Database Files

- **Binary Systems**: 251 unique structures from Ca-P, Li-N, K-Cl, etc.
- **Ternary Systems**: 39 unique structures from K-B-O, Al-Ca-P, etc.
- **Reference Data**: Materials Project phases for convex hull analysis

---

## Local Installation

### Prerequisites

- Python 3.10+
- Git

### Setup

```bash
# Clone repository
git clone https://github.com/YOUR-USERNAME/Web4ElectrideFlow.git
cd Web4ElectrideFlow

# Create environment
conda create -n electride-web python=3.10
conda activate electride-web

# Install dependencies
pip install -r requirements.txt

# Run locally
./start.sh

# Open browser to http://localhost:5000
```

---

## Deployment

### Quick Deploy to Render.com

1. **Push to GitHub**:
   ```bash
   git init
   git add .
   git commit -m "Initial commit"
   git push -u origin main
   ```

2. **Create Render Web Service**:
   - Go to [dashboard.render.com](https://dashboard.render.com/)
   - New → Web Service → Connect your GitHub repo
   - Settings:
     - Build Command: `pip install -r requirements.txt`
     - Start Command: `./start.sh`
     - Add Environment Variable: `SECRET_KEY` (generate with `python3 -c "import secrets; print(secrets.token_hex(32))"`)

3. **Upload Database Files** (if not using Git LFS):
   - Via Render Shell or persistent disk (see deployment guide in repo)

4. **Access at**: `https://your-app-name.onrender.com`

For detailed deployment instructions, see the deployment guide in the repository.

---

## Citation

If you use this database in your research, please cite:

**ElectrideFlow Workflow**:
```bibtex
@article{electride_flow_2026,
  title={High-Throughput Discovery of Inorganic Electrides},
  author={...},
  journal={arXiv preprint arXiv:2601.21077},
  year={2026},
  url={https://arxiv.org/abs/2601.21077}
}
```

**Database**:
```bibtex
@misc{electride_database_2026,
  title={Inorganic Electride Database},
  author={Materials Modelling and Informatics Group},
  year={2026},
  url={https://inorganic-electride.onrender.com}
}
```

### Related Links

- **ElectrideFlow Repository**: https://github.com/MaterSim/ElectrideFlow
- **Research Group**: https://qzhu2017.github.io

---

## License

MIT License - See LICENSE.md for details.

---

## Contact

- **Email**: qzhu8@charlotte.edu
- **Issues**: Open an issue on GitHub for bugs or feature requests

---

**Last updated**: January 2026
