# Calibration and Validation of TKTD and pathogen effect models

Organization of files in this project follows the defaults defined by the [project management package DrWatson.jl](https://github.com/JuliaDynamics/DrWatson.jl), specifically: 

- notebooks: Jupyter notebook files
- src: source code which is used across multiple projects and does not directly produce output (only funcion definitions). Functions defined here may be moved to separate packages later.
- scripts: Project-specific source code which may directly produce output. Most are called from the Jupyter notebooks.
- data: experimental data and simulation output
- plots: Main plots and figures. Some plots (e.g. calibration diagnostics) are stored in `data/sims/` together with the remaining output of the analysis (parameter values etc.)
- manuscript: associated manuscript files (typically .tex)