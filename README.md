# MM_SSS_Lesson_Plan
This is the MolSSI Software Summer School's MM Lesson Plan. This is a private
repository so our (very clever) students do not find and copy at will instead
of writing their own.

This MM module will include:
 - LJ energy and gradient potenital
 - MC Integrator
 - MD Verlet Integrators
 - Pressure Probe
 - Temperature Probe
 - Radial Distribution Functions/plotting tools
 - Velocity-Autocorrelation functions
 - Diffusion calculator
 - Viscosity calculator

## Testing

To run the test suite please first run `pip install -e .` in the base
repository folder. This will register this repository with your local Python so
that `import mm_python` will work in any directory. Tests can then be run
with `py.test -v`. If `pytest` is not found please use `pip install pytest` to
aquire the module.
