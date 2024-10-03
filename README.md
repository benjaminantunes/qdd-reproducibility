# Performance and Reproducibility Assessment of Quantum Dissipative Dynamics Framework

### A Comparative Study of Fortran Compilers, MKL, and FFTW

[Read the full paper here](https://hal.science/hal-04684180/)

---

### Project Structure

- **`QDDpackage_withexamples/`**:  
  Contains the original code and compiling options from QDD.

- **`FFTW_repeat/`**:  
  Contains the code where FFTW is set to be repeatable between runs.

- **`tryRepro/`**:  
  Contains the code and compilation options where we attempted to obtain bitwise identical results between executions, for instance, when enabling or disabling debug mode.  
  - In the basic compilation options, when running the code with or without debug options, results were not repeatable.
