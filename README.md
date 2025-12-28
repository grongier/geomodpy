# geomodpy

geomodpy is a Python package for geostatistics and geomodeling. Its core component is a [wrapper](geomodpy/wrapper) to run GSLIB-like functions compiled in separate executables, including [GSLIB](http://www.gslib.com/) itself.

## Installation

You can directly install geomodpy from GitHub using pip:
   ```
   pip install git+https://github.com/grongier/geomodpy.git
   ```

This includes precompiled Fortran executables but, if you need to compile the Fortran sources yourself, you can follow those steps on Debian-based Linux:

1. Open a terminal and update the package list:  
   ```
   sudo apt-get update
   ```
2. Install *make* and *gfortran*:  
   ```
   sudo apt-get -y install make gfortran
   ```
3. Go to geomodpy:  
   ```
   cd whereisgeomody/geomodpy
   ```
4. Compile:  
   ```
   make
   ```

Those steps on MacOS:

1. Open a terminal and install Xcode’s Command Line Tools:  
   ```
   xcode-select --install
   ```
2. Install Homebrew:  
   ```
   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
   ```
3. Install *make* and *gfortran*:  
   ```
   brew install make gcc
   ```
4. Go to geomodpy:  
   ```
   cd whereisgeomody/geomodpy
   ```
5. Compile:  
   ```
   make
   ```

And those steps on Windows:

1. Install MSYS2: [www.msys2.org](https://www.msys2.org/)  

2. Open *MSYS2 MINGW64*  

3. Update the installation:  
   ```
   pacman -Syu
   ```
4. Install the MinGW-w64 toolchain:
   ```
   pacman -S mingw-w64-x86_64-toolchain
   ```
5. Move to the C drive:
   ```
   cd c:/
   ```
6. Go to geomodpy:  
   ```
   cd whereisgeomody/geomodpy
   ```
7. Compile:  
   ```
   mingw32-make
   ```

## Usage

Here's a basic example of geomodpy's API:

```
from geomodpy.wrapper.gslib import SGSIM

# Instanciate a Gaussian simulation
sgsim = SGSIM()
# Run the simulation
realization = sgsim.run()
# Plot the realization
realization.isel(Realization=0, W=0)['Variable'].plot(x='X', y='Y')
```

For a more complete example, see the Jupyter notebook [using_geomodpy.ipynb](examples/using_geomodpy.ipynb) in [examples](examples).

## Citation

If you use geomodpy in your research, please cite the original article(s) describing the method(s) you used (see the docstrings for the references). An acknowledgment of geomodpy's use is always appreciated.

## Credits

This software was written by:

| [Guillaume Rongier](https://github.com/grongier) <br>[![ORCID Badge](https://img.shields.io/badge/ORCID-A6CE39?logo=orcid&logoColor=fff&style=flat-square)](https://orcid.org/0000-0002-5910-6868)</br> |
| :---: |

## License

Copyright notice: Technische Universiteit Delft hereby disclaims all copyright interest in the program geomodpy written by the Author(s). Prof.dr.ir. S.G.J. Aarninkhof, Dean of the Faculty of Civil Engineering and Geosciences

&#169; 2022-2025, Guillaume Rongier

This work is licensed under a MIT OSS licence, see [LICENSE](LICENSE) for more information.
