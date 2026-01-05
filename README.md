# geomodpy

geomodpy is a Python package for geostatistics and geomodeling. Its core component is a [wrapper](geomodpy/wrapper) to run GSLIB-like functions compiled in separate executables, including [GSLIB](http://www.gslib.com/) itself.

## Installation

### On Windows (64 bit)

You can directly install geomodpy from the wheel of the [release available on GitHub](https://github.com/grongier/geomodpy/releases) using pip:
   ```
   pip install https://github.com/grongier/geomodpy/releases/download/v0.0.1/geomodpy-0.0.1-py3-none-win_amd64.whl
   ```

The wheel already includes pre-compiled Fortran executables, but you can follow those steps to recompile them:

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

### On Linux (Debian-based)

You can directly install geomodpy from the wheel of the [release available on GitHub](https://github.com/grongier/geomodpy/releases) using pip:
   ```
   pip install https://github.com/grongier/geomodpy/releases/download/v0.0.1/geomodpy-0.0.1-py3-none-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl
   ```

The wheel already includes pre-compiled Fortran executables, but you can follow those steps to recompile them:

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

### On macOS

You can directly install geomodpy from the wheel of the [release available on GitHub](https://github.com/grongier/geomodpy/releases) using pip and the following command for Mac computers with Apple silicon:
   ```
   pip install https://github.com/grongier/geomodpy/releases/download/v0.0.1/geomodpy-0.0.1-py3-none-macosx_15_0_arm64.whl
   ```

And the following command for Mac computers with Intel processors:
   ```
   pip install https://github.com/grongier/geomodpy/releases/download/v0.0.1/geomodpy-0.0.1-py3-none-macosx_15_0_x86_64.whl
   ```

The latest versions of macOS might block the pre-compiled Fortran executables already included in the wheel. Running the following command might help:
   ```
   python -m geomodpy.fix_macos
   ```

Otherwise you can follow those steps to recompile them:

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
