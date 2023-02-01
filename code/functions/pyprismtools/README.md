# pyprismtools

This python package implements functions that are used more than once in the analysis scripts of the META-PRISM study.
It typically defines helpful functions for defining colors in plots or for drawing specific types of plots.

## Installation

You may install the package using a distribution archive in the form of a tarball file or a wheel or you may install it
from the source folder in develop mode. Installing in develop mode allows to have changes to the source code taken into
account without having to reinstall the package each time it is modified.

In order to install in develop mode, run

```
pip install -e /path/to/pyprismtools
```

otherwise run

```
python3 -m pip install /path/to/local/archive
```

For this latter command, the distribution archive of the package must have been built beforehand. Please read the next
section to see how to build the package in case you need to do so.

**NOTE**: Be aware of the environment in which you are running the `pip` or `python3` command. If you want to have the
package installed in a specific environment, please activate it (using conda or any other environment manager) before
running the installation commands.q

## Build

For building a distribution archive, you will need the `build` package from **PyPi**. Install it with `python3 -m pip
install build` and build the package with

```
python3 -m build
```

## License

The pyprismtools package is licensed under the terms of the GPL Open Source license and is available for free.

