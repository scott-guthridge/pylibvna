from setuptools import Extension
from setuptools.command.build_py import build_py as _build_py
from numpy import get_include

class build_py(_build_py):
    def run(self):
        self.run_command("build_ext")
        return super().run()

    def initialize_options(self):
        super().initialize_options()
        if self.distribution.ext_modules is None:
            self.distribution.ext_modules = []

        self.distribution.ext_modules.append(
            Extension(name="vna.conv",
                      sources=["src/vna/conv.pyx", "src/vna/conv.pxd"],
                      define_macros=[('NPY_NO_DEPRECATED_API',
                                      'NPY_1_7_API_VERSION')],
                      include_dirs=[get_include(), "/usr/local/include"],
                      libraries=["vna", "m"])
        )
        self.distribution.ext_modules.append(
            Extension(name="vna.data",
                      sources=["src/vna/data.pyx", "src/vna/data.pxd"],
                      define_macros=[('NPY_NO_DEPRECATED_API',
                                      'NPY_1_7_API_VERSION')],
                      include_dirs=[get_include(), "/usr/local/include"],
                      libraries=["vna", "m"])
        )
        self.distribution.ext_modules.append(
            Extension(name="vna.cal",
                      sources=["src/vna/cal.pyx", "src/vna/cal.pxd"],
                      define_macros=[('NPY_NO_DEPRECATED_API',
                                      'NPY_1_7_API_VERSION')],
                      include_dirs=[get_include(), "/usr/local/include"],
                      libraries=["vna", "m"])
        )
