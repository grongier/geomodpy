from setuptools import setup
from setuptools.dist import Distribution
from wheel.bdist_wheel import bdist_wheel


# This is a hacky way to force cibuildwheel to build platform-dependent wheels
# that include platform-dependent binaries


class BinaryDistribution(Distribution):
    def has_ext_modules(self):
        # Tell setuptools this depends on the platform
        return True

    def is_pure(self):
        # Tell setuptools this isn't a pure Python wheel
        return False


class PlatformWheel(bdist_wheel):
    def get_tag(self):
        # Force the Python tag to py3 and the ABI to none
        python, abi, plat = super().get_tag()
        return ("py3", "none", plat)


setup(
    distclass=BinaryDistribution,
    cmdclass={"bdist_wheel": PlatformWheel},
)
