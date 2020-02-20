from setuptools import setup

setup(
    name="flowplot",
    version="0.1dev",
    author="Max Jacobi",
    author_email="mjacobi@theorie.ikp.physik.tu-darmstadt.de",
    packages=['flowplot', 'flowplot/plots', 'flowplot/flow'],
    long_description=open('README.md').read(),
    description="Some scripts to plot flows from WinNet"
)
