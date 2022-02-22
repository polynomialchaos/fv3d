################################################################################
# @file setup.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2022-02-22
# @copyright Copyright (c) 2022 by Florian Eigentler.
#  This work is licensed under terms of the MIT license (<LICENSE>).
################################################################################
from setuptools import setup, find_packages

version_locals = {}
with open('pymeshfv3d/version.py', 'r') as fp:
    exec(fp.read(), globals(), version_locals)

setup(
    name='pyMeshFV3D',
    version=str(version_locals['__version__']),
    author='Florian Eigentler',
    author_email='f.m.eigentler@gmail.com',
    description='Finite volume solver (FV3D) mesh preprocessing',
    packages=find_packages(),
    install_requires=['h5py'],
    entry_points={
        'console_scripts': [
            'pyMeshFV3D=pymeshfv3d.bin.pyMeshFV3D:main',
            'pyMeshVtk=pymeshfv3d.bin.pyMeshVtk:main',
            'pySolutionVtk=pymeshfv3d.bin.pySolutionVtk:main',
        ]
    }
)
