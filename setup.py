from setuptools import setup
import os
 
setup(
    name='exoplasim-legacy',
    version='3.0.0a14',
    packages=['exoplasimlegacy',],
    zip_safe=False,
    install_requires=["numpy","matplotlib","scipy"],
    extras_require = {"netCDF4": ["netCDF4"],
                      "HDF5": ["h5py"]},
    include_package_data=True,
    author='Adiv Paradise',
    author_email='paradise.astro@gmail.com',
    license='GNU General Public License',
    license_files=["LICENSE.TXT",],
    url='https://github.com/alphaparrot/ExoPlaSim',
    description='Exoplanet GCM',
    long_description_content_type='text/x-rst',
    long_description=open('README.rst').read(),
    )
