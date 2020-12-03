from setuptools import setup
import os
 
setup(
    name='exoplasim',
    version='2.0.0.post5',
    packages=['exoplasim',],
    install_requires=["numpy","netCDF4"],
    include_package_data=True,
    author='Adiv Paradise',
    author_email='paradise.astro@gmail.com',
    license='GNU General Public License',
    url='https://github.com/alphaparrot/ExoPlaSim',
    description='Exoplanet GCM',
    long_description=open('README.txt').read(),
    )
