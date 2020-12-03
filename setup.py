from setuptools import setup
import os

from setuptools.command.install import install                                      

class CustomInstall(install):                                                       
    def run(self):                                                                
        os.chdir("exoplasim")
        os.system("./configure.sh")
        os.chdir("..")
        os.system("echo $(pwd)/exoplasim>exoplasim/__init__.py")
        os.system("cat exoplasim/exoplasim.py>>exoplasim/__init__.py")                                                             
        install.run(self)

setup(
    name='exoplasim',
    version='2.0.0',
    packages=['exoplasim',],
    install_requires=["numpy","netCDF4"],
    include_package_data=True,
    author='Adiv Paradise',
    author_email='paradise.astro@gmail.com',
    license='GNU General Public License',
    url='https://github.com/alphaparrot/ExoPlaSim',
    description='Exoplanet GCM',
    long_description=open('README.txt').read(),
    cmdclass={"install":CustomInstall},
    )
