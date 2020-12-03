from setuptools import setup
import os

from setuptools.command.install import install                                       
    

class CustomInstall(install):                                                       
    def run(self):
        def _post_install():
            def find_module_path():
                for p in sys.path:
                    if os.path.isdir(p) and "exoplasim" in os.listdir(p):
                        return os.path.join(p, "exoplasim")
            install_path = find_module_path()
            cwd=os.getcwd()
            os.chdir(install_path)
            os.system("./configure.sh")
            os.system("echo %s>>__init__.py"%install_path)
            os.chdir(cwd)

        atexit.register(_post_install)
        install.run(self)

setup(
    name='exoplasim',
    version='2.0.0-1',
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
