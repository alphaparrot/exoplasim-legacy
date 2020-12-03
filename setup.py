from distutils.core import setup

os.system("./exoplasim/configure.sh")

os.system("echo $(pwd)/exoplasim>exoplasim/__init__.py")
os.system("cat exoplasim/exoplasim.py>>exoplasim/__init__.py")

setup(
    name='exoplasim',
    version='2.0.0',
    packages=['exoplasim',],
    author='Adiv Paradise',
    author_email='paradise.astro@gmail.com',
    license='GNU General Public License',
    url='https://github.com/alphaparrot/ExoPlaSim',
    description='Exoplanet GCM',
    long_description=open('README.txt').read(),
    )
