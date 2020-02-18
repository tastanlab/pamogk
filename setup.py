from distutils.core import setup
import setuptools

# load requirements from file
with open('requirements.txt') as f:
    requirements = f.read().splitlines()

with open('README.md') as f:
    long_description = f.read()

setup(
    name='pamogk',
    version='0.1.0',
    packages=['pamogk',],
    install_requires=requirements,
    license='Creative Commons Attribution-Noncommercial-Share Alike license',
    long_description=long_description,
)
