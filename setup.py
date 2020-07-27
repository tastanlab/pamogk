from distutils.core import setup
import toml

# load requirements from file
import setuptools

with open('README.md') as f:
    long_description = f.read()

with open('Pipfile') as f:
    pipfile = toml.load(f)
    required_packages = [f'{k}{"" if v == "*" else v}' for k, v in pipfile['packages'].items()]
    sep = '\n  '
    pkg_list = sep + sep.join(required_packages)
    print("Required packages:" + pkg_list)

print('Packages to export:', setuptools.find_packages())

setup(
    name='pamogk',
    version='0.1.4',
    packages=setuptools.find_packages(),
    install_requires=required_packages,
    license='Creative Commons Attribution-Noncommercial-Share Alike license',
    author="Fma",
    description="PAMOGK",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/tastanlab/pamogk",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6.1',
)
