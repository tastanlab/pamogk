from distutils.core import setup

# load requirements from file
import setuptools

with open('README.md') as f:
    long_description = f.read()

setup(
    name='pamogk',
    version='0.1.0',
    packages=setuptools.find_packages(),
    install_requires=requirements,
    license='Creative Commons Attribution-Noncommercial-Share Alike license',
    author="Fma",
    description="PAMOGK",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/tastanlab/pamogk",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
