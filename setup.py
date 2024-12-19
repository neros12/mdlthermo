# TEST VERSION

from setuptools import setup, find_packages

setup(
    name="MDL_modules",
    version="0.0.0",
    description="Integrated Computational Package for Thermodynamic Property Prediction",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="SunYoo Hwnag, BumChan Ryu",
    author_email="neros12@naver.com",
    url="https://github.com/neros12/MDL_modules",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=[
        "numpy",
        "pandas",
    ],
)
