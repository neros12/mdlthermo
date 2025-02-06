from setuptools import setup, find_packages

setup(
    name="my_project",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas",
    ],
    description="A sample Python project",
    author="Your Name",
    license="MIT",
    url="https://github.com/yourusername/my_project",
)
