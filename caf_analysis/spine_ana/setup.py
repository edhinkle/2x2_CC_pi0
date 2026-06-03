from setuptools import setup, find_packages

setup(
    name="spine_CC1pi0_analyzer_plotting",
    version="0.1.0",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "matplotlib>=3.5.0",
        "numpy>=1.23.0",
        "uproot>=4.0.0",
    ],
)