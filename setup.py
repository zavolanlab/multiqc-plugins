from setuptools import setup, find_packages

setup(
    name="multiqc_plugins",
    version="0.0.1",
    author="Krish Agarwal",
    author_email="akrish136@gmail.com",
    description="MultiQC plugins for the Zavolan Lab \
        @ University of Basel, Switzerland",
    long_description=__doc__,
    keywords="bioinformatics",
    url="https://github.com/zavolanlab/multiqc-plugins",
    download_url="",  # after releasing a version
    license="Apache 2.0",
    packages=find_packages(),
    include_package_data=True,
    install_requires=["multiqc"],
    entry_points={
        "multiqc.modules.v1": [
            "ALFA = modules.ALFA:MultiqcModule",
        ],
        "multiqc.hooks.v1": ["execution_start = modules.ALFA:ALFA_execution_start"],
    },
    classifiers=[
        "Development Status :: Beta",
        "Environment :: Console",
        "Environment :: Web Environment",
        "Intended Audience :: Science/Research",
        "License :: Apache 2.0",
        "Natural Language :: English",
        "Operating System :: Unix",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
)
