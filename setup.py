from setuptools import setup, find_packages

setup(
    name="multiqc_plugins",
    version="1.3",
    author="Krish Agarwal",
    author_email="akrish136@gmail.com",
    description="MultiQC plugins for the Zavolan Lab \
        @ University of Basel, Switzerland",
    long_description=__doc__,
    keywords="bioinformatics",
    url="https://github.com/zavolanlab/multiqc-plugins",
    download_url="https://github.com/zavolanlab/multiqc-plugins",
    license="Apache 2.0",
    packages=find_packages(),
    include_package_data=True,
    install_requires=["multiqc"],
    entry_points={
        "multiqc.modules.v1": [
            "ALFA = modules.ALFA:MultiqcModule",
            "tin-score = modules.tin_score:MultiqcModule",
            "zpca = modules.zpca:MultiqcModule",
        ],
        "multiqc.hooks.v1": ["execution_start = modules.hook:execution_start"],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Environment :: Web Environment",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Natural Language :: English",
        "Operating System :: Unix",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
)
