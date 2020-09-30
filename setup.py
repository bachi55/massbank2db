from setuptools import setup, find_packages, Extension

mzClust_hclust = Extension('massbank2db.mzClust_hclust', sources=["src/mzClust_hclust.c"])

setup(
    name="massbank2db",
    version="0.4.1",
    license="MIT",
    packages=find_packages(exclude=["tests", "examples", "*.ipynb", "inspect_candidates.py"]),

    # Minimum requirements the package was tested with
    install_requires=[
        "pandas",
        "numpy",
        "scipy",
        "setuptools>=46.1"
    ],

    # Metadata
    author="Eric Bach",
    author_email="eric.bach@aalto.fi",
    description="Build a local SQLite Database from Massbank.",
    url="https://github.com/bachi55/massbank2db",

    # Add hierarchical clustering implementation (in C) as dynamic library.
    ext_modules=[mzClust_hclust]
)
