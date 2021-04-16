from setuptools import setup, find_packages, Extension
from distutils.util import convert_path

mzClust_hclust = Extension('massbank2db.mzClust_hclust', sources=["src/mzClust_hclust.c"])

main_ns = {}
ver_path = convert_path('massbank2db/version.py')
with open(ver_path) as ver_file:
    exec(ver_file.read(), main_ns)

setup(
    name="massbank2db",
    version=main_ns["__version__"],
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
