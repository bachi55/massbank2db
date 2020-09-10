from setuptools import setup, find_packages

setup(
    name="massbank2db",
    version="0.1.1",
    license="MIT",
    packages=find_packages(exclude=["tests", "examples", "*.ipynb", "inspect_candidates.py"]),

    # Minimum requirements the package was tested with
    install_requires=[
        "pandas",
        "numpy",
        "setuptools>=46.1"
    ],

    # Metadata
    author="Eric Bach",
    author_email="eric.bach@aalto.fi",
    description="Build a local SQLite Database from Massbank.",
    url="https://github.com/bachi55/massbank2db",
)