####
#
# The 'massbank2db' package can be used to an build SQLite DB from the Massbank MS/MS repository.
#
#     Copyright (C) 2020 - 2022  Eric Bach <eric.bach@aalto.fi>
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
####

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
    packages=find_packages(exclude=["tests", "examples", "*.ipynb", "generate_massbank_sqlite.py"]),

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
