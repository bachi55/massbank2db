from massbank2db.version import __version__

from os.path import dirname as __os_dirname__
from os.path import join as __os_join__
from ctypes import cdll
from platform import python_version_tuple

__python__version__ = python_version_tuple()
HCLUST_LIB = cdll.LoadLibrary(
    __os_join__(__os_dirname__(__file__),
                "mzClust_hclust.cpython-%s%s-x86_64-linux-gnu.so" % (__python__version__[0], __python__version__[1])))

