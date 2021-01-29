from os.path import dirname as __os_dirname__
from os.path import join as __os_join__
from ctypes import cdll
HCLUST_LIB = cdll.LoadLibrary(__os_join__(__os_dirname__(__file__), "mzClust_hclust.cpython-39-x86_64-linux-gnu.so"))
