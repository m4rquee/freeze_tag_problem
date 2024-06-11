from sys import argv

from src.cp.utils.instance import *


def read_arguments(*types):
    for t, arg in zip(types, argv[1:]): yield t(arg)
