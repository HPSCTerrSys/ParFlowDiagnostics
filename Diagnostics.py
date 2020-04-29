# <editor-fold imports
import heat as ht
import numpy as np
import sys
import os


def printroot(*args, **kwargs):
    """ Only Root Process prints """
    if ht.MPI_WORLD.rank == 0:
        print(*args, **kwargs)
# </editor-fold>


class Diagnostics():  # Make this a subclass of ht.DNDarray?
    def __init__(self, Press, Satur, Mask, Permx, Permy, Permz,
                 Alpha, N, Sres, Poro, Sstorage, Mannings,
                 Slopex, Slopey, Sourcesink, Returnvalue=None):
        self.Press = Press
        self.Satur = Satur
        self.Mask = Mask
        self.Permx = Permx
        self.Permy = Permy
        self.Permz = Permz
        self.Alpha = Alpha
        self.N = N
        self.Sres = Sres
        self.Poro = Poro
        self.Sstorage = Sstorage
        self.Mannings = Mannings
        self.Slopex = Slopex
        self.Slopey = Slopey
        self.Sourcesink = Sourcesink
        self.Returnvalue = Returnvalue

    def TotalSubsurfaceStorage(self):
        pass


if __name__ == '__main__':
    pass
