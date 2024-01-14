"""
Classes to represent a single point charge and a point charge system.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import constants

# FIXME: In order to handle the system of point charges, construct an equation with a term for the field of each point
#        charge. Do this through a for loop where each item is a term, then use np.sum to add them all up at each
#        location in the image. THE WAY THIS IS CURRENTLY IMPLEMENTED WILL NOT WORK. CHANGE TO THIS.


# NOTE: this will only work for a discrete system of charges, and not for a continuous distribution described
#       by some function... can we work in another way to incorporate continuous distributions? ie create a function and 
#       feed it into a class for dealing with this? Would have to find a way to integrate the function (maybe numpy or scipy
#       has something for this)


class PointCharge:
    """
    Class to represent a single point charge.
    """

    def __init__(self, charge: float or int, loc: list = [0, 0]):
        """
        :param charge: Charge of the point charge in Coulombs
        :param loc: Location [x,y] of the point charge
        """
        self.charge = charge
        self.loc = loc
        self.E = []

    def calculate_E_field(self, rmax: float = 1.5, step: float =0.01):
        """
        Calculate the strength of the electric charge
        :param rmax: Maximum distance from the point charge to calculate the field in meters
        :param location: Location of the point charge
        :param step: Step size between each x,y value (set to a smaller value for higher resolution)
        :return: 2D numpy array representing the field values produced by the point charge
        """
        # Define constants for Coulomb's Law: E = k*q/r^2
        eps_0 = constants.epsilon_0
        pi = constants.pi
        k = (4 * pi * eps_0)**-1

        # Find the x,y bounds for the specified radius (rmax)
        xmax = rmax / np.sqrt(2)
        xmin = -xmax
        ymax = rmax / np.sqrt(2)
        ymin = -ymax
        # Define pixel coordinate offset
        x_offset = (xmax - np.abs(xmin)) / 2 # TODO: incorporate this into the plot coords
        y_offset = (ymax - np.abs(ymin)) / 2
        # Calculate the E field for all values of (x,y) in the specified radius
        #for y in np.arange(ymin, ymax, step):
        for y in np.arange(self.loc[1] - ymax, self.loc[1] + ymax, step):
            yprime = y - self.loc[1]
            E_row_vals = []
            #for x in np.arange(xmin, xmax, step):
            for x in np.arange(self.loc[0] - xmax, self.loc[0] + xmax, step):
                xprime = x - self.loc[0]
                r_sq = xprime**2 + yprime**2
                try:
                    E = k * self.charge / (r_sq)
                except ZeroDivisionError:
                    E = 0
                E_row_vals.append(E)
            self.E.append(E_row_vals)
        # TODO: add functionality to place the charge at a location other than the origin -- will have to offset the loc from the origin
        return self.E, [xmin, xmax], [ymin, ymax]

    def plot(self, rmax: float = 1.5, cmap='plasma', figsize: list[int] = [6, 6]):
        """
        Produce a plot of the 2D field of the point charge
        :param cmap: Matplotlib colormap to use. See https://matplotlib.org/stable/gallery/color/colormap_reference.html
        :param figsize: [x, y] size of the figure to produce. Must be a 2-element list
        :return: None # TODO: return the figure obj. instead
        """
        _, xrange, yrange = self.calculate_E_field(rmax)
        plt.figure(figsize = figsize)
        # TODO: handle negative charge values
        plt.imshow(np.log(self.E), cmap=cmap, origin='lower')
        plt.colorbar(label = 'lnE [V/m]')
        plt.xlabel('X [m]', fontsize = 12)
        plt.ylabel('Y [m]', fontsize = 12)
        xticks = np.arange(xrange[0], xrange[1], 0.2) # TODO: fix these ticks to coincide with the actual primed coordinates !
        plt.xticks(np.linspace(0, len(self.E[0]), len(xticks)), [str(np.round(xtick, 2)) for xtick in xticks],
                   fontsize=8, rotation = 90)
        yticks = np.arange(yrange[0], yrange[1], 0.2)
        plt.yticks(np.linspace(0, len(self.E[0]), len(yticks)), [str(np.round(ytick, 2)) for ytick in yticks],
                   fontsize=8)
        plt.show()
        # TODO: we should return the figure obj. so that it can be further manipulated from the notebook


class PointChargeSystem(PointCharge):
    # TODO: will have to incorporate the locations into this...
    def __init__(self, chargesdict):
        """
        Class to handle the fields of a system of point charges.
        :param chargesdict: Dictionary of charges of the form: {'q1':[charge, [x,y]],...'qn':[charge,[x,y]]}
        Example: chargesdict = {'q1': [7, [1, 3]], 'q2': [-4, [2, 4]]}
        """
        super().__init__()
        self.chargesdict = chargesdict
        self.net_E = None

    def calculate_net_E_field(self) -> np.array:
        """
        Calculate the net electric field for the system of point charges in chargesdict.
        :return: 2D numpy array representing the net field values produced by the system of point charges
        """
        for q, loc in self.chargesdict.values():
            E = PointCharge(q, loc).calculate_E_field()
            if not self.net_E:
                self.net_E = np.shape(E)
            self.net_E = np.add(self.net_E, E)
            del E # free up memory as we go
        return self.net_E
