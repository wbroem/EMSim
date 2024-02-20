"""
Classes to represent a single point charge and a point charge system.
"""

import json
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from scipy import constants


class PointCharge:
    """
    Class to represent a single point charge.
    """

    def __init__(self, charge: float | int, loc: list = [0, 0]):
        """
        :param charge: Charge of the point charge in Coulombs
        :param loc: Location [x,y] of the point charge
        """
        self.charge = charge
        self.loc = loc
        self.E = []

    def _find_coord_extrema(self, rmax: float = 1.5):
        """
        Find the x,y bounds for the specified radius (rmax)
        """

        self.xmax = self.loc + rmax
        self.xmin = self.loc - rmax
        self.ymax = self.loc + rmax
        self.ymin = self.loc - rmax

        return self.xmax, self.xmin, self.ymax, self.ymin
    
    @staticmethod # TODO: Should this be a static or class method? 
    def coulombs_law(charge, r_sq):
        """
        Utility function to calculate Coulomb's Law: E = k*q/r^2
        for a given charge and distance squared from that charge.
        :param charge: Charge of the point charge in Coulombs
        :param r_sq: Distance squared from the point charge in m^2
        """

        # Define constants
        eps_0 = constants.epsilon_0
        pi = constants.pi
        k = (4 * pi * eps_0)**-1
        
        try:
            E_loc_val = k * charge / (r_sq)
        except ZeroDivisionError:
            # E -> oo at the point charge location
            E_loc_val = np.inf
        return E_loc_val

    def calculate_E_field(self, step: float = 0.01):
        """
        Calculate the strength of the electric charge
        :param rmax: Maximum distance from the point charge to calculate the field in meters
        :param location: Location of the point charge
        :param step: Step size between each x,y value (set to a smaller value for higher resolution)
        :return: 2D numpy array representing the field values produced by the point charge
        """

        # Calculate the E field for all values of (x,y) in the specified radius
        for y in np.arange(self.ymax, self.ymin, -step):
            _Y = y - self.loc[1]
            E_row_vals = [] #np.array
            for x in np.arange(self.xmax, self.xmin, -step):
                _X = x - self.loc[0]
                r_sq = _X**2 + _Y**2
                E_loc_val = self.coulombs_law(self.charge, r_sq)
                E_row_vals.append(E_loc_val) #np.concatenate
            self.E.append(E_row_vals) #np.concatenate
        self.E = np.array(self.E) #remove this line
        return self.E

    def plot(self, rmax: float = 1.5, step: float = 0.01, cmap='plasma', figsize: list[int] = [6, 6]):
        """
        Produce a plot of the 2D field of the point charge
        :param cmap: Matplotlib colormap to use. See https://matplotlib.org/stable/gallery/color/colormap_reference.html
        :param figsize: [x, y] size of the figure to produce. Must be a 2-element list
        :return: None
        """
        if not self.E:
            self.calculate_E_field(rmax, step)
        plt.figure(figsize = figsize)
        # TODO: handle negative charge values
        plt.imshow(np.log(self.E), cmap=cmap, origin='lower')
        plt.colorbar(label = 'lnE [V/m]')
        plt.xlabel('X [m]', fontsize = 12)
        plt.ylabel('Y [m]', fontsize = 12)
        xticks = np.arange(self.xmin, self.xmax, rmax / 10)
        plt.xticks(np.linspace(0, len(self.E[0]), len(xticks)), [str(np.round(xtick, 3)) for xtick in xticks],
                   fontsize=8, rotation = 90)
        yticks = np.arange(self.ymin, self.ymax, rmax / 10)
        plt.yticks(np.linspace(0, len(self.E[0]), len(yticks)), [str(np.round(ytick, 3)) for ytick in yticks],
                   fontsize=8)



class PointChargeSystem(PointCharge):
    def __init__(self, charge_config: dict | str | Path = {}):
        """
        Class to handle the fields of a system of point charges.
        :param charge_config: Dictionary of charges of the form: 
            {"q1": {"charge": q1, "loc": [x1, y1]}, "q2": {"charge": q2, "loc": [x2, y2]}
            Alternatively, a .json file following this format can be loaded via the 
            load_charge_config() function.
        """
        super().__init__()
        self.charge_config = charge_config
        self.net_E = None
        self.rmax = None

    def load_charge_config(self, config_path) -> dict:
        """
        Load a .json file to specify the charge configuration for this object.
        :param config_path: Path to the .json file to be loaded.
        """
        config_path = Path(config_path)
        if not config_path.exists():
            raise FileNotFoundError(f'The specified path: {config_path} does not exist.')
        ext = config_path.suffix
        if not ext == '.json':
            raise TypeError('The configuration file should be a .json file, not {ext}')
        
        self.charge_config = json.load(config_path)
        return self.charge_config

    def _find_coord_extrema(self, rmax: float = 1.5):
        """
        Determine the dimensions of the E field image; ie find (xmin, xmax, ymin, ymax)
        considering the position and rmax for all of the point charges in the system.
        """
        self.rmax = rmax
        xmax_list = []
        xmin_list = []
        ymax_list = []
        ymin_list = []
        for charge in self.charge_config:
            loc = charge['loc']
            xmax = loc[0] + rmax
            xmax_list.append(xmax)
            xmin = loc[0] - rmax
            xmin_list.append(xmin)
            ymax = loc[1] + rmax
            ymax_list.append(ymax)
            ymin = loc[1] - rmax
            ymin_list.append(ymin)
        self.xmax = max(xmax_list)
        self.xmin = min(xmin_list)
        self.ymax = max(ymax_list)
        self.ymin = min(ymin_list)

        return self.xmax, self.xmin, self.ymax, self.ymin

    def calculate_net_E_field(self, step = 0.01) -> np.array:
        """
        Calculate the net electric field for the system of point charges in chargesdict.
        :return: 2D numpy array representing the net field values produced by the system of point charges
        """
        # FIXME: In order to handle the system of point charges, construct an equation with a term for the field of each point
#        charge. Do this through a for loop where each item is a term, then use np.sum to add them all up at each
#        location in the image. THE WAY THIS IS CURRENTLY IMPLEMENTED WILL NOT WORK.
        self._find_coord_extrema()
        for q, loc in self.charge_config.values(): # change how we read this to match the new config format
            # TODO: need to create E matrix such that it incorporates rmax from each point charge
            E = PointCharge(q, loc).calculate_E_field(step) # TODO: can we use super() instead?
            if not self.net_E:
                self.net_E = np.shape(E)
            self.net_E = np.add(self.net_E, E)
            del E # free up memory as we go
        return self.net_E

    def plot(self, rmax: float = 1.5, step: float = 0.01, cmap='plasma', figsize: list[int] = [6, 6]):
        super().plot()
        # TODO: whatever plot changes we need for this class