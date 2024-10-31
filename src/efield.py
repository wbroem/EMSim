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

    def __init__(self, charge: float | int = 0, loc: list = [0, 0]):
        """
        :param charge: Charge of the point charge in Coulombs
        :param loc: Location [x,y] of the point charge
        """
        self.charge = charge
        self.loc = loc
        self.E = None

    def _find_coord_extrema(self, rmax: float = 1.5):
        """
        Find the x,y bounds for the specified radius (rmax)
        """

        self.xmax = self.loc[0] + rmax
        self.xmin = self.loc[0] - rmax
        self.ymax = self.loc[1] + rmax
        self.ymin = self.loc[1] - rmax

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

    def calculate_E_field(self, charge, loc, step: float = 0.01):
        """
        Calculate the strength of the electric charge
        :param rmax: Maximum distance from the point charge to calculate the field in meters
        :param location: Location of the point charge
        :param step: Step size between each x,y value (set to a smaller value for higher resolution)
        :return: 2D numpy array representing the field values produced by the point charge
        """

        # Find the range of (x,y) values and create an empty matrix of matching size for E-field
        x_range = np.arange(self.xmax, self.xmin, -step)
        y_range = np.arange(self.ymax, self.ymin, -step)
        self.E = np.empty(shape=(len(y_range), len(x_range)))
        # Calculate the E field for all values of (x,y) in the specified radius
        row_number = 0
        for y in y_range:
            _Y = y - loc[1]
            # Start with an empty row
            E_row_vals = np.array([])
            # Calculate the E-field value for each x-value in a row (ie for a given y-value)
            for x in x_range:
                _X = x - loc[0]
                r_sq = _X**2 + _Y**2
                E_loc_val = self.coulombs_law(charge, r_sq)
                E_row_vals = np.append(E_row_vals, E_loc_val)
            # Place the row into the matrix
            self.E[row_number] = E_row_vals
            row_number += 1
        return self.E

    def log_format_E(self):
        """Create a natural log map of the electric field and handle negative or zero values."""
        # if not any values in E are negative:
        # return log(E)
        # else:
        # cast negative values to be -log(val) and positive values to be log(val)

        # Use np.maskedarray !
        neg_ix = np.where(self.E<0)

    def plot(self, rmax: float = 1.5, step: float = 0.01, cmap='plasma', figsize: list[int] = [6, 6]):
        """
        Produce a plot of the 2D field of the point charge
        :param cmap: Matplotlib colormap to use. See https://matplotlib.org/stable/gallery/color/colormap_reference.html
        :param figsize: [x, y] size of the figure to produce. Must be a 2-element list
        :return: None
        """
        if self.E is None:
            self._find_coord_extrema(rmax)
            self.calculate_E_field(self.charge, self.loc, step)
        plt.figure(figsize = figsize)
        #log_E_map = self.log_format_E()
        # plt.imshow(log_E_map, cmap=cmap, origin='lower')
        plt.imshow(self.E, cmap=cmap, origin='lower')
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

    def load_charge_config(self, config_path: str | Path) -> dict:
        """
        Load a .json file to specify the charge configuration for this object.
        :param config_path: Path to the .json file to be loaded.
        """
        if not isinstance(config_path, Path):
            config_path = Path(config_path)
        if not config_path.exists():
            raise FileNotFoundError(f'The specified path: {config_path} does not exist.')
        ext = config_path.suffix
        if not ext == '.json':
            raise TypeError('The configuration file should be a .json file, not {ext}')
        
        with open(config_path, 'r') as file:
            self.charge_config = json.load(file)
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
        for charge in self.charge_config.values():
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
        self._find_coord_extrema()
        if not isinstance(self.charge_config, dict):
            try:
                self.load_charge_config(self.charge_config)
            except Exception as e:
                raise Exception('Unable to load charge config. Use a dictionary or a valid .json config file. \n' + e)
        for charge in self.charge_config.values():
            E = super().calculate_E_field(charge['charge'], charge['loc'], step)
            if self.net_E is None:
                self.net_E = np.zeros(np.shape(E))
            self.net_E = np.add(self.net_E, E)
            del E # free up memory as we go
        return self.net_E

    def plot(self, rmax: float = 1.5, step: float = 0.01, cmap='plasma', figsize: list[int] = [6, 6]):
        if self.net_E is None:
            self._find_coord_extrema(rmax)
            self.calculate_net_E_field(step)
        plt.figure(figsize = figsize)
        # TODO: handle negative charge values
        plt.imshow(np.log(self.net_E), cmap=cmap, origin='lower')
        plt.colorbar(label = 'lnE [V/m]')
        plt.xlabel('X [m]', fontsize = 12)
        plt.ylabel('Y [m]', fontsize = 12)
        xticks = np.arange(self.xmin, self.xmax, (self.xmax - self.xmin) / 15)
        plt.xticks(np.linspace(0, len(self.E[0]), len(xticks)), [str(np.round(xtick, 3)) for xtick in xticks],
                   fontsize=8, rotation = 90)
        yticks = np.arange(self.ymin, self.ymax, (self.ymax - self.ymin) / 15)
        plt.yticks(np.linspace(0, len(self.E[0]), len(yticks)), [str(np.round(ytick, 3)) for ytick in yticks],
                   fontsize=8)
        # TODO: fix tick marks
