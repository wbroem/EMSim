"""
Classes to represent the magnetic field produced by a steady current (static B field)
and a changing current (time-dependent B field).
"""

class Current:
    def __init__(self, charge, speed):
        self.charge = charge
        self.speed = speed

    def calculate_B_field(self):
        pass

    # TODO: so much...


class TimeDependentCurrent:
    def __init__(self):
        pass
    # TODO: How do we want to implement this?