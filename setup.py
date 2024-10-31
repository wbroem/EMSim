from setuptools import setup

REQUIREMENTS = [
    'numpy',
    'scipy',
    'matplotlib',
    'argparse',
]

# TODO: create entry points
setup(entry_points = {
    'efield = efield.main' # edit this to where we put main() for an efield
})