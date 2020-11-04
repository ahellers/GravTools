"""Models for gravity surveys, independent of the used gravity meter."""


class Survey:
    """Gravity survey (instrument-independent)"""

    def __init__(self):
        self.name = ''

    def __str__(self):
        return self.name
