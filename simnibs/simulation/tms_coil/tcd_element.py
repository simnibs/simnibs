class TcdElement:
    """Represents a part of a tcd file
    """
    def to_tcd(self) -> dict:
        """Turns the element into a tcd like dictionary

        Returns
        -------
        dict
            The tcd like dictionary representing the element
        """
        pass

    @classmethod
    def from_tcd_dict(cls, tcd_element: dict):
        """Loads the element from a tcd like dictionary

        Parameters
        ----------
        tcd_element : dict
            The element that was loaded from the tcd like dictionary
        """
        pass
