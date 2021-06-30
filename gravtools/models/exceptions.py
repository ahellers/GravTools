"""
gravtools
=========

Code by Andreas Hellerschmied
andeas.hellerschmid@bev.gv.at

Summary
-------
Contains user defined exceptions.
"""


class Error(Exception):
    """Base class for exceptions."""
    pass


class InvaliFileContentError(Error):
    """Exception raised for invalid file content.

    Parameters
    -----------
    message : str
        Explanation of the error
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class FileTypeError(Error):
    """Exception raised for invalid file types.

    Parameters
    -----------
    message : str
        Explanation of the error
    valid_file_types : list of str
        List of strings declaring the allowed file types.
    """

    def __init__(self, message, valid_file_types=None):
        self.message = message
        if valid_file_types is not None:
            append_str = '\nThe following file types are valid:\n'
            append_str = append_str + '\n'.join(valid_file_types)
            self.message = self.message + append_str
        super().__init__(self.message)
