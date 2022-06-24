"""User defined exceptions.

Copyright (C) 2021  Andreas Hellerschmied <andreas.hellerschmied@bev.gv.at>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
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
