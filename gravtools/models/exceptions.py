class Error(Exception):
    """Base class for exceptions."""
    pass


class InvaliFileContentError(Error):
    """Exception raised for invalid file content.

    Attributes:
        expression -- expression where the error occurred
        message -- explaination of the error
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)
