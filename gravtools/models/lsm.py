"""
gravtools
=========

Code by Andreas Hellerschmied
andeas.hellerschmid@bev.gv.at

Summary
-------
Contains classes for least-squares adjustment of gravimeter-surveys.
"""

import pandas as pd
import numpy as np
import datetime as dt
import os
import copy

from gravtools.settings import SURVEY_DATA_SOURCE_TYPES, STATION_DATA_SOURCE_TYPES, GRAVIMETER_ID_BEV, \
    TIDE_CORRECTION_TYPES, DEFAULT_GRAVIMETER_ID_CG5_SURVEY, REFERENCE_HEIGHT_TYPE, NAME_OBS_FILE_BEV, \
    PATH_OBS_FILE_BEV, BEV_GRAVIMETER_TIDE_CORR_LOOKUP, GRAVIMETER_REFERENCE_HEIGHT_CORRECTIONS_m
from gravtools.const import VG_DEFAULT
from gravtools.models.exceptions import FileTypeError
from gravtools.CG5_utils.cg5_survey import CG5Survey


class LSM_diff:
    """Least-squares adjustment of differential gravimeter observations with weighted constraints.

    Attributes
    ----------
    stat_df : :py:obj:`gravtools.Station.stat_df`, optional (default=None)
            The station dataframe contains all relevant station data.
        setups : dict, optional (default=None)
            The setups dictionary contains all observation data used for the adjustment. The keys of the dictionary
            are the survey names (str) and the items are pandas dataframes containing the observation data (see
            :py:obj:`gravtool.Survey.setup_df`)
    """

    def __init__(self, stat_df=None, setups=None):
        """
        Parameters
        ----------
        stat_df : :py:obj:`gravtools.Station.stat_df`, optional (default=None)
            The station dataframe contains all relevant station data.
        setups : dict, optional (default=None)
            The setups dictionary contains all observation data used for the adjustment. The keys of the dictionary
            are the survey names (str) and the items are pandas dataframes containing the observation data (see
            :py:obj:`gravtool.Survey.setup_df`)

        Notes
        -----
        In order to prevent altering the original input data (in case of assignment via reference), deep copies of the
        input data are stored in objects of this class.

        """
        # Create deep copies of the input data items:
        self.stat_df = stat_df.copy(deep=True)
        self.setups = copy.deepcopy(setups)
        pass

    @classmethod
    def from_campaign(cls, cg5_survey, keep_survey=True):
        """Constructor that generates and populates the survey object from a CG5Survey class object.

        Notes
        -----
        The observation epochs are represented by timezone aware datetime objects with TZ=<UTC>. The TZ is changed
        if necessary when loading data from any source.

        Parameters
        ----------
        keep_survey : bool, optional (default=True)
            If False, this survey is excluded from further processing.
        cg5_survey : :py:obj:`gravtools.CG5_utils.survey.CG5Survey`
            Objects of the class CG5Survey contain all data from CG-5 observation files.

        Returns
        -------
        :py:obj:`.Survey`
            Contains all information of s specific survey independent of the data source.
        """