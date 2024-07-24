"""Modelling gravimeters.

Copyright (C) 2024  Andreas Hellerschmied <andreas.hellerschmied@bev.gv.at>

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

import os

import pandas
import pandas as pd
import json
import numpy as np

import gravtools.settings as settings


class Gravimeters:
    """Manages gravity meters.

    Attributes
    ----------
    gravimeters : dict, optional (default = None)
        The keys are tuples of the form (<gravimeter type>, <gravimeter serial number>) and the values are `Gravimeter`
        objects describing the properties of a single gravity meter.

    """
    def __init__(self):
        """Default constructor."""
        self.gravimeters = {}

    def add_gravimeter(self, gravi, verbose=True):
        """Add a gravimeter object and carry out checks.

        Parameters
        ----------
        gravi: `Gravimeter` object
            Gravimeter object that contains all relevant data.
        verbose: bool, optional (default = `True`)
            Print terminal output if `True`.

        Notes
        -----
        Do not overwrite existing data loaded from a gravimeter file by default values and/or values from an observation
        file.
        """
        # Overwrite existing entry only, if the new data was loaded from a gravimeter file:
        if (gravi.gravimeter_type, gravi.serial_number) in self.gravimeters.keys():
            if gravi.data_source_type == 'survey' and self.gravimeters[(gravi.gravimeter_type, gravi.serial_number)].data_source_type == 'file':
                if verbose:
                    print(f'{gravi.name} was already loaded from a gravimeter file => Not added')
                return
        self.gravimeters[(gravi.gravimeter_type, gravi.serial_number)] = gravi
        if verbose:
            print(f'{gravi.name} added.')

    def add_from_json(self, filename_json: str, verbose=True):
        """Load data of one or multiple gravimeters a json file.

        Parameters
        ----------
        filename_json : str
            Name of the json file holding the gravimeter information.
        verbose: bool, optional (default = `True`)
            Print terminal output if `True`.

        Notes
        -----
        If the gravimeter already exists with data loaded from an observation file (= survey) and default scale factors
        (`self.data_source_type` is `survey`), this information is overwritten by the data loaded from the json file.
        """
        if verbose:
            print(f'Load gravimeter data from json file: {filename_json}')
        with open(filename_json, mode='r') as file:
            gravimeters = json.load(file)

        for gm in gravimeters:
            if verbose:
                print(f"Load gravimeter {gm['type'],} with S/N {gm['serial_number']}")
            gravi = Gravimeter(gravimeter_type=gm['type'],
                               manufacturer=gm['manufacturer'],
                               serial_number=gm['serial_number'],
                               height_offset_m=gm['height_offset_m'],
                               data_source=os.path.basename(filename_json),
                               data_source_type='file',
                               description=gm['description'],
                               scale_parameters=gm['calibration'],
                               code=gm['code'],
                               verbose=verbose)
            self.add_gravimeter(gravi)

    def add_from_survey(self, survey, verbose=True):
        """Get the gravimeter data from a single survey object and set unknown attributes to default values.

        Parameters
        ----------
        survey : `gravtools.models.survey.Survey` object
            GravTools survey object.
        verbose: bool, optional (default = `True`)
            Print terminal output if `True`.

        Notes
        -----
        If there is already an entry for a gravimeter loaded from an external source, e.g. from a json file, the
        existing content will NOT be overwritten by default values! In this case `self.data_source_type` is set to
        `file`.

        The default linear scale factor is defined in `settings.DEFAULT_GRAVIMETER_LINEAR_SCALE_FACTOR`
        """
        gravi = Gravimeter(gravimeter_type=survey.gravimeter_type,
                           manufacturer=settings.DEFAULT_GRAVIMETER_MANUFACTURERS[survey.gravimeter_type],
                           serial_number=survey.gravimeter_serial_number,
                           height_offset_m=settings.DEFAULT_GRAVIMETER_REFERENCE_HEIGHT_OFFSET_M[survey.gravimeter_type],
                           data_source=survey.name,
                           data_source_type='survey',
                           description=settings.DEFAULT_GRAVIMETER_DESCRIPTION,
                           scale_parameters=None,  # Use default values!
                           code=settings.DEFAULT_GRAVIMETER_ONE_LETTER_CODES[survey.gravimeter_type],
                           verbose=verbose)
        self.add_gravimeter(gravi)

    def delete_gravimeter(self, gravimeter_type: str, serial_number: str, verbose=True):
        """Delete the gravimeter with the given type and S/N.

        Parameters
        ----------
        gravimeter_type : str
            Gravimeter type.
        serial_number : str
            Instrument serial number
        verbose: bool, optional (default = `True`)
            Print terminal output if `True`.
        """
        del self.gravimeters[(gravimeter_type, serial_number)]
        if verbose:
            print(f'Deleted Gravimeter {gravimeter_type} with S/N {serial_number}.')

    def apply_linear_scale(self, gravimeter_type: str, serial_number: str, gravity_df: pd.DataFrame, verbose=False) -> pd.DataFrame:
        """Apply linear scale factor on observations of the gravimeter with the given type and S/N.

        Parameters
        ----------
        gravimeter_type : str
            Gravimeter type.
        serial_number : str
            Instrument serial number
        gravity_df : pd.DataFrame
            Pandas DataFrame with two columns containing the reference epochs (col. "epoch_dt") as datetime objects and
            the gravity values (col. "g") that are scaled with the linear scaling factor corresponding to the reference
            epochs.
        verbose: bool, optional (default = `False`)
            Print terminal output if `True`.

        Returns
        -------
        pd.DataFrame : "gravity_df" with two additional columns for the scale factors ("linear_scale") and the scaled
        gravity ("g_scaled").
        """
        if verbose:
            print(f'Linear scaling correction for gravimeter {gravimeter_type} ({serial_number})')
        return self.gravimeters[(gravimeter_type, serial_number)].apply_linear_scaling(gravity_df)

    def get_height_offset(self, gravimeter_type: str, serial_number: str):
        """Returns the height offset between the sensor and the top.

        Parameters
        ----------
        gravimeter_type : str
            Gravimeter type.
        serial_number : str
            Instrument serial number
        """
        return self.gravimeters[(gravimeter_type, serial_number)].height_offset_m

    @property
    def number_of_gravimeters(self):
        """Return the number of gravimeters."""
        return len(self.gravimeters)

    def __str__(self):
        return f'{self.number_of_gravimeters} gravimeters with known properties.'

    def __repr__(self):
        return self.__str__()


class Gravimeter:
    """Describes properties of a single gravity meter.

    Attributes
    ----------
    gravimeter_type : str
        Gravimeter type. The type has to be a key in `settings.GRAVIMETER_TYPES`.
    manufacturer : str
        Manufacturer of the instrument.
    serial_number : str
        Instrument serial number.
    height_offset_m : float
        Height offset between the instrument's reference surface (usually instrument top surface) and the sensor level.
    description : str
        Arbitrary description.
    data_source : str
        Source from which the gravimeter data was loaded. This is either the name of a file, or the name of a survey.
    data_source_type : str
        Type of the gravimeter data source. HAs to be listed as key in `settings.GRAVIMETER_DATA_SOURCE_TYPES`.
    scale_df : :py:obj:`pandas.core.frame.DataFrame`
        Pandas Dataframe containing the instrument's calibration parameters valid for time specific time intervals.

        - start_date : :py:obj:`datetime.date`
            Start of the time span at which the scale is valid.
        - end_date : :py:obj:`datetime.date`
            End of the time span at which the scale is valid.
        - linear_factor : float,
            Linear scaling factor valid for the given time span.
        - comment : str
            Arbitrary comment on the scale factor entry, e.g. to document the source or any issues.

    code : str, optional (default = '')
        Short instrument code.

    Notes
    -----
    All epochs have to be given w.r.t. UTC. Dates have bo TZ.
    """

    def __init__(self, gravimeter_type: str, manufacturer: str, serial_number: str, height_offset_m: float,
                 data_source: str, data_source_type: str, description='', scale_parameters=None, code='', verbose=True):
        """Default constructor.

        Parameters
        ----------
        gravimeter_type : str
            Gravimeter type. The type has to be a key in `settings.GRAVIMETER_TYPES`.
        manufacturer: str
            Manufacturer of the instrument.
        serial_number: str
            Instrument serial number.
        height_offset_m: float
            Height offset between the instrument's reference surface (usually instrument top surface) and the sensor
            level.
        data_source : str
            Source from which the gravimeter data was loaded. This is either the name of a file, or the name of a
            survey.
        data_source_type : str
            Type of the gravimeter data source. Has to be listed as key in `settings.GRAVIMETER_DATA_SOURCE_TYPES`.
        description : str, optional (default = '')
            Arbitrary description.
        scale_parameters : list of dicts, optional (default = `None`)
            Each dict contains the start date, the end date and the linear calibration factor. If no data is provided,
            default values are used.
        code: str, optional (default = '')
            Short instrument code.
        verbose: bool, optional (default = `True`)
            Print terminal output if `True`.
        """
        flag_init_scale_factor_zero = False

        # gravimeter_type
        if not isinstance(gravimeter_type, str):
            raise TypeError('"gravimeter_type" has to by a string.')

        # manufacturer
        if not isinstance(manufacturer, str):
            raise TypeError('"manufacturer" has to by a string.')

        # serial_number
        if not isinstance(serial_number, str):
            raise TypeError('"serial_number" has to by a string.')

        # height_offset_m
        if not isinstance(height_offset_m, float):
            raise TypeError('"height_offset_m" has to by a float.')

        # data_source
        if not isinstance(data_source, str):
            raise TypeError('"data_source" has to by a string.')

        # data_source_type
        if not isinstance(data_source_type, str):
            raise TypeError('"data_source_type" has to by a string.')
        if data_source_type not in settings.GRAVIMETER_DATA_SOURCE_TYPES.keys():
            raise RuntimeError(f'"{data_source_type}" not listed settings.GRAVIMETER_DATA_SOURCE_TYPES!')

        # description
        if not isinstance(description, str):
            raise TypeError('"description" has to by a string.')

        # calibrations
        if scale_parameters is None:
            flag_init_scale_factor_zero = True
        else:
            if isinstance(scale_parameters, list):
                if len(scale_parameters) > 0:
                    # Save calibration factors to dataframe
                    scale_dict = {'start_date': [], 'end_date': [], 'linear_factor': [], 'comment': []}
                    for scale in scale_parameters:
                        scale_dict['start_date'].append(scale['start_date'])
                        scale_dict['end_date'].append(scale['end_date'])
                        scale_dict['linear_factor'].append(scale['linear_factor'])
                        scale_dict['comment'].append(scale['comment'])
                        scale_df = pd.DataFrame(scale_dict)
                else:
                    flag_init_scale_factor_zero = True
            else:
                raise TypeError('"code" has to by a string.')
        if flag_init_scale_factor_zero:
            if verbose:
                print(f'Initialize a linear scale factor of 1.0 for gravimeter {manufacturer} '
                      f'{gravimeter_type}')
            scale_dict = {'start_date': [settings.DEFAULT_CALIBRATION_START_DATE],
                          'end_date': [settings.DEFAULT_CALIBRATION_END_DATE],
                          'linear_factor': [settings.DEFAULT_CALIBRATION_LINEAR_FACTOR],
                          'comment': [settings.DEFAULT_CALIBRATION_COMMENT]}
            scale_df = pd.DataFrame(scale_dict)
        scale_df['start_date'] = pd.to_datetime(scale_df['start_date'], utc=True).dt.date
        scale_df['end_date'] = pd.to_datetime(scale_df['end_date'], utc=True).dt.date

        # Check for overlapping time intervals:
        if len(scale_df) > 1:
            for index, row in scale_df.iterrows():
                tmp_filter = (((scale_df['end_date'] <= row['end_date']) & (scale_df['end_date'] >= row['start_date'])) | ((scale_df['start_date'] <= row['end_date']) & (scale_df['start_date'] >= row['start_date'])))
                num_matches = len(tmp_filter[tmp_filter])
                if num_matches > 1:
                    tmp_df = scale_df.loc[tmp_filter]
                    error_str = tmp_df.to_string(index=False)
                    raise RuntimeError(f'Overlapping calibration factor intervals of gravimeter {gravimeter_type} '
                                       f'(S/N {serial_number}):\n' + error_str)

        # code
        if not isinstance(code, str):
            raise TypeError('"code" has to by a string.')

        # Assign and save data. if no errors occurred:
        self.gravimeter_type = gravimeter_type
        self.manufacturer = manufacturer
        self.serial_number = serial_number
        self.height_offset_m = height_offset_m
        self.data_source = data_source
        self.description = description
        self.data_source_type = data_source_type
        self.scale_df = scale_df
        self.code = code

    @property
    def num_scale_factors(self) -> int:
        """Return the number of scale factors."""
        return len(self.scale_df)

    def get_linear_scale_factor(self, epoch: pandas.Timestamp, raise_error_availability: bool = False,
                                verbose=True) -> float:
        """Returns the linear scale factor for the given epoch.

        Parameters
        ----------
        epoch : Pandas.Timestamp (timezone aware with tz=UTC).
            Epoch for which a scale factor should be returned
        raise_error_availability : optional (default = `False`)
            `True` implies that an error is raised if no scale factor is available for the given epoch. Otherwise, a
            scale factor equal `np.nan` is returned without raising an error.
        verbose: bool, optional (default = `True`)
            Print terminal output if `True`.

        Returns
        -------
        float : Scale factor or `numpy.nan` if not available or in case of multiple matches.

        Notes
        -----
        Only the dates are compared when matching epochs!

        """
        tmp_filter = (self.scale_df['start_date'] <= epoch.date()) & (self.scale_df['end_date'] >= epoch.date())
        num_matches = len(tmp_filter[tmp_filter])
        if num_matches == 0:
            if raise_error_availability:
                raise RuntimeError(f'No linear scale factor available for epoch {epoch.isoformat()} and gravimeter {self.name}.')
            return np.nan
        elif num_matches == 1:
            return self.scale_df.loc[tmp_filter, 'linear_factor'].item()
        else:
            if verbose:
                print(f'Warning: {num_matches} matches for epoch {epoch.isoformat()} in the list of scale factors for '
                      f'the gravimeter {self.name}')
            if raise_error_availability:
                raise RuntimeError(f'Warning: {num_matches} matches for epoch {epoch.isoformat()} in the list of scale '
                                   f'factors for the gravimeter {self.name}')
            return np.nan

    def get_linear_scale_factors(self, ref_epochs: list, raise_error_availability: bool = False) -> list[float]:
        """Returns the linear scale factor for the given epoch.

        Parameters
        ----------
        ref_epochs : List of datetime64 objects (timezone aware with tz=UTC).
            Epochs for which a scale factor should be returned
        raise_error_availability : optional (default = `False`)
            `True` implies that an error is raised if no scale factor is available for the given epoch. Otherwise, a
            scale factor equal `np.nan` is returned without raising an error.

        Returns
        -------
        list : Scale factors.

        """
        linear_scale_factors: list[float] = []
        for ref_epoch in ref_epochs:
            linear_scale_factors.append(self.get_linear_scale_factor(ref_epoch, raise_error_availability))
        return linear_scale_factors

    def apply_linear_scaling(self, gravity_df: pd.DataFrame) -> pd.DataFrame:
        """Applies a linear scaling factor on the input gravity data.

        Parameters
        ----------
        gravity_df : pd.DataFrame
            Pandas DataFrame with two columns containing the reference epochs (col. "epoch_dt") as datetime objects and
            the gravity values (col. "g") that are scaled with the linear scaling factor corresponding to the reference
            epochs.

        Notes
        -----
        An error is raised if noe scaling factor is available for the given epoch.

        Returns
        -------
        pd.DataFrame : "gravity_df" with two additional columns for the scale factors ("linear_scale") and the scaled
        gravity ("g_scaled").
        """
        gravity_df['linear_scale'] = self.get_linear_scale_factors(gravity_df['epoch_dt'].to_list(),
                                                                   raise_error_availability=True)
        gravity_df['g_scaled'] = gravity_df['linear_scale'] * gravity_df['g']
        return gravity_df

    @property
    def name(self):
        """Return the name of the gravimeter by combing type and S/N"""
        return self.gravimeter_type + f' ({self.serial_number})'

    def __str__(self):
        return (f"Gravimeter type {self.gravimeter_type} by {self.manufacturer} has {self.num_scale_factors} "
                f"scale factors.")

    def __repr__(self):
        return self.__str__()


if __name__ == '__main__':
    """Main function for testing"""
    filename = '../../data/gravimeters/gravimeters_template.json'
    meters = Gravimeters()
    meters.add_from_json(filename)
    print(meters)

    epochs = pd.date_range(start='2019-06-23', end='2022-08-15', periods=10, tz='UTC').to_list()
    epoch_pd = pd.to_datetime('2023-12-31 12:33:59', utc=True)
    grav = meters.gravimeters[('CG5', '40601')]

    scale = grav.get_linear_scale_factor(epoch_pd)

    scales = grav.get_linear_scale_factors(epochs)

    print('end')
