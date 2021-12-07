import datetime as dt
import numpy as np
import pickle
import os

from gravtools.models.lsm import LSM
from gravtools.models.lsm_diff import LSMDiff
from gravtools.models.lsm_nondiff import LSMNonDiff
from gravtools.models.mlr_bev_legacy import BEVLegacyProcessing
from gravtools.models.survey import Survey
from gravtools.models.station import Station
from gravtools.settings import ADDITIVE_CONST_ABS_GRTAVITY, GRAVIMETER_TYPES_KZG_LOOKUPTABLE, \
    GRAVIMETER_SERIAL_NUMBER_TO_ID_LOOKUPTABLE, GRAVIMETER_REFERENCE_HEIGHT_CORRECTIONS_m
from gravtools import __version__ as GRAVTOOLS_VERSION


class Campaign:
    """Gravity Campaign dataset.

    A gravity campaign datasets consists of:

    - Campaign name
    - One or more gravity surveys that belong together and
      - Each survey was observed with one gravimeter on a single day
    - Station data (datum and non-datum stations)
    - Reductions and corrections
      - All observations (from surveys) are corrected and reduced in the same way

    Attributes
    ----------
    campaign_name : str
        Name of the campaign.
    output_directory : str
        Path to output directory (all output files are stored there).
    surveys: dict of :py:obj:`.Survey` objects
        Arbitrary number of survey objects.
    stations : :py:obj:`.Station` object
        Data of known stations (datum- and non-datum-stations).
    lsm_runs : list of objects inherited from :py:obj:`gravtools.models.lsm.LSM`
        Each item in the list contains one enclosed LSM object. Each LSM object reflects one dedicated run of an
        least-squares adjustment in order to estimate target parameters.
    ref_delta_t_dt : datetime object
        Reference epoch for relative times within the campaign, e.g. for the determination of the drift polynomials.
        This reference time is equalt to the first (active) observation in the campaign considering all surveys.
    gravtools_version : str
        Version of the gravtools software that was used to create the dataset.
    """

    def __init__(self,
                 campaign_name,
                 output_directory,
                 surveys=None,  # Always use non-mutable default arguments!
                 stations=None,  # Always use non-mutable default arguments!
                 lsm_runs=None,  # Always use non-mutable default arguments!
                 ref_delta_t_dt=None  # Reference time for drift determination
                 ):
        """
        Parameters
        ----------
        campaign_name : str
            Name of the campaign.
        surveys: dict of :py:obj:`.Survey` objects, optional
            Arbitrary number of survey data objects. Default=None which implies that the campaign will be initialized
            without surveys.
        stations: :py:obj:`.Station` object, optional
            Station data (datum- and non-datum-stations). Default=None implies that the campaign will be
            initialized without station data.
        lsm_runs : list of objects inherited from :py:obj:`gravtools.models.lsm.LSM`
            Each item in the list contains one enclosed LSM object. Each LSM object reflects one dedicated run of an
            least-squares adjustment in order to estimate target parameters.

        Raises
        ------
        TypeError
            Wrong input argument type.
        """

        # Check campaign_name:
        if not isinstance(campaign_name, str):
            raise TypeError('The argument "campaign_name" needs to be a string.')
        else:
            if not campaign_name:
                raise ValueError('"campaign_name" should not be empty!')
        self.campaign_name = campaign_name

        # Check output directory:
        if not isinstance(output_directory, str):
            raise TypeError('The argument "output_directory" needs to be a string.')
        else:
            if not output_directory:
                raise ValueError('"output_directory" should not be empty!')
        self.output_directory = output_directory

        # Check surveys:
        if surveys is None:
            surveys = {}
        else:
            if not isinstance(surveys, dict):
                raise TypeError('The argument "survey" needs to be a dict of Survey objects.')
            else:
                for survey_name, survey_obj in surveys.items():
                    if not isinstance(survey_name, str):
                        raise TypeError('The argument "survey" needs to be a string.')
        self.surveys = surveys  # dict: key=Name of Survey, value=Survey object

        # Check stations:
        if stations is None:
            stations = Station()
        else:
            if not isinstance(stations, Station):
                raise TypeError('The argument "stations" needs to be a Station object.')
        self.stations = stations

        # Check lsm_runs:
        if lsm_runs is None:
            lsm_runs = []  # Empty list
        else:
            if not isinstance(lsm_runs, list):
                raise TypeError('The argument "lsm_runs" needs to be a list of LSM-objects.')
            else:
                for items in lsm_runs:
                    if not isinstance(lsm_runs, LSM):
                        raise TypeError('The argument "lsm_runs" needs to be a list of LSM-objects.')
        self.lsm_runs = lsm_runs

        # Check ref_delta_t_dt:
        if ref_delta_t_dt is not None:
            if not isinstance(ref_delta_t_dt, dt.datetime):
                raise TypeError('`ref_delta_t_dt` needs to be a datetime object.')
        self.ref_delta_t_dt = ref_delta_t_dt

        # Version of gravtools:
        self.gravtools_version = GRAVTOOLS_VERSION

    def add_survey(self, survey_add: Survey, verbose=False) -> bool:
        """Add a survey to campaign and specify whether to use it for ths analysis.

        Notes
        -----
        A survey can only be added, f the survey's name is unique within the campaign.

        Parameters
        ----------
        survey_add : :py:obj:`.Survey`
            Contains all information of s specific survey independent of the data source.
        verbose : bool, optional (default=False)
            If True, status messages are printed to the command line.

        Returns
        -------
        bool
            True, if the survey was successfully added; False, if not.
        """
        # Check if a survey with the dame name ("survey_add.name") already exists in this campaign:
        # - Raise warning:
        if survey_add.name in self.surveys.keys():
            if verbose:
                print(f'Warnung: the current campaign already contains a survey named {survey_add.name}.')
                print(' - Survey names need to be unique within a campaign.')
            return False
        else:
            # Add survey:
            self.surveys[survey_add.name] = survey_add
            if verbose:
                print(f"Survey {survey_add.name} added to the campaign.")
            return True

    def remove_survey(self, survey_name: str, verbose=False) -> bool:
        """Remove survey with the specified name from the campaign

        Parameters
        ----------
        survey_name : str
            Name of the survey that will be removed from the campaign.
        verbose : bool, optional (default=False)
            If True, status messages are printed to the command line.

        Returns
        -------
        bool
            True, if the survey was successfully removed; False, if not.
        """
        try:
            del self.surveys[survey_name]
            if verbose:
                print(f'Survey "{survey_name}" removed from campaign.')
        except KeyError:
            if verbose:
                print(f'Survey "{survey_name}" does not exist.')
            return False
        except Exception:
            if verbose:
                print(f'Failed to remove survey {survey_name}.')
            return False
        else:
            return True

    def activate_survey(self, survey_name: str, verbose=False) -> bool:
        """Set the survey with the specified name active.

        Parameters
        ----------
        survey_name : str
            Name of the survey that will be set active.
        verbose : bool, optional (default=False)
            If True, status messages are printed to the command line.

        Returns
        -------
        bool
            True, if the survey was successfully activated; False, if not.
        """
        try:
            if self.surveys[survey_name].keep_survey:
                if verbose:
                    print(f'Survey "{survey_name}" already active.')
                return True
            else:
                self.surveys[survey_name].keep_survey = True
                if verbose:
                    print(f'Survey "{survey_name}" activated.')
        except KeyError:
            if verbose:
                print(f'Survey "{survey_name}" does not exist.')
            return False
        except:
            if verbose:
                print(f'Failed to activate survey "{survey_name}"')
            return False
        else:
            return True

    def deactivate_survey(self, survey_name: str, verbose=False) -> bool:
        """Set the survey with the specified name inactive.

        Parameters
        ----------
        survey_name : str
            Name of the survey that will be set inactive.
        verbose : bool, optional (default=False)
                If True, status messages are printed to the command line.

        Returns
        -------
        bool
            True, if the survey was successfully deactivated; False, if not.
        """
        try:
            if not self.surveys[survey_name].keep_survey:
                if verbose:
                    print(f'Survey "{survey_name}" already inactive.')
                return True
            else:
                self.surveys[survey_name].keep_survey = False
                if verbose:
                    print(f'Survey "{survey_name}" deactivated.')
        except KeyError:
            if verbose:
                print(f'Survey "{survey_name}" does not exist.')
            return False
        except:
            if verbose:
                print(f'Failed to deactivate survey "{survey_name}"')
            return False
        else:
            return True

    def get_survey_names_and_status(self, verbose: bool = False) -> dict:
        """Return list with all survey names and information whether the survey is set active.

        Parameters
        ----------
        verbose : bool, optional (default=False)
            If True, survey names and status are printed to the command line.

        Returns
        -------
        dict
            The keys are the survey names and the values represent the respective status (active=True, inactive=False).
        """
        if verbose:
            print('Surveys and their status:')
        info_dict = {}
        lookup_dict = {True: 'active', False: 'inactive'}
        for surv_name, surv_obj in self.surveys.items():
            if verbose:
                activity_str = lookup_dict[surv_obj.keep_survey]
                print(f' - {surv_name:12s} ({activity_str:8s}): {surv_obj.get_number_of_observations()} observations')
            info_dict[surv_name] = surv_obj.keep_survey
        return info_dict

    @property
    def number_of_surveys(self) -> int:
        """int : Returns the number of surveys in this campaign."""
        return len(self.surveys)

    @property
    def number_of_stations(self) -> int:
        """int : Returns the number of stations in this campaign."""
        return self.stations.get_number_of_stations

    def reduce_observations_in_all_surveys(self, target_ref_height=None, target_tide_corr=None, verbose=False):
        """Reduce the observed gravity by applying the specified corrections.

        Notes
        -----
        - For this reduction vertical gravity gradients are required. They are obtained from the `Station` object.
          Hence, a `Station` object has to be attached to the Campaign object beforehand.

        - All corrections are applied on the survey-level. See py:obj:`.Survey.reduce_observations` for more details.

        Parameters
        ----------
        target_ref_height : string, specifying the target reference height type (default = `None`).
            The target reference height type has to be listed in :py:obj:`gravtools.settings.REFERENCE_HEIGHT_TYPE`.
            Default is `None` indicating that the reference heights of the input data are not changed.
        target_tide_corr : str, specifying the tidal correction type to be applied (default = `None`).
            The target tidal correction type specifies what kind of tidal correction will be applied. Valid types have
            to be listed in :py:obj:`gravtools.settings.TIDE_CORRECTION_TYPES`. Default is `None` indicating that the
            tidal corrections are not considered here (tidal corrections are inherited from input data).
        verbose : bool, optional (default=False)
            If True, status messages are printed to the command line.

        Returns
        -------
        flag_corrections_applied_correctly : bool, default = True
            `False` indicates that an error occurred when applying the observation corrections.
        error_msg : str, default = ''
            Message that describes the error in case an error occurred.

        """
        flag_corrections_applied_correctly = True
        error_msg = ''
        if verbose:
            print(f'## Reduce all observation in this campaign:')
        for survey_name, survey in self.surveys.items():
            if verbose:
                print(f'Survey {survey_name}:')
            if target_ref_height is not None:
                if verbose:
                    print(f' - Get vertical gradients')
                survey.obs_df_populate_vg_from_stations(self.stations, verbose=verbose)
            flag_corrections_applied_correctly, error_msg = survey.reduce_observations(
                target_ref_height=target_ref_height,
                target_tide_corr=target_tide_corr,
                verbose=verbose)
            # Check, if corrections were applied correctly and return error message if an error occurred:
            if not flag_corrections_applied_correctly:
                error_msg = f'Survey: {survey_name}: ' + error_msg
                if verbose:
                    print(error_msg)
                return flag_corrections_applied_correctly, error_msg
                # raise AssertionError('Reduction to reference height failed!')
        return flag_corrections_applied_correctly, error_msg

    def add_stations_from_oesgn_table_file(self, oesgn_filename, verbose=False):
        """Add station from an OESGN table file.

        Parameters
        ----------
        oesgn_filename : string, specifying the path/file of the OESGN file
            Stations in the specified OESGN table file are added to the Campaign.
        verbose : bool, optional (default=False)
            If True, status messages are printed to the command line.
        """
        self.stations.add_stations_from_oesgn_table(filename=oesgn_filename, verbose=verbose)

    def __str__(self):
        return f'Campaign "{self.campaign_name}" with {self.number_of_surveys} surveys ' \
               f'and {self.stations.get_number_of_stations} stations.'

    def synchronize_stations_and_surveys(self, verbose=False):
        """Synchronize information between station and survey data in the campaign.

        The following information is synchronized:

        - The `is_observed` flags in the :py:obj:`.Campaign.stations.stat_df` are set according to the surveys in
          :py:obj:`.Campaign.surveys`. `True` indicated that the station as observed at least once.
        - Populates the vertical gradient columns (``) of the observation DataFrames (:py:obj:`.Campaign.surveys`) with
          values from a Station object (:py:obj:`.Campaign.stations`).
        - Add observed stations

        Notes
        -----
        It is recommended to run this method whenever new stations and/or new surveys are added to the campaign on order
        to synchronize the data.

        Parameters
        ----------
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.
        """
        self.stations.stat_df['is_observed'] = False  # Reset to default.
        self.stations.stat_df['in_survey'] = None  # Reset to default.
        self.sync_observed_stations(verbose=verbose)  # Add stations from surveys.

        # Loop over all surveys to match and synchronize the survey data with the station data:
        for survey_name, survey in self.surveys.items():
            if verbose:
                print(f' - Survey: {survey_name}')
            survey.obs_df_populate_vg_from_stations(self.stations, verbose=verbose)
            self.stations.set_observed_info_from_survey(survey, verbose=verbose)

    def sync_observed_stations(self, verbose=False):
        """Adds all stations that were observed in at least one survey to this campaign's station dataframe.

        Parameters
        ----------
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.
        """
        # Loop over surveys in this campaign:
        for survey_name, survey in self.surveys.items():
            if verbose:
                print(f' - Survey: {survey_name}')
            self.stations.add_stations_from_survey(survey, verbose)

    def calculate_setup_data(self,
                             obs_type='reduced',
                             active_obs_only_for_ref_epoch=True,
                             verbose=False):
        """Calculate accumulated pseudo observations for each active setup in all active surveys.

        Notes
        -----
        Two relative reference epochs are calculated for each setup: (a) w.r.t. the first (active) observation in the
        whole campaign and (b) w.r.t. the first (active) observation in each survey. Both reference time do not differ
        for the first survey on a campaign. Whether active observations only are considered is defined by the input
        parameter `active_obs_only_for_ref_epoch`.

        Parameters
        ----------
        obs_type : str, 'observed' or 'reduced' (default)
            Defines whether the observed (as loaded from an observation file) or the reduced observations from
            `self.obs_df` are used to determine the weighted mean values per setup.
        set_epoch_of_first_obs_as_reference: bool, optional (default=True)
            If `True`, the epoch of the first observation in all surveys in the campaign is used as reference epoch
            for all surveys. `False` implies that the reference epoch is based on the first (active, if
            `active_obs_only_for_ref_epoch` is `True`) observation in each survey and determined individually for
            each survey. Stored in `self.ref_delta_t_dt` (datetime object).
        active_obs_only_for_ref_epoch: bool, optional (default=True)
            `True` implies that the relative reference epochs are determined by considering active observations only.
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.
        """
        # Get reference epoch:
        self.ref_delta_t_dt = self.get_epoch_of_first_observation(active_obs_only_for_ref_epoch)
        # if set_epoch_of_first_obs_as_reference:
        #     # Get epoch of first observation:
        #     self.ref_delta_t_dt = self.get_epoch_of_first_observation(active_obs_only_for_ref_epoch)
        # else:
        #     self.ref_delta_t_dt = None

        # Loop over all surveys in the campaign:
        if verbose:
            print(f'Calculate setup data:')
        for survey_name, survey in self.surveys.items():
            if verbose:
                print(f' - Survey: {survey_name}')
            if survey.keep_survey:
                survey.calculate_setup_data(obs_type=obs_type,
                                            ref_delta_t_campaign_dt=self.ref_delta_t_dt,
                                            active_obs_only_for_ref_epoch=active_obs_only_for_ref_epoch,
                                            verbose=verbose)

    def get_epoch_of_first_observation(self, active_obs_only_for_ref_epoch=True):
        """Returns the epoch of the first (active) observation in this campaign.

        Parameters
        ----------
        active_obs_only_for_ref_epoch: bool, optional (default=True)
            `True` implies that the reference epoch is determined by considering active observations only.

        Returns
        -------
        datetime object
        """
        first_obs_epoch_dt = None
        flag_first_survey_in_campaign = True

        for survey_name, survey in self.surveys.items():

            # Set filter to select active observations only:
            if active_obs_only_for_ref_epoch:
                filter_tmp = survey.obs_df['keep_obs'] == True  # Select active observations only
            else:
                filter_tmp = [True] * len(survey.obs_df)  # Select all observations

            if len(survey.obs_df.loc[filter_tmp, 'obs_epoch']) > 0:
                if flag_first_survey_in_campaign:
                    flag_first_survey_in_campaign = False
                    first_obs_epoch_dt = survey.obs_df.loc[filter_tmp, 'obs_epoch'].min()
                else:
                    if survey.obs_df.loc[filter_tmp, 'obs_epoch'].min() < first_obs_epoch_dt:
                        first_obs_epoch_dt = survey.obs_df.loc[filter_tmp, 'obs_epoch'].min()
        return first_obs_epoch_dt

    def initialize_and_add_lsm_run(self, lsm_method, comment='', write_log=True):
        """Initialize and add an least-squares adjustment run (object) to the campaign.

        Parameters
        ----------
        lsm_method : str
            Defines the adjustment method. Has to b listed in :py:obj:`gravtools.settings.ADJUSTMENT_METHODS`.
        comment : str, optional (default = '')
            Optional comment on the adjustment run.
        write_log : bool, optional (default=True)
            Flag that indicates whether log string should be written or not.
        """
        # Initialize LSM object:
        if lsm_method == 'LSM_diff':
            lsm_run = LSMDiff.from_campaign(self, comment, write_log)
        elif lsm_method == 'LSM_non_diff':
            lsm_run = LSMNonDiff.from_campaign(self, comment, write_log)
        elif lsm_method == 'MLR_BEV':
            lsm_run = BEVLegacyProcessing.from_campaign(self, comment, write_log)
        else:
            raise AssertionError('Unknown LSM method: {lsm_method}')
        # Add LSM object to campaign:
        self.lsm_runs.append(lsm_run)

    @property
    def lsm_run_times(self):
        """Returns a list of lsm-run times/dates that can be used to identify individual runs.

        Returns
        -------
        list : List of string stating the epochs of lsm adjustment runs that can be used to identifiterrowsy individual runs.
        """
        lsm_run_times = []
        for lsm_run in self.lsm_runs:
            lsm_run_times.append(lsm_run.time_str)
        return lsm_run_times

    def delete_lsm_run(self, idx):
        """Delete the LSM run with the specified index in the list.

        Parameters
        ----------
        idx : int
            Index of the LSM object in the list.
        """
        if idx != -1:
            del self.lsm_runs[idx]

    def set_reference_time(self, ref_delta_t_dt):
        """Set refernce time for the determination of relative time spans, e.g. for the drift polynomial.

        Parameters
        ----------
        ref_delta_t_dt : datetime object
            Reference time epoch w.r.t. UTC.
        """
        if isinstance(ref_delta_t_dt, dt.datetime):
            self.ref_delta_t_dt = ref_delta_t_dt
        else:
            raise ValueError('`ref_delta_t_dt` needs to be a datetime object.')

    def write_nsb_file(self, filename: str, lsm_run_index, vertical_offset_mode: str = 'first', verbose=True):
        """Write the results of an LSM run to an nsb file (input for NSDB database).

        Parameters
        ----------
        filename : str
            Name and path of the output nsb file (e.g. /home/johnny/example.nsb)
        lsm_run_index : int
            Index of the lsm run in `campaign.lsm_runs` of which the results are exported to the nsb file.
        vertical_offset_mode : str, optional (default='first')
            Defines how the vertical offsets between instrument top and ground (dhb) and reference marker (dhf),
            respectively, are determined in the case of multiple measurements (setups) on the same point. In the nsb
            file only one dhf/dhb pair per station is allowed. Two options: (1) 'first' indicates that dhb and dhf are
            taken from the first setup at a station. (2) 'mean' indicates that mean values over all setups are taken.
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.
        """
        # Init.:
        nsb_string = ''

        # Get and prepare data:
        # - lsm_run
        lsm_run = self.lsm_runs[lsm_run_index]
        results_stat_df = lsm_run.get_results_stat_df

        # Loop over stations in results dataframe:
        for index, row in results_stat_df.iterrows():
            station_name = row['station_name']
            observed_in_surveys = []
            dhb_list_m = []
            dhf_list_m = []

            # Get surveys at which the station was observed:
            for survey_name, setup_df in lsm_run.setups.items():
                if len(setup_df.loc[setup_df['station_name'] == station_name]) > 0:  # was observed in this setup!
                    observed_in_surveys.append(survey_name)
                    obs_df = self.surveys[survey_name].obs_df
                    setup_ids = obs_df.loc[obs_df['station_name'] == station_name, 'setup_id'].unique()
                    setup_ids = setup_df.loc[setup_df['station_name'] == station_name, 'setup_id'].to_list()
                    # Get list of dhb and dhf:
                    for setup_id in setup_ids:
                        dhb_list_m.append(obs_df.loc[obs_df['setup_id'] == setup_id, 'dhb_m'].values[0])
                        dhf_list_m.append(obs_df.loc[obs_df['setup_id'] == setup_id, 'dhf_m'].values[0])

            if vertical_offset_mode == 'first':
                dhb_m = dhb_list_m[0]
                dhf_m = dhf_list_m[0]
            elif vertical_offset_mode == 'mean':
                dhb_m = np.mean(dhb_list_m)
                dhf_m = np.mean(dhf_list_m)

            # Get gravimeter S/N and gravimeter type of first survey in the list:
            if verbose:
                if len(observed_in_surveys) > 1:
                    print(f'WARNING: station {station_name} was observed in {len(observed_in_surveys)} surveys! Hence, '
                          f'the gravimeter serial number and type may be ambiguous in the nsb file!')
            gravimeter_type = self.surveys[observed_in_surveys[0]].gravimeter_type
            gravimeter_serial_number = self.surveys[observed_in_surveys[0]].gravimeter_serial_number
            date_str = self.surveys[observed_in_surveys[0]].date.strftime('%Y%m%d')

            nsb_string += '{:10s} {:8s}  {:9.0f} {:3.0f} {:1s} {:>4s} {:4.0f} {:4.0f}\n'.format(
                station_name,
                date_str,
                row['g_est_mugal'] + ADDITIVE_CONST_ABS_GRTAVITY,
                row['sd_g_est_mugal'],
                GRAVIMETER_TYPES_KZG_LOOKUPTABLE[gravimeter_type],
                GRAVIMETER_SERIAL_NUMBER_TO_ID_LOOKUPTABLE[gravimeter_serial_number],
                (dhb_m + GRAVIMETER_REFERENCE_HEIGHT_CORRECTIONS_m[gravimeter_type]) * 100,
                (dhf_m + GRAVIMETER_REFERENCE_HEIGHT_CORRECTIONS_m[gravimeter_type]) * 100,
            )

        # Write file:
        with open(filename, 'w') as out_file:
            out_file.write(nsb_string)

    def save_to_pickle(self, filename=None, verbose=True):
        """Save the campaign object to a pickle file at the given path.

        Parameters
        ----------
        filename : str, optional (default=`None`)
            Path and name of the pickle file, e.g. /home/user1/data/camp1.pkl. `None` indicates that the campaign object
            is saved to the default output directory past (`campaign.output_directory`). In this case the file is named
            `<campaign_name>.pkl`.
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.

        Returns
        -------
        str : Name and path of the saved file.
        """
        if filename is None:
            filename = os.path.join(self.output_directory, f'{self.campaign_name}.pkl')
        # Open file:
        if verbose:
            print(f'Export campaign data to {filename}.')
        with open(filename, 'wb') as outfile:
            pickle.dump(self, outfile, protocol=pickle.HIGHEST_PROTOCOL)
        return filename

    @classmethod
    def from_pkl(cls, filename: str, verbose: bool = True):
        """Loads a campaign object from a pickle file.

        Parameters
        ----------
        filename : str
            Name and path of the pickle file containing the campaign data.
        verbose : bool, optional (default=False)
            If `True`, status messages are printed to the command line.

        Returns
        -------
        `Campaign` object
        """
        with open(filename, 'rb') as handle:
            campaign = pickle.load(handle)
        if verbose:
            print(
                f'Loaded campaign "{campaign.campaign_name}" with {campaign.number_of_stations} station(s) and {campaign.number_of_surveys} survey(s).')
        return campaign

    def set_output_directory(self, output_directory):
        """Change the campaign's output directory.

        Parameters
        ----------
        output_directory : str
            New output directory. The specified directory has to exist on the computer/os!
        """
        # Check output directory:
        if not isinstance(output_directory, str):
            raise TypeError('The argument "output_directory" needs to be a string.')
        else:
            if not os.path.isdir(output_directory):
                raise AssertionError(f'The directory "{output_directory}" does not exist!')
        self.output_directory = output_directory


if __name__ == '__main__':
    """Main function, primarily for debugging and testing."""
    filename = '/home/heller/pyProjects/gravtools/out/test.pkl'
    camp = Campaign.from_pkl(filename)
