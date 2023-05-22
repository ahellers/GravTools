"""Example script for creating a corretcion time series by loading data from TSF files.

Copyright (C) 2023  Andreas Hellerschmied <andreas.hellerschmied@bev.gv.at>

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
from datetime import datetime

import numpy as np

from gravtools.tides.tide_data_tfs import TSF
from gravtools.tides.correction_time_series import CorrectionTimeSeries

def main():
    """Main function"""

    # Load file with integer seconds time reference:
    FILENAME_TFS = '../../data/tfs_files/l230406.TSF'

    corr_time_series = CorrectionTimeSeries()
    corr_time_series.load_tfs_file('l230406', FILENAME_TFS, instrument='Theory-Loading')

    corr_time_series.surveys['l230406'].stations['BEV_U3'].tidal_correction.to_df()

    corr_time_series.load_tfs_file('l230406', FILENAME_TFS, instrument='Theory-Loading')


    interp_times = np.array([datetime.strptime('2023-04-06 00:00:30', '%Y-%m-%d %H:%M:%S'),
                    datetime.strptime('2023-04-06 00:01:30', '%Y-%m-%d %H:%M:%S'),
                    datetime.strptime('2023-04-06 00:02:30', '%Y-%m-%d %H:%M:%S')])  # or list()

    interpol_y = corr_time_series.surveys['l230406'].stations['BEV_U3'].tidal_correction.interpolate(interp_times=interp_times)
    print(interpol_y)

    print('end')





if __name__ == '__main__':
    main()