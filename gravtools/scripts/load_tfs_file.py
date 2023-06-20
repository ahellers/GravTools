"""Example script for loading TSoft TSF files with tidal data.

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

from gravtools.tides.tide_data_tfs import TSF

def main():
    """Main function"""

    # Load file with integer seconds time reference:
    FILENAME_TFS = '../../data/tfs_files/TEST_SITE_1.TSF'
    tfs_data = TSF.from_tfs_file(FILENAME_TFS)
    print(tfs_data.filename)
    print(tfs_data)
    print(tfs_data.get_channel_df(1))

    # Load file with integer sub-seconds time reference:
    FILENAME_TFS = '../../data/tfs_files/TEST_SITE_1_subsecond.TSF'
    tfs_data = TSF.from_tfs_file(FILENAME_TFS)
    print(tfs_data.filename)
    print(tfs_data)
    print(tfs_data.get_channel_df('last'))
    print(tfs_data.get_channel_df('TEST_SITE_1:Theory Sol. Earth:WDD'))

    print('end')

if __name__ == '__main__':
    main()