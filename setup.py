import setuptools
import sys
import os


# Add current work directory to python path to get the version tag:
sys.path.append(os.getcwd())
import gravtools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="gravtools",
    version=gravtools.__version__,
    author="Andreas Hellerschmied",
    author_email="andreas.hellerschmied@bev.gv.at",
    description="Gravgui"
                "ity surveys utilities.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    # url="https://gitlab.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6.1',

    install_requires=[
        'cycler==0.11.0',
        'joblib==1.1.0',
        'kiwisolver==1.3.1',
        'matplotlib==3.3.4',
        'numpy==1.19.5',
        'pandas==1.1.5',
        'Pillow==8.4.0',
        'pyparsing==3.0.6',
        'PyQt5==5.15.6',
        'PyQt5-Qt5==5.15.2',
        'PyQt5-sip==12.9.0',
        'pyqtgraph==0.11.1',
        'python-dateutil==2.8.2',
        'pytz==2021.3',
        'scikit-learn==0.24.2',
        'scipy==1.5.4',
        'six==1.16.0',
        'threadpoolctl==3.0.0',
    ],

    entry_points={
        "console_scripts": [
            "gt=gravtools.scripts.run_gui:run_gui",
        ]
    }
)