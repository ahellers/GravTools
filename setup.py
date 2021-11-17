import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="gravtools",
    version="0.0.2",
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
    python_requires='>=3.6',

    install_requires=[
        'cycler==0.11.0',
        'fonttools==4.28.1',
        'joblib==1.1.0',
        'kiwisolver==1.3.2',
        'matplotlib==3.5.0',
        'numpy==1.21.4',
        'packaging==21.2',
        'pandas==1.3.4',
        'Pillow==8.4.0',
        'pyparsing==3.0.6',
        'PyQt5==5.15.6',
        'PyQt5-Qt5==5.15.2',
        'PyQt5-sip==12.9.0',
        'pyqtgraph==0.12.3',
        'python-dateutil==2.8.2',
        'pytz==2021.3',
        'scikit-learn==1.0.1',
        'scipy==1.7.2',
        'setuptools-scm==6.3.2',
        'six==1.16.0',
        'threadpoolctl==3.0.0',
        'tomli==1.2.2',
    ],

    entry_points={
        "console_scripts": [
            # "schwaus=bev_legacy.command_line:schwaus",
            "gt=gravtools.scripts.run_gui:run_gui",
        ]
    }
)