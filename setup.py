import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="gravtools",
    version="0.0.1",
    author="Andreas Hellerschmied",
    author_email="andreas.hellerschmied@bev.gv.at",
    description="Gravity surveys utilities.",
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
    entry_points={
        "console_scripts": [
            "schwaus=gravtools.scripts:schwaus",
        ]
    }
)