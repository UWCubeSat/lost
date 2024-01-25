from setuptools import setup

# TODO: add error checking and/or automated build if /lost is missing
# TODO: handle multiple architectures

setup(
    name='lost',
    version='0.0.0',
    packages=['lost'],
    package_dir={'lost': 'lostpy'},
    include_package_data=True,
    package_data={'lost': [
            'lost',  # the LOST executable
            'bright-star-catalog.tsv',  # the stars are needed...
            'tmp/tmp_database.dat',  # dummy database so pip uninstall deletes generated database
    ]},
    install_requires=[
        'imageio',
    ],
    zip_safe=False
)