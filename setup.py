from setuptools import setup

setup(
    name='LocusSolver',
    version='0.1.0',
    description='Merge GFF3 annotations',
    url='',
    author='Christopher Neely',
    author_email='christopher.neely1200@gmail.com',
    license='GNU GPL 3',
    packages=['src', "src.filters", "src.util"],
    python_requires='>=3.8',
    install_requires=[
        "bcbio-gff==0.6.6",
        "biopython==1.76",
        "plumbum==1.6.6"
    ],
    scripts=[
        "src/locus_solver",
    ]
)
