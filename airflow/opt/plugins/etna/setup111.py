from setuptools import find_packages, setup

setup(
    name='etna',
    version="0.0.1",
    description='Mount Etna airflow library code',
    long_description="",
    long_description_content_type='text/markdown',
    entry_points={
        "apache_airflow_provider": [
            "provider_info=etna.__init__:get_provider_info"
        ]
    },
    license='Apache License 2.0',
    packages=['etna', 'etna.auth', 'etna.executors', 'etna.hooks', 'etna.operators', 'etna.utils'],
    install_requires=['apache-airflow>=2.2'],
    setup_requires=['setuptools', 'wheel'],
    python_requires='~=3.7',
)