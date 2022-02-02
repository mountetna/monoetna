from setuptools import setup

setup(
    name='etna',
    version="0.0.1",
    description='Mount Etna airflow library code',
    long_description="",
    long_description_content_type='text/markdown',
    entry_points={
        "apache_airflow_provider": [
            "provider_info=etna.get_provider_info:get_provider_info"
        ]
    },
    license='Apache License 2.0',
    packages=[],
    # packages=[
    #     'etna', 'etna.xcom', 'etna.auth', 'etna.executors', 'etna.hooks',
    #     'etna.operators', 'etna.utils', 'etna.operators', 'etna.dags',
    #     'etna.etls'
    # ],
    # install_requires=['apache-airflow>=2.2'],
    # setup_requires=['setuptools', 'wheel'],
    python_requires='>=3.7'
)