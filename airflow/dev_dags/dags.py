from airflow.decorators import dag

from etna import run_in_swarm


@dag()
def weird_dag():
    run_in_swarm('test_swarm_operator', 'test-service-swarm', ['bach', '-c', "echo 1"], output_json=True)
