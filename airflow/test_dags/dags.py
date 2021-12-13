from airflow.decorators import dag, task

@dag()
def weird_dag():
    @task
    def task_1():
        return 123

    @task
    def task_2(value):
        return 234 + value

    task_2(task_1())
