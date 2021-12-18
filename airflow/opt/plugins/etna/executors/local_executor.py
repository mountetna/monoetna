from airflow.executors.local_executor import LocalExecutor as OrigLocalExecutor
from docker import APIClient

from etna.operators.swarm_operator import swarm_cleanup

class LocalExecutor(OrigLocalExecutor):
    def start(self):
        try:
            cli = APIClient(base_url="unix://var/run/docker.sock")
            swarm_cleanup(cli)
        except Exception as e:
            self.log.error("Failed to perform swarm cleanup %s.", str(e))
        else:
            self.log.info("Swarm cleanup completed successfully")

        super(LocalExecutor, self).start()
