from .etna_base import EtnaClientBase

class Janus(EtnaClientBase):
    def generate_token(
        self, project_name: str, read_only=True, token_type="task"
    ) -> bytes:
        response = self.session.post(
            self.prepare_url("api", "tokens", "generate"),
            json=dict(
                token_type=token_type, project_name=project_name, read_only=read_only
            ),
        )
        return response.content

    def nonce(self) -> bytes:
        response = self.session.get(self.prepare_url("api", "tokens", "nonce"))
        return response.content
