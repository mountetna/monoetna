import contextlib
import tempfile
from typing import Optional, ContextManager

import os

import stat

from airflow.exceptions import AirflowException
from airflow.models import Connection


@contextlib.contextmanager
def prepared_key_from(connection: Connection, assert_password_non_empty=False) -> ContextManager[Optional[str]]:
    """
    Attempts to create an ssh compatible key file from the password field of a connection.
    Use this with a context eg:
        with prepared_key_from(conn) as key_file_path:
          # only use the key file path here.

    NOTE: The file that backs the key will be deleted once the context is left, so be mindful.
    The file will only be readable by this process user.
    """
    # The 'password field' does not support multiline currently, so when a key is pasted in it,
    # it does weird things to its formatting.  This function handles unpacking the connection.password
    # into a workable key file to use.
    if connection.password:
        # prepare a key.
        with tempfile.NamedTemporaryFile() as file:
            password = connection.password
            words = password.split(' ')
            lines = []
            header = False

            for word in words:
                if header:
                    lines[-1] += " " + word
                else:
                    lines.append(word)

                if word.startswith('--'):
                    header = True
                if word.endswith('--'):
                    header = False

            file.write('\n'.join(lines).encode('utf8'))
            os.chmod(file.name, stat.S_IRUSR | stat.S_IWUSR)
            file.flush()
            yield file.name
    else:
        if assert_password_non_empty:
            raise AirflowException(f"RSA key was not set for connection '{connection.conn_id}'")
        yield None
