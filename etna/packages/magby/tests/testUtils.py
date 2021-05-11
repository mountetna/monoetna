from betamax.cassette import cassette
from betamax import Betamax
from betamax_serializers import pretty_json
import os

from ..magby.Magma import *


def prepCassette(session: Session, cassette_dir: str):
    vcr = Betamax(session)
    vcr.register_serializer(pretty_json.PrettyJSONSerializer)
    if not os.path.exists(cassette_dir):
        os.mkdir(cassette_dir)
    with vcr.configure() as config:
        config.cassette_library_dir = cassette_dir
        config.default_cassette_options['serialize_with'] = 'prettyjson'
        config.before_record(callback=sanitizeToken)

    return vcr


def sanitizeToken(interaction: Dict, current_cassette: cassette) -> None:
    headers = interaction.data['request']['headers']
    token = headers.get('Authorization')[0]
    if token is None:
        return
    current_cassette.placeholders.append(
        cassette.Placeholder(placeholder='<AUTH_TOKEN>', replace=token)
    )