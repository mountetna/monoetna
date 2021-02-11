from betamax.cassette import cassette
from betamax import Betamax
from betamax_serializers import pretty_json
from requests import Session


def prepCassette(session: Session()):
    vcr = Betamax(session)
    vcr.register_serializer(pretty_json.PrettyJSONSerializer)
    with vcr.configure() as config:
        config.cassette_library_dir = './magby/tests/fixtures/cassettes'
        config.default_cassette_options['serialize_with'] = 'prettyjson'
        config.before_record(callback=sanitizeToken)

    return vcr


def sanitizeToken(interaction, current_cassette):
    headers = interaction.data['request']['headers']
    token = headers.get('Authorization')[0]
    if token is None:
        return
    current_cassette.placeholders.append(
        cassette.Placeholder(placeholder='<AUTH_TOKEN>', replace=token)
    )