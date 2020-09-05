default['hosts'] = {
    'magma' => 'https://magma.ucsf.edu',
    'metis' => 'https://metis.ucsf.edu',
    'timur' => 'https://timus.ucsf.edu',
    'janus' => 'https://janus.ucsf.edu',
    'polyphemus' => 'https://polyphemus.ucsf.edu',
}

default['waiver'] = {
    'release_bucket' => '',
    'restrict_bucket' => '',
}

default['janus_token'] = {
    'name': 'JANUS_TOKEN',
    'algo': 'RS256',
}

default['users'] = [
    {
        name: 'root',
    }
]

default['docker'] = {
    'data_root' => '/var/run/docker',
    'host' => 'unix:///var/run/docker.sock',
    'default_tag' => 'master',
    'swarm' => {
        'token' => nil,
        'leader_addr' => nil,
    },
    'registry' => {
        'url' => 'https://hub.docker.io',
        'username' => nil,
        'password' => nil,
        'email' => nil,
    },
}

default['slack_notifications'] = {
    'watchtower' => {
        'hook_url' => '',
        'channel' => '#bioinformatics-ping'
    }
}

default['polyphemus']['janus_token'] = nil

default['rollbar'] = {}
