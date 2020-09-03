default['users'] = [
    {
        name: 'root',
    }
]

default['docker'] = {
    'data_root' => '/var/run/docker',
    'host' => 'unix:///var/run/docker.sock',
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