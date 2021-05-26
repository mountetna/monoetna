
activate_control_app
plugin :yabeda
plugin :yabeda_prometheus

bind "tcp://0.0.0.0:3000"
threads 3, 16
