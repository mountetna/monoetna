world=World!
template testerb.erb "$dest"/abc.conf

function onApply() {
  enableConfig
  # docker kill --signal="USR1" httpd
  cat /config/cur/abc.conf
}