world=World!

applyTemplate testerb.erb "$dest"/abc.conf
enableConfig
# docker kill --signal="USR1" httpd
cat /config/cur/abc.conf