# Prometheus

Scrapping (pull) based _sample_ monitoring and alerting service from google.

## But I need realtime streams of data.

We can decrease the scrapping interval if it is needed.  That said, Prometheus is not going to provide real time data.  
It is focused on easy maintenance and reliability first and foremost.  Many of its decisions are counterintuitive at
first glance, but actually guarantee the least systems impact even as you scale the number of metrics we collect
over time.  It is efficient, reliable, and low maintenance.

# But I need non sampled data!

Prometheus *IS LOSEFUL*.  It doesn't collect discrete data points, but rather aggregatable samples.  
Most of our collectors exist in memory and will drop data on process restart.
Most collectors uses histograms to approximate percentile distributions.  
That is because prometheus's primary object is maintainability and low systems impact.  
That said, if you need high fidelity, loseless data, all you need to do is buffer that
data somewhere and then provide it as a scrapping target to prometheus.  Or, not use prometheus.

Overtime, I can provide examples of alternative strategies, but I highly recommend considering wether a sampled
approach is sufficient, due to its low impact and easy of implementation.

## Ok, how do I setup prometheus locally?

It will run as part of `make up` just fine.  But you'll likely need to make some network changes to enable access.

First, ensure you have a hosts entry for `prometheus.development.local` and `grafana.development.local` if you intend  
to use those UIs.  Secondly, you'll need to put your local docker
into swarm mode.  Thirdly, you'll need to enable docker metrics in the daemon.json file.

To do all of this correctly, you need to decide on which network you want to make available to your swarm.  Remember,
prometheus is connecting inside docker and does not have access to the host's loopback interface.  This network is likely
the one your ethernet or wireless card is connected to, since docker does expose those networks to inner containers.


```
# This is for wireless devices on linux, you may also try eth for instance
[home@excalibur:~/monoetna]$ ifconfig -s | grep wlp
wlp0s20f  1500  1082411      0      0 0        188175      0      0      0 BMRU
[home@excalibur:~/monoetna]$ ifconfig wlp0s20f3 | grep inet
        inet 10.0.0.115  netmask 255.255.255.0  broadcast 0.0.0.0
        inet6 fe80::362e:b7ff:fe9d:e4b9  prefixlen 64  scopeid 0x20<link>
        inet6 2601:646:c401:3bd0::e07e  prefixlen 128  scopeid 0x0<global>
        inet6 2601:646:c401:3bd0:102b:1ad9:93a7:f4c5  prefixlen 64  scopeid 0x0<global>
        inet6 2601:646:c401:3bd0:362e:b7ff:fe9d:e4b9  prefixlen 64  scopeid 0x0<global>
[home@excalibur:~/monoetna]$ docker swarm init --advertise-addr 10.0.0.115
```

Don't worry about the token it displays, in this case only our local system will be part of our swarm.  Make sure
the `advertise-addr` is the inet address of the device you are going to be relying on for connection from inside docker.

Ok, cool.  Now more than likely, you'll also need to open up some ports since accessing this address from the device
ip will hit your firewall (it's not the loopback, again).  You'll want to open TCP ports 7100, 9323, and 9100.

Now, for your docker daemon.json file, you'll want to merge in these values:

```
sudo vim /etc/docker/daemon.json
  1 {
  2   "metrics-addr" : "0.0.0.0:9323",
  3   "experimental" : true
  4 }
  
systemctl restart docker # Restart docker!
```

Ok cool, remember to `make up` to bring up your services once docker has finished starting.

Last step, you'll need to install the node exporter as well.  Fortunately I have a script to do this for you.
Just run ./bin/node_exporter in a separate tab and keep this process open in the background.  You can kill it anytime
you like, it's not crucial.
