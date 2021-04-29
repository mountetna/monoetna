# Prometheus

Scrapping (pull) based monitoring and alerting service from google.

## So... not push?

Right.  We don't _push_ metrics in real time into prometheus.  We collect and buffer metrics in our services and then have
prometheus collect them over time.

This has some downsides to keep in mind

1.  No guarantee of delivery.  Most metrics are just buffered in memory and thus process restart or bug can lose them.
2.  Not real time.  Prometheus is configured with a collection target which aims for efficiency of the system, not real time collection.
3.  Routing.  Prometheus needs access to each target, which means configuration involves things like service discovery and networking configuration.
4.  Sample-oriented.  Prometheus is focused on gauges, windows, histograms, and sample based approaches to data.  This means by default most data is not perfectly accurate to infinite granularity.

So, why?  Why not use something like statsd?

1.  Higher guarantee of delivery, including retry and file system buffering.
2.  Real time -- send metrics and show it immediately!
3.  Only concerned about routing to the central system, service discovery only necessary for one system.
4.  No need to buffer or sample, collect every metric always!

This sounds heavily in favor of push systems, but there is one.  major.  flaw.  with that approach.

Reliability.

A systems monitoring tool has ONE job.  To be up, and alerting, when the rest of the system is having issues.
There is no point to having real time, highly granular data that is experiencing network congestion, latency, 
and storage issues in unpredictable ways.  Real time push systems must always constantly scale above and beyond the
system it is responsible for or risk not having the most important data available at the moment the system is failing.

PULLing buffered, sampled data on regular intervals is systematically more reliable.  There is no unpredictable growth.
If data grows, sampling simply squashes that data into consistent bucketed sizes.  If the system's network is unreliable,
pulling ensures that only a small amount of reliability is required to deliver consistently small payloads.  Pulling
allows the system to prioritize some metrics over others.  Pulling allows the controller to estimate data lag by distance
from last known pull, no matter the network conditions.

## Ok, how do I setup prometheus locally?

It will run as part of `make up` just fine.  But you'll likely need to make some network changes to enable access.

First, ensure you have a hosts entry for `prometheus.development.local`.  Secondly, you'll need to put your local docker
into swarm mode.  Thirdly, you'll need to enable docker metrics in the daemon.json file.

To do all of this correctly, you need to decide on which network you want to make this available to prometheus.  Remember,
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
  
systemctl restart docker
```

Ok cool, remember to `make up` to bring up your services once docker has finished starting.

Last step, you'll need to install the node exporter as well.  Fortunately I have a script to do this for you.
Just run ./bin/node_exporter in a separate tab and keep this process open in the background.  You can kill it anytime
you like, it's not crucial.
