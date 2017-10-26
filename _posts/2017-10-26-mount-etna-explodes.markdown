---
layout: post
title:  "Mount Etna Explodes!"
date:   2017-10-26 23:13:05 +0000
categories: etna 
---
But don't worry, we're comfortable working amidst volcanic ejecta!

This is our inaugural post. In this space we intend to document the travails
and successes of our effort to create a fantastic bioinformatics platform and
facilitate excellent research at UCSF (and elsewhere).

For the uninitiated, Mount Etna is a series of interacting applications that
present a *data system* which organizes all of the different kinds of data
associated with scientific research projects.

The key players:

- *Magma* - the data warehouse, magma disgorges red-hot data to all interested parties using an expressive query interface, and gobbles up new data using friendly loaders.
- *Timur* - the data browser, where data may be viewed, downloaded, and analyzed. Timur collects data from various Mount Etna services, interfaces with calculation services, provides custom data views tailored to specific biological applications, and provides general query and plotting tools for exploring your dataset.
- *Janus* - an authentication service that holds user and project permissions.
- *Polyphemus* - a workhorse who manages dispatch of work on Magma data.
- *Metis* - a file service which serves up large binary files, but also supports interesting partitions of those files (e.g. BAM slicing).
