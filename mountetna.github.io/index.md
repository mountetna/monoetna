---
layout: default
---

## Introduction to Mount Etna

Mount Etna is a set of applications that together provide tools for storing,
organizing, viewing and analyzing research data sets.

Data sets in Mount Etna are broken up into "projects", with each application in
the system providing access to different pieces of the dataset.

The principal actors are:

* [Janus]({{ site.baseurl }}{% link janus.md %}) - the **authentication service**, controls role-based access to each project

* [Magma]({{ site.baseurl }}{% link magma.md %}) - the **data warehouse**, organizes each project's data into models and records.

* [Metis]({{ site.baseurl }}{% link metis.md %}) - the **file service**, holds
  project-related files. Records in Magma may link to binary files stored on
  Metis.

* [Timur]({{ site.baseurl }}{% link timur.md %}) - the **data browser**, allows viewing, searching and analyzing project data from Magma.

* [Vulcan]({{ site.baseurl }}{% link vulcan.md %}) - performs **data analysis** on project data from Magma and Metis.
