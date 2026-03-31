# ADR-0006: Project Dashboard

## Status

Proposed

## Context

Currently the experience of using the data library applications is diffuse;
there are several applications covering various scopes of work, but the user
is left to navigate each of these applications without assistance. The new
user, especially, may feel adrift when beginning with a blank project, and may
not know how to navigate to the necessary functions in each application.

By contrast there is a clear workflow that has emerged on how, generally, data
should move in the library.

1. Create the project - usually initiated by a library administrator, adding an
   admin user to the project.
2. Add users to the project
3. Create models in Magma
4. Create a naming grammar in Gnomon
5. Create buckets for holding data in Metis
6. Ingest data to Metis via metis_client
7. Create data loaders in Polyphemus
8. Run data loaders to link records into Magma
9. Query, search, etc. data
10. Analyze data using Vulcan workflows

## Decision

We will create a project dashboard, which will appear for the moment in the
Timur data browsing app, and will be the primary landing point directed to by
Vesta landing page app.

The dashboard will summarize project state across each of the library
applications, and will also direct the user on how to complete outstanding
tasks to facilitate the above progression.

The app components will each present data and actions appropriate to the user's
role in the project. For viewers, the component will summarize the project
data. For editors, it will link to tools or documentation to add or edit data.
For admins, it will link to tools or documentation to configure the project.

HEADER

PROJECT          Project
TEXT             Dashboard

                 Map

Janus
  viewer: Show count of users in project broken down by role (m admins, n editors, o viewers), link to janus
  admin: Show "Add new users" link to Janus project page if user count is 0

Timur:
  viewer: Show count of defined models. Show count of total project records if models exist, link to search
  admin: Show "Add new models" link to map page if model count is 0, show link to modeling docs

Gnomon:
  viewer: Show count of rules / models with identifiers. Show count of created identifiers.
  editor: Show "Create grammar" link to Gnomon project page, show link to naming docs if grammar is absent. Show "Add identifiers" link to gnomon if gnomon mode allows.

Metis:
  viewer: Show count of buckets. Show total count of files. Show total byte count. -> download links
  editor: Show "Upload data" link to metis_client docs if file count is less than 10.
  admin: Show "Create data buckets" link to metis if bucket count is 0

Polyphemus:
  viewer: Show count of data loaders. Show last runtime for any loader.
  editor: Show "Create data loader" link to polyphemus if data loader count is 0.

Vulcan:
  viewer: Show count of workspaces created. Show count of project workflows.
  editor: Show "Create workspace" link to vulcan if workspace count is 0.

## Consequences

Users will have a more streamlined experience into the library, with a single
reference point they can consult to maintain and improve the state of the
project.
