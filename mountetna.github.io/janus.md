---
layout: default
group: app
---

# Janus
{:.no_toc}
## Authentication Service
{:.no_toc}

Janus is an authentication and identity service for Etna applications. It is based on the [etna](https://github.com/mountetna/etna) gem.

Janus implements the basic user and project structure of Etna applications.

* TOC
{:toc}

## Users

A Janus user is primarily identified by an email address. You may add users via the `bin/janus add_user` command.

See below on how users may authenticate.

## Projects

A 'project' is the entity that produced a particular data set. You may add projects using the `bin/janus add_project` command. Generally Etna applications communicate using the `project_name` and rarely `The Full Name of the Project`.

## Permissions

Users are granted specific permissions on each project. A permission consists of:
  * user
  * project
  * role - either of `[ 'administrator', 'editor', 'viewer' ]`
  * restricted - true if the user can see the project's restricted data

You may add permissions using the `bin/janus permit` command.

# Identification

Janus provides identity by yielding a JSON web token (JWT). Client applications may verify this token using Janus's public key, most likely through the Etna::Auth rack middleware.

The token format is: `<header>.<params>.<signature>`

Each section is a base64-encoded JSON hash. The params Janus reports are: `{ email, first, last, perm }` - the latter encodes the user's project permissions.

# Authenticating

There are three ways to get a token:

## Password login

The endpoint `/login` can be configured to display an HTML form for password entry. If successful, Janus will set a cookie with the token in the response. This endpoint is mostly useful for developers.

## Shibboleth login

The `/login` endpoint can also be configured as a Shibboleth-protected endpoint for authentication. If successful, Janus will set a cookie with the token in the response. This endpoint is most suitable for browser applications. To use this method for authentication, set `auth_method: shibboleth` in the configuration.

## Public-key login

Machine users who cannot use browsers can use a registered public key to generate a token.

The endpoint `/api/tokens/nonce` returns a cryptographic nonce. The user signs the nonce, base64-encodes the signature, and concatenates the result to the nonce.

The endpoint `/api/tokens/generate` returns a valid token if the `Authorization` header is set to the appropriate value.

The token generation api can be accessed via the `etna` gem, using the command `etna administrate tokens generate --email your@email.id`.

Here is a bash script that will successfully generate a Janus token on most systems (you need 'openssl', 'wget' and 'base64' utilities).

```
#!/bin/bash

# The base URL for janus
JANUS_URL=$1
# Your secret key file in PEM format
PEM=$2
# Your email address
EMAIL=$(echo -n $3 | base64 -w 0)

NONCE=$(wget -q -O - $JANUS_URL/api/tokens/nonce)

SIG=$(echo -n $NONCE.$EMAIL | openssl dgst -sha256 -sign $PEM | base64 -w 0)

AUTH=$NONCE.$EMAIL.$SIG

TOKEN=$(wget --post-data='' -q -O - --header="Authorization: Signed-Nonce $AUTH" $JANUS_URL/api/tokens/generate )

echo $TOKEN
```

### Task tokens

Janus tokens are usually limited in lifespan as they cannot be easily
invalidated. However, this makes them unsuitable for many tasks (e.g. uploading
to Metis) since they may expire before the task can be completed. Therefore
Janus will also issue a "task" token with a longer lifespan and limited
permissions for a single project.

To acquire a task token, you may pass the argument `token_type=task` to the
endpoint `/api/tokens/generate`, along with a `project_name` argument.
Alternatively, you may use the `etna` gem command `etna administrate token
generate --email your@email.id --task --project-name your_project`.

Task tokens may be generated either using a signed nonce or using a regular
token.

# Configuration

Janus is an Etna application and puts all of its configuration into a `config.yml` YAML file.

Example:

`./config.yml`

```
---
:development:
  :db:
    :database: janus
    :host: localhost
    :adapter: postgres
    :encoding: unicode
    :username: developer
    :password: <%= @psql_developer_password %>
    :pool: 5
    :timeout: 5000

    # We recommend using a search path that is not 'public' (for postgres only).
    :search_path: [private]

  # How Janus should generate passwords (using Etna::SignService).
  :pass_algo: sha256
  :pass_salt: <%= @janus_password_salt %>

  # Token generation options.
  :token_domain: ucsf.edu
  :token_life: 86400
  :token_algo: RS256
  :token_name: <%= @janus_token_name %>

  :log_level: debug
  :log_file: '<%= @log_file %>'

  # Janus private key (used to generate the user token).
  :rsa_private: |
<%= @rsa_private %>

  # Janus public key. This is required in the config.yml files of other
  # servers/apps that use the etna gem. We keep it here for reference.
  :rsa_public: |
<%= @rsa_public %>
```

## Generating keys

Janus relies on a RSA public and private key pair. You may generate keys in PEM
format using the command `bin/janus generate_key_pair <key_size>`.

# User and Project Setup

## Via Command line

### Creating a user

You may add a new user with the `add_user` command. The primary identifier for
a user is an email address. They may also have a first and last name. You may optionally set a password here.

### Setting a public key

Some users will want to set a public key to allow them to generate a janus token via the `/generate` endpoint (see above). You may set this key using the `add_user_key` command and a public key file. Keys must be in PEM format and must be RSA keys.

To create a public / private key pair, you can use the following commands:

```
$ openssl genrsa -out privatekey.pem 2048
$ openssl rsa -in privatekey.pem -out publickey.pem -pubout -outform PEM
```

And the resulting file (`publickey.pem`) should look like:

```
-----BEGIN PUBLIC KEY-----
MIIBIjANBgkqhkiG9w0BAQEFAAOCAQ8AMIIBCgKCAQEArzKDJ1+08l3y12nMzYz3
h3xo4n/4FZ+jrnWTsPa4F6OxYo9J1E0PBjHLbSIvbIMjNcf5OTywzjzDS8Hh6P20
5fPD5pdUcGfF1u7WkFYGEYI3X6oCRZeSUR1x/u0lLHEwWDdKGZaIeyic4gWGEaF4
G55FSt01FszgOnt1Gq/JCtnlblKzp4jjvNWcMkLeXWjxP91R5XzRDbMjx/dM8Ke9
U2+ZRoibtEm5n0FpLVC0A31lm3XfLAGi7H1Ynle/oz8OafTGTbTJ7ipsC354AVQ6
+qg2GQdJj6tyIda9/2TB0u1faqcYbuWO6uz3HuEPO0MS94GdEYktr6UisRsZ6F/Z
OQIDAQAB
-----END PUBLIC KEY-----
```

You can then paste that key into the Janus `Settings` page, in the textarea for your SSH key.

### Creating a project

You may add a new project with the `add_project` command.  The project_name is `snake_cased` and is the primary referrent for the project throughout Etna applications.  Most Etna applications will not acknowledge a project if there is no corresponding Janus project entry.

### Adding permissions

Each user has a permission for a project. You may add a permission using the `permit` command. Each permission consists of a role (`administrator`, `editor`, or `viewer`) and whether or not the user can see `restricted` data.

## Via web portal

Project permissions may be managed for each user by browsing to `JANUS_HOST/project/:project_name`, e.g. https://janus.example.org/project/my_project. This page includes forms to alter roles and privileged data access for each user in your project, disable a user's access to your project, or add a new user to the project with a given role.
