#!/usr/bin/env bash

set -e

VERIFIED_SHIB_DIR=/opt/shibboleth/
mkdir -p $VERIFIED_SHIB_DIR

cd $VERIFIED_SHIB_DIR
curl -o InCommon-metadata.xml https://md.incommon.org/InCommon/InCommon-metadata.xml
curl -o inc-md-cert.pem https://md.incommon.org/certs/inc-md-cert.pem

# Verify the metadata's signing
xmlsec1 --verify --pubkey-cert-pem inc-md-cert.pem --id-attr:ID "urn:oasis:names:tc:SAML:2.0:metadata:EntitiesDescriptor" InCommon-metadata.xml




