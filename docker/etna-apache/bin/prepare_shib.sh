#!/usr/bin/env bash

set -e

VERIFIED_SHIB_DIR=/opt/shibboleth/
mkdir -p $VERIFIED_SHIB_DIR

cd $VERIFIED_SHIB_DIR
curl -o InCommon-metadata.xml https://mdq.incommon.org/entities

curl -fsSL -o inc-md-cert-mdq.pem https://mdq.incommon.org/certs/inc-md-cert-mdq.pem

FP1="sha256 Fingerprint=60:49:74:D6:1F:E0:D7:F4:D6:3D:6C:8D:B9:8A:85:7E:64:2A:B9:B4:70:E3:E8:5D:D5:4D:66:3D:04:96:F9:00"
FP2=$(openssl x509 -sha256 -noout -fingerprint -in inc-md-cert-mdq.pem)

# Verify signing key
[[ $FP1 == $FP2 ]] || exit 1

# Verify the metadata's signing
xmlsec1 --verify --pubkey-cert-pem inc-md-cert-mdq.pem --id-attr:ID "urn:oasis:names:tc:SAML:2.0:metadata:EntitiesDescriptor" InCommon-metadata.xml
