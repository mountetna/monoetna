use std::alloc::System;
use std::borrow::Borrow;
use std::collections::{BTreeMap, HashMap};
use std::io::Bytes;
use std::path::PathBuf;
use std::sync::Arc;
use std::time::{Duration, SystemTime, UNIX_EPOCH};
use jwt_simple::prelude::*;
use rsa::{Hash, PaddingScheme, RsaPrivateKey};
use serde::{Deserialize, Serialize};

use reqwest::{Request, Response};
use reqwest::header::{HeaderMap, HeaderName, HeaderValue, InvalidHeaderValue};

use reqwest_middleware::{Middleware, Next};
use reqwest_retry::policies::ExponentialBackoff;
use reqwest_retry::RetryTransientMiddleware;
use async_trait::async_trait;
use anyhow::{anyhow, Error};
use jwt_simple::JWTError;
use jwt_simple::prelude::RS256PublicKey;
use lazy_static::lazy_static;
use rsa::pkcs1::FromRsaPrivateKey;
use task_local_extensions::Extensions;
use crate::etna::{EtnaEndpoint, EtnaServices, EtnaTransaction};

lazy_static! {
    static ref janus_public_key: &'static str = "-----BEGIN PUBLIC KEY-----
MIICIjANBgkqhkiG9w0BAQEFAAOCAg8AMIICCgKCAgEA3cmYfKvOUxKLgv/TQ9aG
hMA4TCs2gQosmYejdE8I4Y39dSHtUnERX3ztMoesZVdhXU1GR7LWTfuOIze/H9q+
x/IDx/KBBoiYV+oZ6tCRVV74lNJJpAEyRT5hAW9yW2/OTG2R6wlXc1oxjRfi+M38
mMzNN/mR54Yzgr4Rch1NYCd/ZOt0zrBKuNg7xu17N1wDjc7XDVgCn4btgpwO3g8D
wsiAXwnpbXgi1FP+7oiHiHPS07+CuKDafvLXqAVB2VymIR5muyWyayrxGjxyulZY
TiDEPew/1TSBS+F8c4ZmAuH7fj+XX4ofWboZFrKl3R8V/041TdZX5nfOiLL3guDT
V1fcc4gKfED4wI8bui0KDZqSTNYXsHNEm2r6y5VRhRwpRqFYT6YvAtiHrS9/U0L8
rl8qoQWGG+bqepaBcQZ6GuPzWcTQfookI0LMQHpexAfBef1IQS7JJ7CrRgLD8ryS
eolIQsOIX1fVTPUlY2x7liDaC49dnoZNDjMTy6bOiWGXzek0iF7qZd12nbDW5vUQ
iJ3tahGjB3AbekXFGXNBOByW/qZrj0EHWWJMysZ2lDlTnAPtBmwwcF9YTZQA7XGO
sPDtWvOVFzFxTEIPok1Cv/4dZDcVX4raeLaftV4PNY38o19E+ensWKeFUoA0a+wC
VnOJnAt45Hq/mNZEBP2I7/MCAwEAAQ==
-----END PUBLIC KEY-----";
}

#[derive(Default)]
struct GetNonce {}

#[async_trait]
impl EtnaTransaction<String> for GetNonce {
    async fn run(&self, services: &EtnaServices) -> anyhow::Result<String> {
        let client = services.client()?;
        let response = client.get(services.janus.format_url("/api/tokens/nonce")).send().await?.error_for_status()?;
        Ok(response.text().await?)
    }
}

#[derive(Serialize, Deserialize)]
struct GenerateTaskToken {
    token_type: String,
    project_name: String,
    read_only: bool,
}

impl Default for GenerateTaskToken {
    fn default() -> Self {
        GenerateTaskToken {
            token_type: "token".to_string(),
            project_name: "".to_string(),
            read_only: false,
        }
    }
}

#[async_trait]
impl EtnaTransaction<String> for GenerateTaskToken {
    async fn run(&self, services: &EtnaServices) -> anyhow::Result<String> {
        let client = services.client()?;
        let response = client.post(services.janus.format_url("/api/tokens/generate"))
            .json(&self).send().await?.error_for_status()?;

        Ok(response.text().await?)
    }
}

pub fn sign_nonce(
    nonce: &str,
    private_key_pem: &str,
    email: &str
) -> anyhow::Result<HeaderMap> {
    let key = RsaPrivateKey::from_pkcs1_pem(private_key_pem)?;
    let digest = format!("{}.{}", nonce, base64::encode(email));
    let sig = base64::encode(key.sign(PaddingScheme::PKCS1v15Sign { hash: Some(Hash::SHA2_256)},
             digest.as_bytes())?);
    auth_headers("Signed-Nonce", format!("{}.{}", digest, sig).as_str())
}

pub fn auth_headers(lead: &str, value: &str) -> anyhow::Result<HeaderMap> {
    let mut map = HeaderMap::new();
    map.insert(
        HeaderName::from_static("Authorization"),
        HeaderValue::from_str(&format!("{} {}", lead, value))?
    );
    Ok(map)
}

#[derive(Serialize, Deserialize)]
pub struct Claims {
    perm: String,
    pub email: String,
    pub name: String,
    flags: Option<String>,
    pub task: Option<bool>,
}


fn select_public_key(public_key_override: &Option<String>) -> anyhow::Result<RS256PublicKey> {
    if let Some(public_key) = public_key_override {
        RS256PublicKey::from_pem(public_key)
    } else {
        RS256PublicKey::from_pem(&janus_public_key)
    }
}

fn parse_token<T: RSAPublicKeyLike>(token: &str, public_key: &T) -> anyhow::Result<JWTClaims<Claims>> {
    public_key.verify_token::<Claims>(token, None)
}

#[derive(Default)]
pub struct AuthConfig {
    email: String,
    private_key: String
}

impl AuthConfig {
    async fn rotate_token(&self, services: &EtnaServices, project_name: &str) -> anyhow::Result<String> {
        let nonce = GetNonce::default().run(services).await?;
        let signed_auth = sign_nonce(&nonce, &self.private_key, &self.email)?;

        let task = GenerateTaskToken {
            project_name: project_name.to_string(),
            ..GenerateTaskToken::default()
        };

        task.run(&services.with_headers(Arc::new(signed_auth))).await
    }
}

pub struct EtnaAuthManager {
    project_name: String,
    services: EtnaServices,
    cur_token: Option<String>,
    janus_pub_key_override: Option<String>,
    auth_config: AuthConfig,
}

impl EtnaAuthManager {
    fn new(
        project_name: String,
        services: EtnaServices,
        auth_config: AuthConfig,
        janus_pub_key_override: Option<String>,
    ) -> Self {
        EtnaAuthManager {
            project_name,
            services,
            janus_pub_key_override,
            auth_config,
            cur_token: None
        }
    }

    fn token_needs_rotation(&self, cur_token: &str) -> anyhow::Result<bool> {
        let key = select_public_key(&self.janus_pub_key_override)?;
        parse_token(cur_token, &key).map(|claims| {
            SystemTime::now().duration_since(UNIX_EPOCH).ok()
                .and_then(|now| claims.expires_at.map(|exp| exp.as_secs() - now.as_secs() < 60))
                .unwrap_or(true)
        })
    }

    async fn run<T>(&mut self, transaction: &dyn EtnaTransaction<T>) -> anyhow::Result<T> {
        if let Some(cur_token) = &self.cur_token {
            if self.token_needs_rotation(cur_token)? {
                self.cur_token = None
            }
        }

        if self.cur_token.is_none() {
            let token = self.auth_config.rotate_token(&self.services, &self.project_name).await?;
            self.services = self.services.with_headers(Arc::new(auth_headers("Etna", &token)?));
            self.cur_token = Some(token);
        };

        transaction.run(&self.services).await
    }
}
