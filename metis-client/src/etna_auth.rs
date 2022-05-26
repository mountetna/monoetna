use std::borrow::Borrow;
use std::collections::{BTreeMap, HashMap};
use std::io::Bytes;
use std::time::Duration;
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
use crate::etna_client::EtnaClient;

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

#[derive(Serialize, Deserialize)]
struct TaskTokenRequest {
    token_type: String,
    project_name: String,
    read_only: bool,
}

async fn get_nonce_auth(
    janus: &EtnaClient,
    private_key_pem: &str,
    email: &str
) -> Result<HeaderMap, anyhow::Error> {
    let client = janus.client();
    let response = client.get(janus.format_url("/api/tokens/nonce")).send().await?;
    let nonce = response.text().await?;

    let (digest, sig) = sign_it(&nonce, private_key_pem, email)
        .map_err(|rsa_err| { anyhow!(rsa_err) })?;

    auth_headers("Signed-Nonce", format!("{}.{}", digest, sig).as_str())
}

async fn get_task_token_auth(
    janus: &EtnaClient,
    project_name: &str,
) -> Result<HeaderMap, anyhow::Error> {
    let client = janus.client();
    let response = client.post(janus.format_url("/api/tokens/generate"))
        .json(&TaskTokenRequest {
            token_type: "token".to_string(),
            project_name: project_name.to_string(),
            read_only: false,
        }).send().await?;
    let token = response.text().await?;
    auth_headers("Etna", token.as_str())
}

fn sign_it(
    nonce: &str,
    private_key_pem: &str,
    email: &str
) -> Result<(String, String), rsa::errors::Error> {
    let key = RsaPrivateKey::from_pkcs1_pem(private_key_pem)?;
    let digest = format!("{}.{}", nonce, base64::encode(email));
    let sig = key.sign(PaddingScheme::PKCS1v15Sign { hash: Some(Hash::SHA2_256)},
             digest.as_bytes())?;
    Ok((
        digest,
        base64::encode(sig),
    ))
}

pub fn auth_headers(lead: &str, value: &str) -> Result<HeaderMap, anyhow::Error> {
    let mut map = HeaderMap::new();
    map.insert(
        HeaderName::from_static("Authorization"),
        HeaderValue::from_str(&format!("{} {}", lead, value))
            .map_err(|e| anywho!(e))?
    );
    Ok(map)
}

pub struct Claims {
    perm: String,
    pub email: String,
    pub name: String,
    flags: Option<String>,
    pub task: Option<bool>,
}

fn select_public_key(public_key_override: &Option<&str>) -> Result<RS256PublicKey, anyhow::Error> {
    RS256PublicKey::from_pem(
        public_key_override.unwrap_or(&janus_public_key)
    )
}

fn parse_token<T: RSAPublicKeyLike>(token: &str, public_key: &T) -> Result<Claims, anyhow::Error> {
    public_key.verify_token::<Claims>(token)
}

pub fn retry_up_to_secs(secs: u64) -> Box<dyn Middleware> {
    Box::new(RetryTransientMiddleware::new_with_policy(
        ExponentialBackoff::builder()
            .build_with_total_retry_duration(Duration::from_secs(secs))
    ))
}
