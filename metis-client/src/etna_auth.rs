use std::borrow::Borrow;
use std::collections::BTreeMap;
use std::io::Bytes;
use std::time::Duration;
use jwt_simple::prelude::*;

use reqwest::{Request, Response};
use reqwest::header::{HeaderName, HeaderValue};

use reqwest_middleware::{Middleware, Next};
use reqwest_retry::policies::ExponentialBackoff;
use reqwest_retry::RetryTransientMiddleware;
use async_trait::async_trait;
use anyhow::{anyhow, Error};
use jwt_simple::JWTError;
use jwt_simple::prelude::RS256PublicKey;
use lazy_static::lazy_static;
use task_local_extensions::Extensions;

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

pub struct EtnaTokenAuth {
    token_header: HeaderValue,
    pub token: String,
}

#[async_trait]
impl Middleware for EtnaTokenAuth {
    async fn handle(
        &self,
        req: Request,
        extensions: &mut Extensions,
        next: Next<'_>,
    ) -> reqwest_middleware::Result<Response> {
        let mut new_req = req.try_clone().ok_or(
        reqwest_middleware::Error::Middleware(anyhow!(
                "Could not clone request to add new headers".to_string()
            ))
        )?;
        // HeaderValue::from_bytes()
        new_req.headers_mut().insert(
            HeaderName::from_static("Authorization"),
            self.token_header.clone()
        );
        next.run(new_req, extensions).await
    }
}

pub fn apply_auth(token: &Option<&str>, key: &Option<&[u8]>) {

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
