use std::borrow::Borrow;
use std::ops::Deref;
use std::sync::{Arc, Mutex};
use std::time::Duration;
use reqwest::{Method, RequestBuilder, Url};
use reqwest::header::HeaderMap;
use reqwest_retry::policies::{ExponentialBackoff, ExponentialBackoffBuilder};
use reqwest_retry::RetryTransientMiddleware;
use async_trait::async_trait;
use task_local_extensions::Extensions;
use crate::janus::EtnaAuthState;

pub type EtnaClient = reqwest_middleware::ClientWithMiddleware;

#[async_trait]
pub trait EtnaTransaction<T = ()> {
    async fn run(&self, client: &EtnaClient, extensions: &mut Extensions, services: &EtnaServices) -> anyhow::Result<T>;
}

#[derive(Clone)]
pub struct EtnaEndpoint {
    host: String,
}

impl EtnaEndpoint {
    pub fn for_host(host: &str) -> Self {
        EtnaEndpoint {
            host: host.to_string(),
        }
    }

    pub fn for_service(service: &str) -> Self {
        EtnaEndpoint {
            host: format!("{}.ucsf.edu", service),
        }
    }

    pub fn format_url(&self, path: &str) -> String {
        if path.starts_with('/') {
            format!("https://{}{}", self.host, path)
        } else {
            format!("https://{}/{}", self.host, path)
        }
    }
}

#[derive(Clone)]
pub struct EtnaServices {
    pub janus: Arc<EtnaEndpoint>,
    pub metis: Arc<EtnaEndpoint>,
    pub magma: Arc<EtnaEndpoint>,
}

pub fn base_etna_client_builder(
    headers: HeaderMap,
) -> anyhow::Result<reqwest_middleware::ClientBuilder> {
    let builder = reqwest_middleware::ClientBuilder::new(
        reqwest::Client::builder()
            .default_headers(headers)
            .build()?
    );

    let builder = builder.with(
        RetryTransientMiddleware::new_with_policy(
            ExponentialBackoff::builder()
                .build_with_total_retry_duration(Duration::from_secs(60 * 5))
        )
    );

    Ok(builder)
}

impl Default for EtnaServices {
    fn default() -> Self {
        EtnaServices {
            janus: Arc::new(EtnaEndpoint::for_service("janus")),
            metis: Arc::new(EtnaEndpoint::for_service("metis")),
            magma: Arc::new(EtnaEndpoint::for_service("magma")),
        }
    }
}