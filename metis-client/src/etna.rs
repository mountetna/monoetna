use std::borrow::Borrow;
use std::ops::Deref;
use std::sync::{Arc, Mutex};
use std::time::Duration;
use reqwest::{Client, Method, RequestBuilder, Url};
use reqwest::header::HeaderMap;
use reqwest_retry::policies::{ExponentialBackoff, ExponentialBackoffBuilder};
use reqwest_retry::RetryTransientMiddleware;
use async_trait::async_trait;

type InnerClient = reqwest_middleware::ClientWithMiddleware;

#[async_trait]
pub trait EtnaTransaction<T = ()> {
    async fn run(&self, services: &EtnaServices) -> anyhow::Result<T>;
}

pub struct EtnaEndpoint {
    host: String,
}

impl EtnaEndpoint {
    pub fn for_host(host: &str, headers: Arc<HeaderMap>) -> Self {
        EtnaEndpoint {
            host: host.to_string(),
        }
    }

    pub fn for_service(service: &str, headers: Arc<HeaderMap>) -> Self {
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

impl Clone for EtnaEndpoint {
    fn clone(&self) -> Self {
        EtnaEndpoint {
            host: self.host.clone(),
        }
    }
}

pub struct EtnaServices {
    pub janus: EtnaEndpoint,
    pub metis: EtnaEndpoint,
    pub magma: EtnaEndpoint,
    headers: Arc<HeaderMap>,
}

impl Clone for EtnaServices {
    fn clone(&self) -> Self {
        EtnaServices {
            janus: self.janus.clone(),
            metis: self.metis.clone(),
            magma: self.magma.clone(),
            headers: self.headers.clone(),
        }
    }
}

impl EtnaServices {
    pub fn with_headers(&self, headers: Arc<HeaderMap>) -> Self {
        EtnaServices {
            headers,
            ..self.clone()
        }
    }

    pub fn client(&self) -> anyhow::Result<InnerClient> {
        let builder = reqwest_middleware::ClientBuilder::new(
            reqwest::Client::builder()
                .default_headers(self.headers.deref().clone())
                .build()?);

        let builder = builder.with(
            RetryTransientMiddleware::new_with_policy(
                ExponentialBackoff::builder()
                    .build_with_total_retry_duration(Duration::from_secs(60 * 5))
            )
        );

        Ok(builder.build())
    }
}

impl Default for EtnaServices {
    fn default() -> Self {
        let headers = Arc::new(HeaderMap::default());

        EtnaServices {
            janus: EtnaEndpoint::for_service("janus", headers.clone()),
            metis: EtnaEndpoint::for_service("metis", headers.clone()),
            magma: EtnaEndpoint::for_service("magma", headers.clone()),
            headers,
        }
    }
}