use std::time::Duration;
use reqwest::{Client, Method, RequestBuilder, Url};
use reqwest::header::HeaderMap;
use reqwest_retry::policies::{ExponentialBackoff, ExponentialBackoffBuilder};
use reqwest_retry::RetryTransientMiddleware;

type InnerClient = reqwest_middleware::ClientWithMiddleware;

pub struct EtnaClient {
    host: String,
    headers: HeaderMap,
}

impl EtnaClient {
    pub fn format_url(&self, path: &str) -> String {
        if path.starts_with('/') {
            format!("https://{}{}", self.host, path)
        } else {
            format!("https://{}/{}", self.host, path)
        }
    }

    pub fn set_headers(&mut self, headers: HeaderMap) -> Result<(), anyhow::Error> {
        self.headers.extend(headers);
        self.client = self.rebuild_client()?;
        Ok(())
    }

    pub fn new(host: &str) -> Self {
        EtnaClient {
            host: host.to_string(),
            headers: HeaderMap::new(),
        }
    }

    pub fn client(&self) -> anyhow::Result<InnerClient> {
        let builder = reqwest_middleware::ClientBuilder::new(
            reqwest::Client::builder()
                .default_headers(self.headers.clone())
                .build()?);

        builder.with(
            RetryTransientMiddleware::new_with_policy(
                ExponentialBackoff::builder()
                    .build_with_total_retry_duration(Duration::from_secs(10))
            )
        );

        Ok(builder.build())
    }
}

impl Clone for EtnaClient {
    fn clone(&self) -> Self {
        EtnaClient {
            host: self.host.clone(),
            headers: self.headers.clone(),
        }
    }
}