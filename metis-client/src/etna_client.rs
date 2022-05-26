use std::time::Duration;
use reqwest::{Method, RequestBuilder, Url};
use reqwest::header::HeaderMap;
use reqwest_retry::policies::{ExponentialBackoff, ExponentialBackoffBuilder};
use reqwest_retry::RetryTransientMiddleware;

type InnerClient = reqwest_middleware::ClientWithMiddleware;

pub struct EtnaClient {
    host: String,
    client: InnerClient,
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

    pub fn client(&'_ self) -> &'_ InnerClient {
       &self.client
    }

    pub fn set_headers(&mut self, headers: HeaderMap) -> Result<(), anyhow::Error> {
        self.headers.extend(headers);
        self.client = self.rebuild_client()?;
        Ok(())
    }

    fn rebuild_client(&self) -> Result<InnerClient, anyhow::Error> {
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
            client: self.client.clone(),
            headers: self.headers.clone(),
        }
    }
}