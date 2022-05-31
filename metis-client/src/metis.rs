use std::borrow::Borrow;
use std::path::PathBuf;
use crate::etna::{EtnaEndpoint, EtnaServices};
use crate::etna::EtnaTransaction;
use async_trait::async_trait;

struct Upload {
    project: Option<String>,
    bucket: Option<String>,
    remote_path: Option<String>,
    local_path: Vec<PathBuf>,
}

// #[async_trait]
// impl EtnaTransaction for Upload {
//     async fn run(&self, client: &EtnaEndpoint) -> anyhow::Result<()> {
//     }
// }

struct DownloadChunk {
    download_url: String,
    chunk_start: u64,
    chunk_length: u64,
}

#[async_trait]
impl EtnaTransaction<bytes::Bytes> for DownloadChunk {
    async fn run(&self, services: &EtnaServices) -> anyhow::Result<bytes::Bytes> {
        let client = services.client()?;
        let response = client.get(&self.download_url)
            .header(
                "Range",
                format!(
                    "bytes={}-{}",
                    self.chunk_start,
                    self.chunk_start + self.chunk_length
                )
            )
            .send().await?.error_for_status()?;

        Ok(response.bytes().await?)
    }
}
