use std::borrow::Borrow;
use std::path::PathBuf;
use std::str::FromStr;
use anyhow::anyhow;
use crate::etna::{EtnaClient, EtnaEndpoint, EtnaServices};
use crate::etna::EtnaTransaction;
use async_trait::async_trait;
use serde::{Deserialize, Serialize};
use task_local_extensions::Extensions;

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

pub struct DownloadChunk {
    download_url: String,
    chunk_start: u64,
    chunk_length: u64,
    expected_md5: Option<String>,
}

#[derive(Debug,Deserialize)]
pub struct File {
    file_name: String,
    file_hash: String,
    updated_at: String,
    file_path: Option<String>,
    folder_id: Option<u64>,
    project_name: String,
    bucket_name: String,
    download_url: Option<String>,
}

#[derive(Debug,Deserialize)]
pub struct Folder {
    id: u64,
    folder_path: String,
    project_name: String,
    bucket_name: String,
    folder_name: String,
    updated_at: String,
    folder_id: u64,
}

#[derive(Debug,Deserialize)]
pub struct FolderContents {
    folders: Option<Vec<Folder>>,
    files: Option<Vec<File>>,
}

#[async_trait]
impl EtnaTransaction<bytes::Bytes> for DownloadChunk {
    async fn run(&self, client: &EtnaClient, ext: &mut Extensions, services: &EtnaServices) -> anyhow::Result<bytes::Bytes> {
        let response = client.get(&self.download_url)
            .header(
                "Range",
                format!(
                    "bytes={}-{}",
                    self.chunk_start,
                    self.chunk_start + self.chunk_length
                )
            )
            .send_with_extensions(ext).await?.error_for_status()?;

        if let Some(expected_md5) = &self.expected_md5 {
            if let Some(md5) = response.headers().get("Content-MD5") {
                if md5.as_bytes() != expected_md5.as_bytes() {
                    return Err(anyhow!("Server md5 changed during processing! Found {}, expected {}", md5.to_str()?, expected_md5));
                }
            }
        }


        Ok(response.bytes().await?)
    }
}

#[derive(Deserialize,Debug)]
struct AuthorizedDownload {
    download_url: String,
}

#[derive(Clone,Debug,Serialize)]
pub enum MetisPath {
    FilePath(MetisFilePath),
    FolderPath(MetisFolderPath),
}

#[derive(Clone,Debug,Serialize)]
pub struct MetisFilePath {
    project_name: String,
    bucket_name: String,
    file_path: String,
}

#[derive(Clone,Debug,Serialize)]
pub struct MetisFolderPath {
    project_name: String,
    bucket_name: String,
    folder_path: String,
}

pub struct AuthorizeDownload(pub MetisFilePath);
pub struct ListFolder(pub MetisFolderPath);

#[derive(Clone,Debug)]
pub enum PathNodeType {
    Folder,
    File
}

pub struct DownloadMeta {
    pub download_url: String,
    pub md5: Option<String>,
    pub size: u64,
}

#[async_trait]
impl EtnaTransaction<FolderContents> for ListFolder {
    async fn run(&self, client: &EtnaClient, extensions: &mut Extensions, services: &EtnaServices) -> anyhow::Result<FolderContents> {
        let req = if self.0.folder_path.is_empty() {
            client.get(services.metis.format_url(&format!("{}/list/{}", self.0.project_name, self.0.bucket_name)))
        } else {
            client.get(services.metis.format_url(&format!("{}/list/{}/{}", self.0.project_name, self.0.bucket_name, self.0.folder_path)))
        };

        let response = req.send_with_extensions(extensions).await?;
        Ok(response.json::<FolderContents>().await?)
    }
}

#[async_trait]
impl EtnaTransaction<DownloadMeta> for AuthorizeDownload {
    async fn run(&self, client: &EtnaClient, ext: &mut Extensions, services: &EtnaServices) -> anyhow::Result<DownloadMeta> {
        let url = services.metis.format_url("/authorize/download");
        let response = client.post(url).json(&self.0).send_with_extensions(ext).await?;
        let authorized = response.json::<AuthorizedDownload>().await?;

        let response = client.get(&authorized.download_url)
            .header(
               "Range",
               "bytes=0-0",
            )
            .send_with_extensions(ext).await?.error_for_status()?;

        let mut result = DownloadMeta {
            download_url: authorized.download_url,
            size: 0,
            md5: None,
        };

        if let Some(v) = response.headers().get("Content-MD5") {
            result.md5 = Some(v.to_str()?.to_string());
        }

        if let Some(v) = response.headers().get("Content-Length") {
            result.size = u64::from_str(v.to_str()?)?
        } else {
            return Err(anyhow!("No valid Content-Length found in metis download headers!"))
        }

        Ok(result)
    }
}