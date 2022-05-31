use std::path::{Path, PathBuf};
use async_trait::async_trait;
use crate::janus::EtnaAuthManager;
use crate::worker::Work;

pub enum PathNodeType {
    Folder,
    File
}

pub struct PullNode {
    local_path: PathBuf,
    project: String,
    bucket: String,
    remote_path: String,
    assumed_type: PathNodeType,
}

#[async_trait]
impl Work for PullNode {
    async fn do_work(auth_manager: &mut EtnaAuthManager) -> anyhow::Result<()> {
        todo!()
    }
}