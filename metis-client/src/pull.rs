use std::borrow::Borrow;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use async_trait::async_trait;
use crate::etna::{EtnaServices, EtnaTransaction};
use crate::janus::EtnaAuthContext;
use crate::metis::{AuthorizeDownload, DownloadMeta, FolderContents, ListFolder, MetisFilePath, MetisPath};
use crate::pull::PlannedPullWork::DownloadFile;
use crate::worker::{EtnaWork, EtnaWorker};
use md5::{Md5, Digest};

#[derive(Clone,Debug)]
pub struct PullNode {
    local_path: PathBuf,
    metis_path: MetisPath,
}

pub enum PlannedPullWork {
    RemoveLocalPath(PathBuf),
    DownloadFile(MetisFilePath, PathBuf),
    MakeLocalFolder(PathBuf),
}

#[async_trait]
impl EtnaWork for PullNode {
    async fn do_work(&self, worker: &mut EtnaWorker) -> anyhow::Result<Vec<&dyn EtnaWork>> {
        todo!()
    }
}

async fn plan_pull(worker: &mut EtnaWorker, pull: &PullNode) -> anyhow::Result<Vec<PlannedPullWork>> {
    let mut result: Vec<PlannedPullWork> = Vec::new();

    match &pull.metis_path {
        MetisPath::FilePath(path) => {
            let md = fs::metadata(&pull.local_path).ok();
            if md.as_ref().map(|r|r.is_dir()).unwrap_or(false) {
                result.push(PlannedPullWork::RemoveLocalPath(pull.local_path.clone()))
            }

            if md.as_ref().map(|r|r.is_file()).unwrap_or(false) {
                let meta = worker.execute::<DownloadMeta>(&AuthorizeDownload(path.clone())).await?;
                if let Some(remote_md5) = meta.md5 {
                    let mut buff: Vec<u8> = Vec::with_capacity(1024 * 1024);
                    let mut hasher = Md5::new();
                    let mut fd = fs::File::open(&pull.local_path)?;

                    loop {
                        let read_bytes = fd.read(buff.as_mut_slice())?;
                        if read_bytes < buff.len() {
                            buff.truncate(read_bytes);
                            hasher.update(buff.as_slice());
                            break
                        } else {
                            hasher.update(buff.as_slice());
                        }
                    }

                    let digest = hasher.finalize();
                    let hex_result = format!("{:02x?}", digest);

                    if hex_result != remote_md5 {
                        result.push(PlannedPullWork::DownloadFile(path.clone(), pull.local_path.clone()))
                    }
                }
            } else {
                result.push(PlannedPullWork::DownloadFile(path.clone(), pull.local_path.clone()))
            }
        },
        MetisPath::FolderPath(path) => {
            let list = ListFolder(path.clone());
            let folder_contents = worker.execute::<FolderContents>(&list).await?;
        }
    }

    Ok(result)
}