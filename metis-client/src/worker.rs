use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use anyhow::Context;
use async_trait::async_trait;
use reqwest::header::HeaderMap;
use task_local_extensions::Extensions;
use crate::janus::EtnaAuthContext;
use tempfile::tempfile;
use crate::etna::{base_etna_client_builder, EtnaClient, EtnaServices, EtnaTransaction};

pub struct EtnaWorker {
    client: EtnaClient,
    services: EtnaServices,
    extentions: Extensions,
}

impl EtnaWorker {
    pub fn tmp_file() -> anyhow::Result<File> {
        Ok(tempfile::tempfile()?)
    }

    pub async fn execute<T>(&mut self, trans: &dyn EtnaTransaction<T>) -> anyhow::Result<T> {
        trans.run(&self.client, &mut self.extentions, &self.services).await
    }
}

#[async_trait]
pub trait EtnaWork {
    async fn do_work(&self, context: &mut EtnaWorker) -> anyhow::Result<Vec<&dyn EtnaWork>>;
}