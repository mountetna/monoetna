
use async_trait::async_trait;
use crate::janus::EtnaAuthManager;

#[async_trait]
pub trait Work {
    // Future: Abstract out an 'etna work context', right now that is just the auth manager
    // itself.
    async fn do_work(auth_manager: &mut EtnaAuthManager) -> anyhow::Result<()>;
}