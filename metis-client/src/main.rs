mod etna_auth;
mod etna_client;

use std::fmt::Display;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Duration;
use anyhow::anyhow;
use clap::{Parser, Subcommand, Args};
use reqwest_middleware::ClientBuilder;

#[derive(Debug, Parser)]
#[clap(name = "metis-client")]
#[clap(about = "A robust metis syncing client", long_about = None)]
#[clap(args_conflicts_with_subcommands = true)]
struct Cli {
    #[clap(subcommand)]
    command: Option<Commands>,

    #[clap(flatten)]
    push: Push,
}

#[derive(Debug, Subcommand)]
enum Commands {
    Push(Push),
    Pull(Pull),
}

#[derive(Debug, Args)]
struct Push {
    #[clap(long)]
    project: Option<String>,
    #[clap(long)]
    bucket: Option<String>,
    #[clap(long)]
    remote_path: Option<String>,
    #[clap(long, parse(from_os_str))]
    rsa_file: Option<Vec<PathBuf>>,

    #[clap(long, parse(from_os_str))]
    path: Option<Vec<PathBuf>>,
}

#[derive(Debug, Args)]
struct Pull {
    #[clap(long)]
    project: Option<String>,
    #[clap(long)]
    bucket: Option<String>,
    #[clap(long)]
    remote_path: Option<String>,
    #[clap(long, parse(from_os_str))]
    rsa_file: Option<Vec<PathBuf>>,

    #[clap(required = true, parse(from_os_str))]
    path: Vec<PathBuf>,
}

fn default_key_path() -> Option<PathBuf> {
    if let Some(user_dirs) = directories::UserDirs::new() {
        return Some(user_dirs.home_dir().join(".ssh").join(""));
    }

    None
}

fn try_rsa_pem(file_override: &Option<Vec<PathBuf>>) -> anyhow::Result<Option<String>> {
    if let Some(file_override) = file_override {
        return Ok(Some(fs::read_to_string(file_override)?));
    }

    if let Some(key_path) = default_key_path() {
        if let Ok(contents) = fs::read_to_string(key_path) {
            return Ok(Some(contents))
        }
    }

    Ok(None)
}

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    let opts: Cli = Cli::parse();

    let env_token = std::env::var("TOKEN").ok();

    match opts.command.unwrap_or(Commands::Push(opts.push)) {
        Commands::Push(push) => {
            let private_key = try_rsa_pem(&push.rsa_file)?;
        }
        Commands::Pull(pull) => {
            let private_key = try_rsa_pem(&pull.rsa_file)?;
        }
    }

    Ok(())
}