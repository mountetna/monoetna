mod etna_auth;
mod etna_client;

use std::path::PathBuf;
use std::time::Duration;
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

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let opts: Cli = Cli::parse();

    let client = reqwest::Client::builder()
        .https_only(true)
        .use_rustls_tls()
        .connect_timeout(Duration::from_secs(6))
        .build()?;

    let env_token = std::env::var("TOKEN").ok();

    match opts.command.unwrap_or(Commands::Push(opts.push)) {
        Commands::Push(push) => {
            ClientBuilder::new(client);
        }
        Commands::Pull(pull) => {
            ClientBuilder::new(client);
        }
    }
    // let resp = reqwest::get("https://httpbin.org/ip")
    //     .await?
    //     .json::<HashMap<String, String>>()
    //     .await?;
    // println!("{:?}", opts);
    // println!("{:#?}", resp);
    Ok(())
}