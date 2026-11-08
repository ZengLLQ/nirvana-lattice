use clap::Parser;
use nirvana_lattice::cli::{Cli, Commands};

#[test]
fn test_cli_parsing_subcommands() {
    let cli = Cli::parse_from(["nirvana-lattice", "info"]);
    assert!(matches!(cli.command, Commands::Info));

    let cli = Cli::parse_from(["nirvana-lattice", "inspect"]);
    assert!(matches!(cli.command, Commands::Inspect));

    let cli = Cli::parse_from(["nirvana-lattice", "ntt-benchmark", "--iterations", "10"]);
    if let Commands::NttBenchmark { iterations } = cli.command {
        assert_eq!(iterations, 10);
    } else {
        panic!("Expected NttBenchmark");
    }

    let cli = Cli::parse_from(["nirvana-lattice", "keygen", "-k", "4", "-l", "3"]);
    if let Commands::Keygen { k, l, out } = cli.command {
        assert_eq!(k, 4);
        assert_eq!(l, 3);
        assert!(out.is_none());
    } else {
        panic!("Expected Keygen");
    }
}

#[test]
fn test_cli_info_execution() {
    let cli = Cli::parse_from(["nirvana-lattice", "info"]);
    assert!(cli.run().is_ok());
}

#[test]
fn test_cli_ntt_benchmark_execution() {
    let cli = Cli::parse_from(["nirvana-lattice", "ntt-benchmark", "--iterations", "20"]);
    assert!(cli.run().is_ok());
}

#[test]
fn test_cli_inspect_execution() {
    let cli = Cli::parse_from(["nirvana-lattice", "inspect"]);
    assert!(cli.run().is_ok());
}

#[test]
fn test_cli_self_verify_execution() {
    let cli = Cli::parse_from(["nirvana-lattice", "verify"]);
    assert!(cli.run().is_ok());
}
