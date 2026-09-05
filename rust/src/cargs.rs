//! MAFFT's `arguments()` option loop, transliterated once.
//!
//! Every C binary in `core/` repeats the same loop, and it has three details
//! that a natural rewrite gets wrong -- each of which changes observable
//! behaviour, so they are the reason this exists rather than a hand-rolled
//! parser per binary:
//!
//! ```c
//! while( --argc > 0 && (*++argv)[0] == '-' )
//! {
//!     while ( (c = *++argv[0]) )
//!     {
//!         switch( c )
//!         {
//!             case 'i': inputfile = *++argv; --argc; goto nextoption;
//!             default:  fprintf( stderr, "illegal option %c\n", c );
//!                       argc = 0;
//!                       break;
//!         }
//!     }
//!     nextoption: ;
//! }
//! if( argc != 0 ) { fprintf( stderr, "options: Check source file !\n" ); exit( 1 ); }
//! ```
//!
//! 1. **`break` in `default:` leaves the switch, not the scan.** An unknown
//!    option does not stop the loop -- the rest of the cluster is still walked
//!    and each unknown character is named on its own line.
//! 2. **`argc = 0` is not "stop cleanly".** The loop condition's own `--argc`
//!    then takes it to -1, which is what the check below the loop turns into
//!    `exit(1)`. An unknown option is fatal, but only after it is reported.
//! 3. **`argv` advances inside the condition**, so it advances on the
//!    iteration whose `'-'` test *fails* too. That leaves `argc` non-zero and
//!    makes a positional argument fatal rather than silently ignored -- which
//!    is why `mafft input.fa` is an error and `mafft -i input.fa` is not.

/// What a binary's option loop knows about one option character.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Opt {
    /// A bare flag: `-D`, `-Y`.
    Flag,
    /// Consumes the following argument: `-i file`, `-o 3`.
    Value,
}

/// One option as the loop saw it, in the order the loop saw them.
pub struct Parsed {
    /// The option character.
    pub c: char,
    /// The following argument, for [`Opt::Value`] options. `None` when the
    /// option was the last thing on the command line -- the C reads one past
    /// the end there, and the caller decides what that means for it.
    pub value: Option<String>,
}

/// Walk MAFFT's option loop and return the options it recognised.
///
/// Diagnostics and the exit go through here rather than through the caller, so
/// every binary reports an unknown option identically to the C.
///
/// Args:
///     args: The arguments with the program name already removed.
///     classify: What each option character is, or `None` if the binary does
///         not know it.
///
/// Returns:
///     The recognised options in command-line order. Does not return at all
///     when the C would have exited.
pub fn parse(args: &[String], classify: impl Fn(char) -> Option<Opt>) -> Vec<Parsed> {
    let mut argc = args.len() as i64 + 1;
    let mut p: i64 = -1;
    let mut out = Vec::new();
    loop {
        argc -= 1;
        if argc <= 0 {
            break;
        }
        p += 1;
        let Some(arg) = args.get(p as usize) else {
            break;
        };
        if !arg.starts_with('-') {
            break;
        }
        let cluster: Vec<char> = arg.chars().skip(1).collect();
        let mut k = 0usize;
        while k < cluster.len() {
            let c = cluster[k];
            match classify(c) {
                Some(Opt::Value) => {
                    p += 1;
                    out.push(Parsed {
                        c,
                        value: args.get(p as usize).cloned(),
                    });
                    argc -= 1;
                    // `goto nextoption`: the rest of the cluster is abandoned.
                    break;
                }
                Some(Opt::Flag) => {
                    out.push(Parsed { c, value: None });
                    k += 1;
                }
                None => {
                    eprintln!("illegal option {c}");
                    argc = 0;
                    k += 1;
                }
            }
        }
    }
    if argc != 0 {
        eprintln!("options: Check source file !");
        std::process::exit(1);
    }
    out
}

/// `myatoi()`: `atoi` behind a NULL check.
///
/// The check is the whole reason the C has this wrapper -- an option given as
/// the last word on the command line hands it the pointer one past the end,
/// and it exits rather than dereferencing. `None` is that case.
///
/// Args:
///     s: The argument text, or `None` when the option had no following
///         argument.
///
/// Returns:
///     The leading integer, or 0 for text that does not start with one. Does
///     not return for `None`.
pub fn myatoi(s: Option<&str>) -> i64 {
    let Some(s) = s else {
        eprintln!("Error in myatoi()");
        std::process::exit(1);
    };
    let b = s.as_bytes();
    let mut i = 0usize;
    let mut sign = 1i64;
    while i < b.len() && (b[i] == b' ' || b[i] == b'\t' || b[i] == b'\n') {
        i += 1;
    }
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        if b[i] == b'-' {
            sign = -1;
        }
        i += 1;
    }
    let mut v = 0i64;
    while i < b.len() && b[i].is_ascii_digit() {
        v = v * 10 + i64::from(b[i] - b'0');
        i += 1;
    }
    sign * v
}

#[cfg(test)]
mod tests {
    use super::*;

    fn classify_countlen(c: char) -> Option<Opt> {
        match c {
            'i' => Some(Opt::Value),
            _ => None,
        }
    }

    #[test]
    fn a_value_option_takes_the_next_argument() {
        let args = vec!["-i".to_string(), "in.fa".to_string()];
        let got = parse(&args, classify_countlen);
        assert_eq!(got.len(), 1);
        assert_eq!(got[0].c, 'i');
        assert_eq!(got[0].value.as_deref(), Some("in.fa"));
    }

    #[test]
    fn no_arguments_yields_nothing_and_does_not_exit() {
        assert!(parse(&[], classify_countlen).is_empty());
    }

    #[test]
    fn a_cluster_is_walked_left_to_right_until_a_value_option() {
        // `-DP` sets both; `-Di file` stops the cluster at `i`.
        let classify = |c: char| match c {
            'D' | 'P' => Some(Opt::Flag),
            'i' => Some(Opt::Value),
            _ => None,
        };
        let args = vec!["-DP".to_string()];
        let got = parse(&args, classify);
        assert_eq!(got.iter().map(|o| o.c).collect::<Vec<_>>(), vec!['D', 'P']);

        let args = vec!["-Di".to_string(), "in.fa".to_string()];
        let got = parse(&args, classify);
        assert_eq!(got.iter().map(|o| o.c).collect::<Vec<_>>(), vec!['D', 'i']);
        assert_eq!(got[1].value.as_deref(), Some("in.fa"));
    }

    // A value option in last position is deliberately NOT unit-tested here.
    // The C's loop leaves `argc` at -1 in that case and exits, and so does
    // this one -- `std::process::exit` would take the whole test binary with
    // it, so the harness cannot observe it. That path is covered where it can
    // be, as a subprocess: `verify_utils.sh`'s "countlen -i (no arg)" and
    // "f2cl -l 0" cases compare the message and the status against the C.

    #[test]
    fn myatoi_stops_at_the_first_non_digit() {
        assert_eq!(myatoi(Some("12")), 12);
        assert_eq!(myatoi(Some("12abc")), 12);
        assert_eq!(myatoi(Some("-7")), -7);
        assert_eq!(myatoi(Some("abc")), 0);
        assert_eq!(myatoi(Some("")), 0);
        assert_eq!(myatoi(Some("  34")), 34);
    }
}
