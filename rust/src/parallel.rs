//! Row-wise fan-out over a half matrix, with the result independent of how
//! many threads ran.
//!
//! The three all-pairs matrices in this crate -- the k-tuple matrix in `dist`,
//! its length-corrected second pass, and the matrix `disttbfast` recomputes
//! over the aligned rows -- all have the same shape: cell `(i, j)` is a pure
//! function of `i`, `j` and data nobody writes while the matrix is being
//! filled. Nothing accumulates across cells and nothing reads a cell another
//! cell writes. So the only thing a thread count can change is *when* each
//! cell is computed, and the bytes do not depend on that.
//!
//! That is the whole argument for threading these and not the progressive
//! step, which walks persistent scratch and whose node order is load-bearing.
//! Here the determinism is structural, not something a differential has to
//! establish per release -- though one is run anyway, because "structural" is
//! an argument and this repository prefers measurements.
//!
//! The C threads the same loops differently: one mutex-guarded `(i, j)` cursor
//! that every worker pulls the next pair from (`disttbfast.c:1274-1296`). That
//! is dynamic load balancing, and it also means the C's own assignment of
//! pairs to threads differs run to run. It reaches the same answer for the
//! same reason this does. Static round-robin by row is used instead because it
//! needs no lock and no shared cursor, and because a fixed assignment is
//! reproducible under a debugger.
//!
//! No dependency and no `unsafe`: `std::thread::scope` lets the workers borrow
//! the shared inputs, and dealing out `&mut Vec<T>` row handles gives each
//! worker exclusive access to the rows it fills.

/// Fills every row of `mtx` by calling `body(&mut scratch, i, &mut mtx[i])`.
///
/// Each row goes to exactly one worker, so `body` may write its row freely.
/// It must not depend on which worker runs it or in what order: that is the
/// caller's side of the bargain and what keeps the output thread-count
/// independent.
///
/// `scratch` is called once per worker, not once per row, because these loops
/// carry per-worker buffers that the C allocates once per thread too. A worker
/// that gets no rows still builds one; the row counts are only balanced, not
/// equal, and sizing that buffer is cheap next to the matrix.
///
/// `nthreads` of 0 or 1 runs inline on the calling thread with no spawn at all,
/// which is the default path -- the wrapper passes `-C 0-0` unless someone
/// asks for threads.
pub fn fill_rows<T, S>(
    mtx: &mut [Vec<T>],
    nthreads: usize,
    scratch: impl Fn() -> S + Send + Sync,
    body: impl Fn(&mut S, usize, &mut Vec<T>) + Send + Sync,
) where
    T: Send,
{
    if nthreads <= 1 || mtx.len() <= 1 {
        let mut s = scratch();
        for (i, row) in mtx.iter_mut().enumerate() {
            body(&mut s, i, row);
        }
        return;
    }

    // Round-robin rather than contiguous blocks. These are half matrices, so
    // row `i` holds `njob - i` cells and a block split would hand the first
    // worker nearly all the work. Dealing every nth row gives each worker the
    // same spread of row lengths.
    let nthreads = nthreads.min(mtx.len());
    let mut buckets: Vec<Vec<(usize, &mut Vec<T>)>> =
        (0..nthreads).map(|_| Vec::new()).collect();
    for (i, row) in mtx.iter_mut().enumerate() {
        buckets[i % nthreads].push((i, row));
    }

    // Rebound as references so each worker's `move` closure captures the
    // reference and not the closure itself.
    let scratch = &scratch;
    let body = &body;
    std::thread::scope(|sc| {
        for bucket in buckets {
            sc.spawn(move || {
                let mut s = scratch();
                for (i, row) in bucket {
                    body(&mut s, i, row);
                }
            });
        }
    });
}

/// How many workers to run, given what `-C` asked for.
///
/// `-C 0-0` -- what the wrapper passes with no `--thread` -- means sequential,
/// and so does any negative value. Everything else is taken at face value:
/// MAFFT's own `--thread -1` resolves to a core count in the wrapper
/// (`mafft.tmpl:1210-1246`) long before a component sees it, so a component
/// that second-guesses the number would be overriding a decision already made.
pub fn workers(nthread: i32) -> usize {
    if nthread <= 0 {
        1
    } else {
        nthread as usize
    }
}
