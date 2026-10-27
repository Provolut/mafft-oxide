//! `commonIP`: the traceback matrix the aligners share.
//!
//! The C keeps one of these for the whole process -- `extern TLS int
//! **commonIP` in `dp.h`, grown monotonically by whichever aligner needs it
//! biggest (`Salignmm.c:1422`, `Lalign11.c:426`, `Galign11.c:527` and others
//! all assign the same `commonAlloc1`). The port gave `AlignState` and
//! `LalignState` one each, and on `sample.linsi` massif shows both live at the
//! peak at 1,798,280 B apiece: 14% of that case's heap, twice.
//!
//! Sharing it cannot change a single output value. Neither aligner reads a
//! cell it has not written in the same call -- that is what makes the C's
//! sharing safe too -- and neither re-zeroes the buffer on entry today, so
//! handing over dirty contents is already the behaviour. What changes is that
//! the buffer is sized to the larger of the two demands rather than allocated
//! twice.
//!
//! Lent out rather than borrowed in place: both call sites already do
//! `std::mem::take(&mut self.ijp)` and put it back at the end, so taking it
//! out of the pool for the duration costs nothing and means no borrow is held
//! across a DP. A caller that asks while it is already out gets its own buffer
//! instead of a panic -- there is no such path today (`align.rs` never calls
//! into `lalign.rs` or the reverse), and this way adding one would cost memory
//! rather than correctness.

use std::cell::RefCell;

#[derive(Default)]
struct Pool {
    /// The matrix, or empty while it is lent out.
    mtx: Vec<Vec<i32>>,
    /// Dimensions requested so far, as the C's `commonAlloc1`/`2`: the running
    /// maximum, never shrinking.
    alloc1: usize,
    alloc2: usize,
}

thread_local! {
    static POOL: RefCell<Pool> = RefCell::new(Pool::default());
}

/// The shared matrix, able to hold `orlgth1+10` rows of `orlgth2+10`.
///
/// Returns:
///     The pooled matrix, or a fresh one of the same size if the pool's copy
///     is already out on loan.
pub(crate) fn take(orlgth1: usize, orlgth2: usize) -> Vec<Vec<i32>> {
    POOL.with(|p| {
        let mut p = p.borrow_mut();
        let ll1 = orlgth1.max(p.alloc1);
        let ll2 = orlgth2.max(p.alloc2);
        if p.mtx.is_empty() {
            // First call, or someone else holds it. Record the dimensions
            // either way, so the next grow decision sees the running maximum.
            p.alloc1 = ll1;
            p.alloc2 = ll2;
            return vec![vec![0i32; ll2 + 10]; ll1 + 10];
        }
        if orlgth1 > p.alloc1 || orlgth2 > p.alloc2 {
            // Free the old matrix *before* allocating the replacement.
            // `p.mtx = vec![..]` evaluates its right-hand side first, so the
            // two exist at once and peak RSS pays for both -- and growth
            // happens at exactly the moment the working set is largest, so
            // that overlap lands on the peak rather than beside it. massif
            // caught it on `sample.nwns2`: 2,095,760 B and 1,798,280 B live in
            // the same snapshot, 54% of that case's heap, for one matrix's
            // worth of information.
            //
            // The same overlap was in both `ensure_ijp` implementations before
            // they shared a pool; it is not new here, only visible here.
            p.mtx = Vec::new();
            p.mtx = vec![vec![0i32; ll2 + 10]; ll1 + 10];
            p.alloc1 = ll1;
            p.alloc2 = ll2;
        }
        std::mem::take(&mut p.mtx)
    })
}

/// Hand the matrix back.
///
/// Keeps the larger of what the pool holds and what is returned, so a caller
/// that was given its own buffer because the pooled one was out does not
/// shrink the pool on the way back.
pub(crate) fn put(mtx: Vec<Vec<i32>>) {
    POOL.with(|p| {
        let mut p = p.borrow_mut();
        if mtx.len() >= p.mtx.len() {
            p.alloc1 = mtx.len().saturating_sub(10);
            p.alloc2 = mtx.first().map_or(0, |r| r.len()).saturating_sub(10);
            p.mtx = mtx;
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The pool is thread-local, so each test gets its own and they cannot
    /// interfere -- but within a test the order matters, hence one test.
    #[test]
    fn the_pool_grows_monotonically_and_is_reused() {
        let a = take(100, 200);
        assert_eq!(a.len(), 110);
        assert_eq!(a[0].len(), 210);
        put(a);

        // A smaller request must not shrink it: the C's commonAlloc never
        // decreases, and a shrink here would reallocate on every alternation
        // between a big aligner and a small one.
        let b = take(10, 20);
        assert_eq!(b.len(), 110, "a smaller request reallocated the pool");
        assert_eq!(b[0].len(), 210);
        put(b);

        // A bigger one grows it.
        let c = take(300, 50);
        assert_eq!(c.len(), 310);
        assert_eq!(c[0].len(), 210, "the other dimension must not shrink");
        put(c);
    }

    #[test]
    fn a_second_taker_gets_its_own_buffer_rather_than_an_empty_one() {
        let first = take(50, 60);
        // Nothing returns it, so the pool is empty. A second caller must still
        // get something it can index, not a zero-row Vec that panics.
        let second = take(50, 60);
        assert_eq!(second.len(), 60);
        assert_eq!(second[0].len(), 70);
        put(first);
        put(second);
    }
}
