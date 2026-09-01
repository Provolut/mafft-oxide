#!/bin/bash
# Upstream's own test/script, but per-case and reporting instead of exit-on-first.
BUILD="$1"; T="$2/test"
cd "$T" || exit 1
export PATH="$BUILD/bin:$PATH"
pass=0; fail=0; failed=""
run() { # name, expected, args...
  local name="$1" exp="$2"; shift 2
  mafft "$@" sample > /tmp/mt.$$ 2>/dev/null
  if diff -q /tmp/mt.$$ "$exp" >/dev/null 2>&1; then pass=$((pass+1)); else fail=$((fail+1)); failed="$failed $name"; fi
}
run fftns2            sample.fftns2
run fftnsi            sample.fftnsi            --maxiterate 100
run dpparttree        sample.dpparttree        --dpparttree
run gins1             sample.gins1             --globalpair
run ginsi             sample.ginsi             --globalpair --maxiterate 100
run ginsi.allowshift  sample.ginsi.allowshift  --allowshift --globalpair --maxiterate 100
run lins1             sample.lins1             --localpair
run linsi             sample.linsi             --localpair --maxiterate 100
run parttree          sample.parttree          --parttree
rm -f /tmp/mt.$$
echo "pass=$pass fail=$fail${failed:+  FAILED:$failed}"
