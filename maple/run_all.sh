#!/bin/sh
# Run every Maple example and print its checks.
#
#   sh run_all.sh                 # uses the `maple` found in PATH
#   MAPLE=/opt/maple/bin/maple sh run_all.sh
#
# Each script prints one [pass]/[FAIL] line per identity and a final
# ALL CHECKS PASSED / SOME CHECKS FAILED.  The exit status is 0; read the
# output.  A recorded run is in transcript.txt.
#
# The names sort in the order of the paper: ex_3_4, ex_5_6, ex_6_3, ex_6_5,
# ex_6_7, then the algorithm demos alg_2, alg_3, alg_diagnostics.

cd "$(dirname "$0")" || exit 1
MAPLE=${MAPLE:-maple}

"$MAPLE" -q <<'EOF'
printf("%s\n\n", kernelopts(version));
EOF

for f in ex_*.mpl alg_*.mpl; do
    printf '\n======================================================================\n'
    printf '%s\n' "$f"
    printf '======================================================================\n'
    "$MAPLE" -q < "$f"
done
