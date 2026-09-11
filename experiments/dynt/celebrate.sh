#!/usr/bin/env bash
# celebrate.sh -- play a sound + banner when an exclusion run finishes clean.
#
# Usage:
#   celebrate.sh milestone "S_64 > 10^146"        # a chime for a partial win
#   celebrate.sh certified "S_64 = N (148 digits)" # the FANFARE: full [<L>,N] excluded
#
# Two tiers on purpose. Excluding everything below a PARTIAL bound (e.g. 10^146)
# tightens the interval -> milestone chime. Excluding everything below the
# incumbent N itself proves S_64 = N -> the fanfare. Only the second is "we
# proved 64". macOS `afplay`; falls back to the terminal bell everywhere else.
tier="${1:-milestone}"; msg="${2:-}"
play() { for s in "$@"; do f="/System/Library/Sounds/$s.aiff"; [ -f "$f" ] && { afplay "$f"; return; }; done; printf '\a'; }
case "$tier" in
  certified)
    play Hero Glass; play Hero Glass; play Hero Glass   # triple fanfare
    echo "================================================================"
    echo "  ***  S_64 CERTIFIED  ***  $msg"
    echo "  Every 64-factor Carmichael below the incumbent excluded."
    echo "================================================================" ;;
  *)
    play Glass Hero Funk
    echo ">>> milestone: $msg" ;;
esac
