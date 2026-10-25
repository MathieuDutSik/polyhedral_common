#!/bin/bash
# Run a command under a memory cap and a wall-time cap.
#
#   scripts/run_capped.sh <maxGB> <max_seconds> <logfile> <command> [args...]
#
# The command's stdout+stderr go to <logfile>. The monitor polls the WHOLE
# process tree below the command and sums the resident sizes, so launching
# through a wrapper (/usr/bin/time, a shell, mpirun) cannot blind the cap --
# which is exactly the mistake that let a GENUS_Enumerate run eat all the
# memory on 2026-10-25: the monitor of that run watched the wrapper's RSS.
#
# max_seconds = 0 means no time limit: only the memory cap applies.
#
# Exit codes: the command's own code on normal termination;
#   137 when the memory cap was hit, 124 when the time cap was hit.
# A line MEMORY_CAPPED / TIME_CAPPED / peak RSS is appended to <logfile>.
#
# macOS has no enforceable ulimit -v, so polling (default every 10 s) is the
# mechanism; a very fast allocation spike can overshoot the cap by what it
# allocates within one poll interval.
set -u
if [ $# -lt 4 ]; then
    echo "usage: $0 <maxGB> <max_seconds> <logfile> <command> [args...]" >&2
    exit 2
fi
MAXGB=$1; MAXSEC=$2; LOG=$3; shift 3
MAXKB=$((MAXGB * 1024 * 1024))
POLL=${RUN_CAPPED_POLL:-10}

tree_pids() {
    ps -axo pid=,ppid= | awk -v root="$1" '
        {pid[NR]=$1; ppid[NR]=$2}
        END {
            keep[root]=1; added=1
            while (added) {
                added=0
                for (i=1; i<=NR; i++)
                    if (!(pid[i] in keep) && (ppid[i] in keep)) {
                        keep[pid[i]]=1; added=1
                    }
            }
            for (p in keep) print p
        }'
}

tree_rss_kb() {
    ps -axo pid=,ppid=,rss= | awk -v root="$1" '
        {pid[NR]=$1; ppid[NR]=$2; rss[NR]=$3}
        END {
            keep[root]=1; added=1
            while (added) {
                added=0
                for (i=1; i<=NR; i++)
                    if (!(pid[i] in keep) && (ppid[i] in keep)) {
                        keep[pid[i]]=1; added=1
                    }
            }
            s=0
            for (i=1; i<=NR; i++) if (pid[i] in keep) s+=rss[i]
            print s
        }'
}

kill_tree() {
    for p in $(tree_pids "$1"); do
        kill -9 "$p" 2>/dev/null
    done
}

"$@" > "$LOG" 2>&1 &
PID=$!
trap 'kill_tree $PID; exit 130' INT TERM
START=$(date +%s)
PEAK=0
while kill -0 "$PID" 2>/dev/null; do
    RSS=$(tree_rss_kb "$PID")
    [ -n "$RSS" ] && [ "$RSS" -gt "$PEAK" ] && PEAK=$RSS
    if [ -n "$RSS" ] && [ "$RSS" -gt "$MAXKB" ]; then
        echo "MEMORY_CAPPED: tree RSS ${RSS} KB > ${MAXKB} KB" >> "$LOG"
        kill_tree "$PID"
        exit 137
    fi
    if [ "$MAXSEC" -gt 0 ]; then
        ELAPSED=$(( $(date +%s) - START ))
        if [ "$ELAPSED" -gt "$MAXSEC" ]; then
            echo "TIME_CAPPED: ${ELAPSED} s > ${MAXSEC} s" >> "$LOG"
            kill_tree "$PID"
            exit 124
        fi
    fi
    sleep "$POLL"
done
wait "$PID"
RC=$?
echo "RUN_CAPPED: exit ${RC}, peak tree RSS ${PEAK} KB, elapsed $(( $(date +%s) - START )) s" >> "$LOG"
exit $RC
