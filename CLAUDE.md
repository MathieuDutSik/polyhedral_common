# CLAUDE.md

## `basic_common_cpp/` and `permutalib/` are git submodules — do not edit in place

Both `basic_common_cpp/` and `permutalib/` inside this repo are git submodules of their respective upstream projects. Any edit made directly to a file under either directory will be silently overwritten the next time the submodule is updated and will not propagate to the upstream.

To change a file owned by one of these submodules (e.g. `Temp_common.h`, `Basic_random.h`, `TypeConversion.h`, `Boost_bitset_kernel.h`, `hash_functions.h` in `basic_common_cpp/`; any header under `permutalib/`):

1. Edit the file in the upstream working copy of that submodule.
2. Commit and push the change there.
3. From this repo, run `./update.sh` (which calls `git submodule update --remote`) to advance the submodule pointer, then commit that pointer bump.

If unsure whether a path lives inside a submodule, `cat .gitmodules` lists them.

## Print statements must be gated by named `#ifdef` blocks

Free-standing `std::cerr << ...`, `os << ...`, etc. that produce output unconditionally must not be added to the codebase. Every diagnostic line must live inside a named preprocessor gate, in one of three categories with different intents. The suffix is the **context** — usually the source-file area (`SHVEC`, `LATTICE_STAB_EQUI_CAN`, `TSPACE_FUNCTIONS`, …) — never a generic name like `DEBUG`.

### `#ifdef DEBUG_<CONTEXT>`
For development-time tracing. May print verbosely. **Must not abort, throw, or otherwise affect the execution path** when the macro is undefined or defined. Example:
```cpp
#ifdef DEBUG_SHVEC
  os << "SHVEC: entering enumeration loop, norm_bound=" << B << "\n";
#endif
```

### `#ifdef SANITY_CHECK_<CONTEXT>`
For invariant checking. **Prints only when the invariant is violated**, then signals a programming error (`throw TerminalException{1}` or equivalent). When everything is fine, the block is silent. Example:
```cpp
#ifdef SANITY_CHECK_INVARIANT_VECTOR_FAMILY
  if (!IsSymmetricMatrix(M)) {
    std::cerr << "InvariantVectorFamily: input is not symmetric\n";
    throw TerminalException{1};
  }
#endif
```

### `#ifdef TIMINGS_<CONTEXT>`
For wall-clock measurement. Prints elapsed times around the operation being measured. Should rely on the in-tree timing helpers (`HumanTime`, `MicrosecondTime`, etc.) so output is consistent. Example:
```cpp
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  MicrosecondTime time;
#endif
  result = ExpensiveCall(...);
#ifdef TIMINGS_LATTICE_STAB_EQUI_CAN
  os << "LATTICE_STAB_EQUI_CAN: ExpensiveCall took " << time << "\n";
#endif
```

Existing macros in use can be discovered with `grep -rh '#ifdef \(DEBUG\|SANITY_CHECK\|TIMINGS\)_' src_*/ | sort -u`. Reuse an existing context name when adding to that area rather than inventing a new one.

## Builds must be warning-clean

A build of polyhedral_common code must produce no compiler warnings. The only warnings tolerated are those emitted from upstream third-party headers (Boost, Eigen, GMP, GLPK, cddlib, MPI, etc.) that we cannot patch. Warnings in code we own — including the submodules `basic_common_cpp/` and `permutalib/` (which should be patched upstream) — are bugs and must be fixed before merge.

When adding code, run a clean build and check the compiler output; do not assume "it worked last time" carries over.

## Reuse existing utilities before writing new ones

The repo has accumulated a substantial library of helpers across `basic_common_cpp/src_basic`, `basic_common_cpp/src_number`, `basic_common_cpp/src_matrix`, `src_basic`, `src_matrix`, `src_number`, `src_poly`, etc. Before writing a new helper, search for an existing one:

- For predicates and matrix utilities: `grep -rn 'IsSymmetric\|IsPositive\|RankMat\|ReadMatrix\|WriteMatrix' src_matrix/ basic_common_cpp/src_matrix/`
- For number-theory or type conversion: `basic_common_cpp/src_number/`
- For combinatorial / bitset / hash: `basic_common_cpp/src_comb/`, `basic_common_cpp/src_basic/`
- Use the Serena MCP `find_symbol` / `search_for_pattern` tools rather than re-reading whole files.

If a near-fit exists but isn't quite right, prefer extending or templatizing it over duplicating its body. Two copies of the same logic in different files is a bug — when the original is updated, the copy rots silently.

## CI is calendar-scheduled, one workflow per day

`.github/workflows/` follows a deliberately peculiar scheme that avoids saturating CI runners and keeps long-running test matrices spread out across the month. The rules:

1. **Each CI workflow runs at most once per day.** All triggers are `schedule: cron:` based; pushes do not trigger CI. (Workflows can also be fired by hand via `workflow_dispatch`.)
2. **The numeric index in the filename is the day of the month the workflow runs.** A file named `ci_NN_<slug>.yml` (or `ci_NNA_<slug>.yml` / `ci_NNB_<slug>.yml`) fires on day `NN`. For example `ci_14_schedule_full_sweep.yml` runs on the 14th: `cron: "0 0 14 * *"`.
3. **`A` / `B` suffix splits the day across alternating months.** When two workflows share a day, the `A` variant runs in odd months and `B` in even months:
   - `ci_02A_int_automorphy.yml`     → `cron: "0 0 2 1,3,5,7,9,11 *"` (Jan, Mar, May, Jul, Sep, Nov)
   - `ci_02B_rat_automorphy.yml`     → `cron: "0 0 2 2,4,6,8,10,12 *"` (Feb, Apr, Jun, Aug, Oct, Dec)

   So an `A` workflow fires every two months; same for `B`; together they cover every month, alternating.

### Practical consequences

- **Adding a new workflow.** Pick an unused day, name the file `ci_NN_<slug>.yml`, and set `cron: "0 0 NN * *"`. If the day is already taken by exactly one workflow that runs every month, convert that workflow to `ci_NNA_*` (odd months) and add yours as `ci_NNB_*` (even months) — update its cron accordingly. If both A and B slots are occupied, pick a different day.
- **You will not see your push trigger CI.** A change merged today will first be exercised by whichever workflow's day comes next. To verify immediately, trigger by hand with `gh workflow run ci_NN_<slug>.yml` (or the Actions UI's "Run workflow" button) — this works because every workflow also declares `workflow_dispatch:`.
- **Reading logs after a failure.** The user fetches CI failure logs via `gh run download` (or downloads the `L_<name>` zip from the Actions UI). When asked to diagnose a CI failure, expect the log location to be provided explicitly; do not assume it is under any particular path.
- **Bi-monthly cadence on `A` / `B` workflows.** An issue introduced just after, say, `ci_10A`'s February run will not be caught by `ci_10A` again until April. If a change is suspected to affect an `A`/`B` workflow, trigger it manually rather than waiting two months.
