# Devcontainer manual

A short, practical guide to using the `.devcontainer/` in this repository for LLM-assisted development on polyhedral_common.

## What this is

A pre-configured, containerized dev environment bundling:

- Node.js 22 and Python 3.
- The LLM CLIs `claude` (Anthropic Claude Code) and `codex` (OpenAI Codex CLI), each pre-wired to the same three MCP servers:
  - **serena** — semantic code navigation (via `uv`)
  - **sequentialthinking** — structured reasoning helper
  - **memory** — persistent knowledge graph

The environment is CLI-first: no IDE required. Everything is driven by the `@devcontainers/cli` tool and `docker`.

## Prerequisites

Install once, on your host:

```bash
# Docker Desktop (macOS/Windows) or docker-ce (Linux). Verify:
docker --version

# The devcontainer CLI:
npm install -g @devcontainers/cli
devcontainer --version
```

You do **not** need `claude` or `codex` on the host — they run inside the container.

## First run — build and start

From the repository root:

```bash
# One-time: create the named volume that will persist OAuth credentials,
# npm globals, uv cache, MCP memory graph, and bash history.
docker volume create poly_llm_home

# Build the image and start the container. Idempotent — safe to re-run.
devcontainer up --workspace-folder .
```

The first `devcontainer up` performs the full apt install + Node feature + npm globals dance. Expect several minutes. Subsequent runs reuse the built image and are nearly instant.

## Opening a shell

```bash
devcontainer exec --workspace-folder . bash
```

You land in `/workspace`, which is a bind mount of the repo root on your host. Edits inside the container are edits on your host filesystem — same inode, same file, same `git status`. Open as many shells as you like in parallel.

The container's welcome banner prints on every attach (via `postAttachCommand`), reminding you of the available CLIs and MCP servers.

## First-time OAuth login

Do this once per named volume. Credentials persist across `devcontainer up`.

Inside a container shell:

```bash
claude /login       # opens a URL; paste into your host browser, complete OAuth
codex login         # same pattern; will prompt for ChatGPT subscription
```

For the OAuth callback to reach the CLI, the callback port needs to be forwarded from host to container. That is handled two ways:

- **Linux**: change `runArgs` in `devcontainer.json` to `["--network=host"]` for the simplest OAuth flow. All ports the CLI picks are directly reachable.
- **macOS / Windows** (Docker Desktop): the default `runArgs` forwards ports `1455`, `8000`, and `54545` — the commonly seen values. If your CLI version picks a different port, add more `-p host:container` entries to `runArgs` and rebuild the container (`--remove-existing-container`).

Once authenticated, `~/.claude.json` and `~/.codex/auth.json` land in the named volume — every future container reuses them.

## Daily workflow

The bind mount means you can work from either side.

**Variant A — host-side commits:**

```bash
git switch -c my-feature
devcontainer exec --workspace-folder . bash
# inside: edit, build, test
exit
git add -p && git commit && git push
```

**Variant B — all-inside commits:**

```bash
devcontainer exec --workspace-folder . bash
# inside:
git switch -c my-feature
vim ... ; cmake --build build ; ./tests
git add -p && git commit && git push
```

Variant B needs git to know your identity *inside* the container. Either mount `~/.gitconfig` and `~/.ssh` through (add to `devcontainer.json` `mounts`) or configure git inside once and rely on the named volume to persist it.

Typical build / test cycle inside the container:

```bash
cd /workspace
cmake -B build -S . -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build build -j
cd CI_tests/<some-test-dir> && ./gap.sh Test<Name>.g
```

## Session control

There is no `devcontainer down` subcommand. Use `docker` directly:

```bash
# Find the container for the current workspace
docker ps --filter "label=devcontainer.local_folder=$(pwd)"

# Stop it (container survives, restart is instant)
docker stop $(docker ps -q --filter "label=devcontainer.local_folder=$(pwd)")

# Resume: exactly the same command as first time
devcontainer up --workspace-folder .

# Force a full rebuild after editing the Dockerfile or devcontainer.json
devcontainer up --workspace-folder . --remove-existing-container

# Remove the container (state in the named volume is preserved)
docker rm -f $(docker ps -aq --filter "label=devcontainer.local_folder=$(pwd)")
```

Full teardown, for when you want to wipe everything:

```bash
docker rm -f $(docker ps -aq --filter "label=devcontainer.local_folder=$(pwd)")
docker rmi $(docker images -q --filter "reference=vsc-polyhedral*")
docker volume rm poly_llm_home     # WARNING: wipes OAuth credentials
```

Optional convenience: drop a `Makefile` at the repo root with `dc-up`, `dc-shell`, `dc-stop`, `dc-rebuild`, `dc-nuke` targets — much friendlier than the raw docker incantations.

## Persistence model

Three kinds of state, three lifetimes.

| State | Where it lives | Persists across `--rm` | Persists across image rebuild |
|---|---|---|---|
| Source code | Bind-mounted from `${localWorkspaceFolder}` | yes (it's on your host) | yes |
| `/home/vscode` — OAuth, npm globals, uv cache, MCP memory, bash history | Named docker volume `poly_llm_home` | yes | yes |
| Everything else in `/` (tool binaries, `/tmp`, etc.) | Container's writable layer | no | no |

The named volume is copy-on-first-use: docker seeds it from the image's `/home/vscode` the first time you run `devcontainer up`, then never touches the image copy again. This has one consequence worth remembering:

**If you rebuild the image and expect updated MCP configs to appear, they will not.** The volume still holds the old versions. To pick up new defaults, delete the specific files inside the container — `rm ~/.claude.json ~/.codex/config.toml` — and re-attach; the `postStartCommand` re-seeds them from `/usr/local/share/poly_llm/` on the next start.

## MCP servers

All three are registered under both CLIs. Config templates live in `/usr/local/share/poly_llm/`; the entrypoint copies them to `~/.claude.json` and `~/.codex/config.toml` on every container start unless already present.

| MCP server | Command inside the container | Purpose |
|---|---|---|
| serena | `uvx --from git+https://github.com/oraios/serena serena-mcp-server` | Symbolic navigation (`find_symbol`, `get_symbols_overview`) — cuts token usage on large codebase reads. |
| sequentialthinking | `npx -y @modelcontextprotocol/server-sequential-thinking` | Structured multi-step reasoning helper. |
| memory | `npx -y @modelcontextprotocol/server-memory` | Persistent knowledge graph, useful across sessions. |

To disable one temporarily, remove or comment its entry from `~/.claude.json` (for Claude) or `~/.codex/config.toml` (for Codex). Since these files live in the named volume, the change survives container restarts.

## API-key mode (optional)

If you prefer pay-per-token billing over subscriptions:

```bash
docker rm -f $(docker ps -aq --filter "label=devcontainer.local_folder=$(pwd)")
export ANTHROPIC_API_KEY=sk-ant-...
export OPENAI_API_KEY=sk-...
devcontainer up --workspace-folder .
```

`devcontainer.json` already declares `remoteEnv` passthroughs for both variables — the CLIs will pick them up. You can mix modes freely (API key set → API mode; unset → subscription/OAuth).

## Troubleshooting

**OAuth completes in the browser but the CLI hangs.** The callback port isn't reaching the container. On macOS/Windows, check which port the CLI printed; if it's not in `1455 / 8000 / 54545`, add another `-p N:N` to `runArgs` and rebuild.

**`claude /login` opens `http://localhost:...` and the browser can't connect.** You may need `-p` forwarding as above, or (Linux) switch `runArgs` to `["--network=host"]`.

**Old MCP configs after image rebuild.** See the note under [Persistence model](#persistence-model). `rm ~/.claude.json ~/.codex/config.toml` inside the container, then exit and `devcontainer up` again.


**`git push` fails from inside the container.** SSH keys or gitconfig aren't shared. Either commit from the host, or add these `mounts` entries to `devcontainer.json` and rebuild:

```json
"source=${localEnv:HOME}/.gitconfig,target=/home/vscode/.gitconfig,type=bind,readonly",
"source=${localEnv:HOME}/.ssh,target=/home/vscode/.ssh,type=bind,readonly"
```

**Container filesystem fills with build artifacts.** They live in the container's writable layer, not the volume. `docker system df` shows usage; `docker builder prune` and rebuilding the container clear it.

**Serena/uvx is slow on first call.** The Dockerfile pre-warms the uv cache, but if the pre-warm failed (see build log for `WARN: Serena pre-warm failed`), the first `uvx` call downloads the Serena package from GitHub. One-time cost.

## Related files

- [`Dockerfile`](Dockerfile) — image definition (base + apt + uv + entrypoint script + MCP config templates)
- [`devcontainer.json`](devcontainer.json) — features, mounts, run args, lifecycle hooks
- [`../CLAUDE.md`](../CLAUDE.md) / [`../AGENTS.md`](../AGENTS.md) — project rules read by the LLM CLIs
