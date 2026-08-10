# Devcontainer image for polyhedral_common LLM-assisted development.
#
# Builds on Microsoft's reference C++ devcontainer base (ships with gcc,
# g++, clang, cmake, gdb, a non-root `vscode` user at UID 1000, sudo,
# locales, etc.), then layers:
#   * polyhedral_common build deps: boost, eigen, gmp,
#     openmpi, libtbb, lrslib
#   * nauty built from https://github.com/MathieuDutSik/nauty (tag 2.9.3)
#     into /opt/nauty, with $NAUTY_PATH / $NAUTY_INCLUDE / $NAUTY_LINK
#     exported so the per-directory Makefiles under src_*/ work directly
#   * debug & dev tools: valgrind, clangd, clang-tidy, clang-format,
#     ripgrep, fd-find, jq, tmux, etc.
#   * uv + Serena MCP (pre-warmed)
#   * MCP config templates (deployed by .devcontainer/devcontainer.json
#     lifecycle hooks; see /usr/local/share/poly_llm/)
#
# Node.js, the LLM CLIs (claude, codex), and the two npm-based MCPs
# (sequentialthinking, memory) are installed by the devcontainer
# Features mechanism + postCreateCommand — see devcontainer.json.

FROM mcr.microsoft.com/devcontainers/cpp:1-ubuntu-24.04

LABEL maintainer="mathieu.dutour@gmail.com"

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=Etc/UTC \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# ---------- Extra dev tooling on top of the MS C++ base ----------
RUN apt-get update -y && apt-get install -y --no-install-recommends \
        clangd clang-format clang-tidy lldb \
        valgrind strace bear \
        ripgrep fd-find jq tmux htop ncdu tree \
        less man-db manpages-dev bash-completion \
        emacs-nox nano \
        autoconf automake libtool \
        pipx \
    && rm -rf /var/lib/apt/lists/*

# fd is installed as `fdfind` on Debian/Ubuntu; alias to the conventional name
RUN ln -sf "$(command -v fdfind)" /usr/local/bin/fd

# ---------- polyhedral_common build dependencies ----------
RUN apt-get update -y && apt-get install -y --no-install-recommends \
        libboost-all-dev \
        libeigen3-dev \
        libgmp-dev \
        libtbb-dev \
        openmpi-bin libopenmpi-dev \
        lrslib \
    && rm -rf /var/lib/apt/lists/*

# ---------- nauty (built from the project's fork) ----------
# Ubuntu 24.04 no longer ships libnauty2-dev. Build it here so both CMake
# and the per-directory Makefiles (which reference $NAUTY_INCLUDE /
# $NAUTY_LINK) find headers and libraries out of the box.
# The top-level CMakeLists.txt still fetches its own copy via
# ExternalProject_Add — that lives inside the CMake build tree and is
# independent of this system-wide install.
RUN git clone --depth 1 --branch 2.9.3 \
        https://github.com/MathieuDutSik/nauty.git /tmp/nauty \
 && cd /tmp/nauty \
 && ./configure --prefix=/opt/nauty \
 && make -j"$(nproc)" \
 && make install \
 && rm -rf /tmp/nauty

# Env vars consumed by the per-directory Makefiles under src_*/.
ENV NAUTY_PATH=/opt/nauty
ENV NAUTY_INCLUDE="-I/opt/nauty/include"
ENV NAUTY_LINK="-L/opt/nauty/lib -lnauty -Wl,-rpath,/opt/nauty/lib"

# Standard env-var contract the per-directory Makefiles rely on. Setting
# these here keeps `make` in any src_* directory buildable without extra
# shell setup. `LIBBLISS_*` are intentionally left unset — bliss is not
# installed; the Makefile passes -I with empty argument (GCC tolerates it)
# and no bliss link fragment, which works for all targets that don't
# actually reference bliss headers.
ENV GMP_INCDIR=/usr/include
ENV GMP_CXX_LINK="-lgmp -lgmpxx"
ENV BOOST_INCDIR=/usr/include
ENV BOOST_LINK="-lboost_serialization"
ENV EIGEN_PATH=/usr/include/eigen3

# Sanity check: every header / pkg-config that CMakeLists.txt probes must
# resolve here. Failing fast at build time beats failing at devcontainer-up.
RUN pkg-config --exists eigen3 \
    && command -v cmake >/dev/null \
    && command -v mpic++ >/dev/null \
    && command -v g++ >/dev/null \
    && command -v clang++ >/dev/null \
    && echo "polyhedral_common build env OK"

# ---------- MCP wiring (templates baked into the image) ----------
# Stored under /usr/local/share/poly_llm so a host-mounted named volume
# on /home/vscode (the standard persistence pattern) can't shadow them.
# A postStartCommand in devcontainer.json copies them into the user's
# home only when absent — user overrides are preserved.
# Auth here comes from CLAUDE_CODE_OAUTH_TOKEN (env, set on the host), which
# satisfies API auth but does NOT mark first-run onboarding complete. Without
# these flags the interactive TUI still walks the theme / "Select login method"
# / trust-folder wizard on every fresh container. Pre-seeding them skips the
# wizard so the token is used directly.
RUN mkdir -p /usr/local/share/poly_llm && cat > /usr/local/share/poly_llm/claude.json <<'JSON'
{
  "hasCompletedOnboarding": true,
  "bypassPermissionsModeAccepted": true,
  "theme": "dark",
  "projects": {
    "/workspace": {
      "hasTrustDialogAccepted": true,
      "hasCompletedProjectOnboarding": true
    }
  },
  "mcpServers": {
    "serena": {
      "command": "uvx",
      "args": ["--from", "git+https://github.com/oraios/serena", "serena-mcp-server", "--context", "ide-assistant"]
    },
    "sequentialthinking": {
      "command": "npx",
      "args": ["-y", "@modelcontextprotocol/server-sequential-thinking"]
    },
    "memory": {
      "command": "npx",
      "args": ["-y", "@modelcontextprotocol/server-memory"]
    }
  }
}
JSON

RUN cat > /usr/local/share/poly_llm/codex.config.toml <<'TOML'
[mcp_servers.serena]
command = "uvx"
args = ["--from", "git+https://github.com/oraios/serena", "serena-mcp-server", "--context", "ide-assistant"]

[mcp_servers.sequentialthinking]
command = "npx"
args = ["-y", "@modelcontextprotocol/server-sequential-thinking"]

[mcp_servers.memory]
command = "npx"
args = ["-y", "@modelcontextprotocol/server-memory"]
TOML

# Restoration script: copies MCP config templates into HOME if missing.
# Driven by postStartCommand in devcontainer.json so it runs after any
# bind-mounts / named volumes are applied.
RUN cat > /usr/local/bin/poly_llm-restore.sh <<'SH'
#!/bin/bash
set -e
HOME_DIR="${HOME:-/home/vscode}"
mkdir -p "${HOME_DIR}/.codex"
if [ ! -f "${HOME_DIR}/.claude.json" ]; then
  cp /usr/local/share/poly_llm/claude.json "${HOME_DIR}/.claude.json"
fi
if [ ! -f "${HOME_DIR}/.codex/config.toml" ]; then
  cp /usr/local/share/poly_llm/codex.config.toml "${HOME_DIR}/.codex/config.toml"
fi
# Enforce always-on "YOLO" mode: this container is a disposable sandbox, so
# every claude session bypasses permission prompts. defaultMode drives the
# session; bypassPermissionsModeAccepted skips the one-time warning dialog.
# Merged (not cp'd) on every start because claude owns these files and
# rewrites them; the merge is idempotent and survives claude's own writes.
python3 - "${HOME_DIR}" <<'PY' || true
import json, os, sys
home = sys.argv[1]
def load(p):
    try:
        return json.load(open(p))
    except Exception:
        return {}
sp = os.path.join(home, ".claude", "settings.json")
os.makedirs(os.path.dirname(sp), exist_ok=True)
s = load(sp)
s.setdefault("permissions", {})["defaultMode"] = "bypassPermissions"
json.dump(s, open(sp, "w"), indent=2)
cp = os.path.join(home, ".claude.json")
c = load(cp)
c["bypassPermissionsModeAccepted"] = True
json.dump(c, open(cp, "w"), indent=2)
PY
SH
RUN chmod +x /usr/local/bin/poly_llm-restore.sh

# Welcome banner shown on every interactive shell (driven by postAttachCommand).
RUN cat > /usr/local/bin/poly_llm-welcome.sh <<'SH'
#!/bin/bash
cat <<'BANNER'
==========================================================
 polyhedral_common LLM dev container

 LLM CLIs:           claude     codex
 Auth:               host-side, no login here (CLAUDE_CODE_OAUTH_TOKEN
                     env + bind-mounted ~/.codex; see devcontainer.json)
 Permissions:        claude runs in bypass / YOLO mode (sandbox) — no
                     prompts; this is intended, do not run untrusted code
 MCP servers:        serena  sequentialthinking  memory
 Project root:       /workspace
==========================================================
BANNER
SH
RUN chmod +x /usr/local/bin/poly_llm-welcome.sh

# ---------- uv + Serena MCP (user-scope, baked) ----------
# uv installs to $HOME/.local/bin; ensure PATH covers it for non-login shells too.
USER vscode
ENV PATH="/home/vscode/.local/bin:${PATH}"

RUN curl -LsSf https://astral.sh/uv/install.sh | sh \
 && (uvx --from git+https://github.com/oraios/serena serena-mcp-server --help >/dev/null 2>&1 \
     || echo "WARN: Serena pre-warm failed; first MCP call will install on demand.")

# ---------- Shell quality of life ----------
RUN echo 'alias ll="ls -lah"' >> /home/vscode/.bashrc \
 && echo 'export EDITOR=vim' >> /home/vscode/.bashrc

# Workdir is set by devcontainer.json (workspaceFolder); leave a sensible default.
WORKDIR /workspace
