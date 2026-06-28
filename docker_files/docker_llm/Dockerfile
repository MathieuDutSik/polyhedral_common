# This dockerfile provides a feature-rich development environment for
# working on polyhedral_common with LLM agents (Claude Code and Codex CLI).
#
# It bundles:
#   - C/C++ toolchain: gcc/g++/clang + clangd/clang-format/clang-tidy
#   - Debug tools: valgrind, gdb, lldb, strace
#   - Project build deps: cmake, boost, eigen3, gmp, glpk, cddlib,
#     openmpi, libnauty2, libtbb
#   - Node.js 22 LTS, Python 3 + uv
#   - LLM CLIs:
#       * @anthropic-ai/claude-code   (claude)
#       * @openai/codex                (codex)
#   - MCP servers, pre-wired into both Claude and Codex configs:
#       * serena              (semantic code navigation; via uv)
#       * sequentialthinking  (structured reasoning helper)
#       * memory              (persistent knowledge graph)
#   - Dev convenience: git, ripgrep, fd-find, jq, tmux, vim, emacs,
#     bash-completion, ssh-client
#
# It is used in the following way:
# ---Copy this "Dockerfile" to a directory, e.g. "docker_llm".
#    Do NOT rename the file.
# ---Compile the dockerfile with the command
#    docker build -t polyhedral_llm:dev .
# ---For building with no-cache option, use
#    docker build --no-cache -t polyhedral_llm:dev .
# ---Create one named volume to persist the dev user's HOME across runs
#    (OAuth credentials, MCP memory graph, uv cache, bash history). A
#    named volume is copy-on-first-use: docker seeds it from the image's
#    /home/dev on first run, then keeps it forever.
#      docker volume create poly_llm_home
# ---Run the container with your polyhedral_common checkout bind-mounted
#    at /workspace and the home volume attached:
#    docker run --rm -it \
#      --network host \                          # Linux: simplest OAuth callback
#      -v poly_llm_home:/home/dev \
#      -v "$(pwd)":/workspace \
#      polyhedral_llm:dev
#
#    macOS / Windows: Docker Desktop doesn't fully honor --network host.
#    Drop that flag and forward the OAuth callback port instead. The
#    exact port depends on the CLI version; common values are 1455 and
#    54545 for Claude Code and 1455 / 8000 for Codex. A safe shotgun is:
#      -p 1455:1455 -p 8000:8000 -p 54545:54545
#
# ---First-time login (do once; the named volume persists it):
#      claude /login         # opens a URL; paste it into your host browser
#      codex login           # likewise; will prompt for ChatGPT subscription
#    Subsequent `docker run` invocations reuse the stored credentials.
#
# ---No API keys are baked into the image. If you prefer pay-per-token
#    billing instead of subscriptions, add
#    `-e ANTHROPIC_API_KEY -e OPENAI_API_KEY` to the docker run line.
#
# Inside the container, run `claude` or `codex` from /workspace. The
# three MCPs (serena, sequentialthinking, memory) are auto-registered.
#
# Ask Mathieu Dutour Sikiric (mathieu.dutour@gmail.com) for questions.

FROM ubuntu:24.04

LABEL maintainer="mathieu.dutour@gmail.com"

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=Etc/UTC \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# ---------- Base system, toolchains, dev tools ----------
RUN apt-get update -y && apt-get install -y --no-install-recommends \
        ca-certificates curl wget gnupg lsb-release \
        sudo tzdata locales \
        bash-completion man-db manpages-dev less tree jq htop ncdu tmux \
        vim emacs-nox nano joe \
        git openssh-client \
        ripgrep fd-find \
        build-essential pkg-config make autoconf automake libtool \
        gcc g++ gfortran \
        clang clangd clang-format clang-tidy lldb \
        gdb valgrind strace \
        cmake bear \
        python3 python3-pip python3-venv pipx \
    && rm -rf /var/lib/apt/lists/*

# fd is installed as `fdfind` on Debian/Ubuntu; symlink to the conventional name
RUN ln -sf "$(command -v fdfind)" /usr/local/bin/fd

# ---------- polyhedral_common build dependencies ----------
RUN apt-get update -y && apt-get install -y --no-install-recommends \
        libboost-all-dev \
        libeigen3-dev \
        libgmp-dev \
        libglpk-dev \
        libcdd-dev \
        libtbb-dev \
        libnauty2-dev \
        openmpi-bin libopenmpi-dev \
        lrslib \
    && rm -rf /var/lib/apt/lists/*

# ---------- Node.js 22 LTS (for Claude Code, Codex, npm-based MCPs) ----------
RUN curl -fsSL https://deb.nodesource.com/setup_22.x | bash - \
    && apt-get install -y --no-install-recommends nodejs \
    && rm -rf /var/lib/apt/lists/*

# ---------- LLM CLIs and npm-based MCP servers (system-wide) ----------
# Installing the MCP servers globally means `npx -y <pkg>` resolves locally
# with no network round-trip on first invocation.
RUN npm install -g --silent \
        @anthropic-ai/claude-code \
        @openai/codex \
        @modelcontextprotocol/server-sequential-thinking \
        @modelcontextprotocol/server-memory

# ---------- Non-root dev user with sudo ----------
# Ubuntu 24.04 ships a default `ubuntu` user at UID 1000; remove it before
# creating `dev` so bind-mounts from a host UID-1000 user own files correctly.
RUN userdel --remove ubuntu 2>/dev/null || true \
    && groupdel ubuntu 2>/dev/null || true \
    && useradd --create-home --shell /bin/bash --uid 1000 dev \
    && echo 'dev ALL=(ALL) NOPASSWD:ALL' > /etc/sudoers.d/dev \
    && chmod 0440 /etc/sudoers.d/dev

USER dev
WORKDIR /home/dev

# ---------- uv + Serena MCP (user-scope, baked into the image) ----------
# uv installs to ~/.local/bin; ensure that's on PATH for every shell + for
# any MCP launcher (Claude / Codex spawn `uvx` and `serena-mcp-server`).
ENV PATH="/home/dev/.local/bin:${PATH}"

# Pre-warm the uv cache with Serena's git source so the first MCP call has
# no download delay. We use `uvx` rather than `uv tool install` to mirror
# Serena's documented invocation pattern and avoid coupling to whatever
# PyPI name the project picks internally.
RUN curl -LsSf https://astral.sh/uv/install.sh | sh \
 && (uvx --from git+https://github.com/oraios/serena serena-mcp-server --help >/dev/null 2>&1 \
     || echo "WARN: Serena pre-warm failed; first MCP call will install on demand.")

# ---------- MCP wiring (templates) ----------
# Stored under /usr/local/share/poly_llm so a host bind-mount onto
# ~/.codex or ~/.claude can't shadow them. The entrypoint copies these
# into the dev user's home only if absent — preserving any user override.
USER root
RUN mkdir -p /usr/local/share/poly_llm && cat > /usr/local/share/poly_llm/claude.json <<'JSON'
{
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

# Entrypoint: restores MCP configs after any bind-mounts have taken
# effect, then exec's the requested command (default: login bash).
RUN cat > /usr/local/bin/poly_llm-entrypoint.sh <<'SH'
#!/bin/bash
set -e
HOME_DIR="${HOME:-/home/dev}"
mkdir -p "${HOME_DIR}/.codex"
if [ ! -f "${HOME_DIR}/.claude.json" ]; then
  cp /usr/local/share/poly_llm/claude.json "${HOME_DIR}/.claude.json"
fi
if [ ! -f "${HOME_DIR}/.codex/config.toml" ]; then
  cp /usr/local/share/poly_llm/codex.config.toml "${HOME_DIR}/.codex/config.toml"
fi
exec "$@"
SH

RUN chmod +x /usr/local/bin/poly_llm-entrypoint.sh

USER dev

# ---------- Shell quality of life ----------
RUN echo '. /etc/bash_completion' >> /home/dev/.bashrc \
    && echo 'export PS1="\\[\\e[36m\\]docker_llm\\[\\e[0m\\]:\\w\\$ "' >> /home/dev/.bashrc \
    && echo 'alias ll="ls -lah"' >> /home/dev/.bashrc \
    && echo 'export EDITOR=vim' >> /home/dev/.bashrc

WORKDIR /workspace
ENTRYPOINT ["/usr/local/bin/poly_llm-entrypoint.sh"]
CMD ["/bin/bash", "-l"]
