#!/usr/bin/env bash
set -e

# Deploy wsl2-ubuntu dotfiles: installs zsh, oh-my-zsh, plugins, pyenv, then symlinks configs.
# Run as your regular user: ./deploy.sh

if [ "$EUID" -eq 0 ]; then
    echo "ERROR: Do not run this script as root or with sudo."
    exit 1
fi

WSL_DIR="$(cd "$(dirname "$0")" && pwd)"

symlink() {
    local src="$1"
    local dst="$2"
    if [ -e "$dst" ] && [ ! -L "$dst" ]; then
        echo "  Backing up existing $dst -> ${dst}.bak"
        mv "$dst" "${dst}.bak"
    fi
    mkdir -p "$(dirname "$dst")"
    ln -sfn "$src" "$dst"
    echo "  Linked $dst -> $src"
}

# ── 1. System packages ────────────────────────────────────────────────────────
echo "Installing system packages..."
sudo apt-get update -qq
sudo apt-get install -y \
    zsh \
    git \
    curl \
    wget \
    build-essential \
    libssl-dev \
    zlib1g-dev \
    libbz2-dev \
    libreadline-dev \
    libsqlite3-dev \
    libncursesw5-dev \
    xz-utils \
    tk-dev \
    libxml2-dev \
    libxmlsec1-dev \
    libffi-dev \
    liblzma-dev \
    fortune-mod \
    fortunes \
    python3-pygments

# ── 2. oh-my-zsh ──────────────────────────────────────────────────────────────
if [ ! -d "$HOME/.oh-my-zsh" ]; then
    echo "Installing oh-my-zsh..."
    RUNZSH=no CHSH=no sh -c "$(curl -fsSL https://raw.githubusercontent.com/ohmyzsh/ohmyzsh/master/tools/install.sh)"
else
    echo "oh-my-zsh already installed, skipping."
fi

ZSH_CUSTOM="${ZSH_CUSTOM:-$HOME/.oh-my-zsh/custom}"

# ── 3. Powerlevel10k ──────────────────────────────────────────────────────────
if [ ! -d "$ZSH_CUSTOM/themes/powerlevel10k" ]; then
    echo "Installing Powerlevel10k..."
    git clone --depth=1 https://github.com/romkatv/powerlevel10k.git "$ZSH_CUSTOM/themes/powerlevel10k"
else
    echo "Powerlevel10k already installed, skipping."
fi

# ── 4. zsh-syntax-highlighting ────────────────────────────────────────────────
if [ ! -d "$ZSH_CUSTOM/plugins/zsh-syntax-highlighting" ]; then
    echo "Installing zsh-syntax-highlighting..."
    git clone --depth=1 https://github.com/zsh-users/zsh-syntax-highlighting.git "$ZSH_CUSTOM/plugins/zsh-syntax-highlighting"
else
    echo "zsh-syntax-highlighting already installed, skipping."
fi

# ── 5. zsh-autosuggestions ────────────────────────────────────────────────────
if [ ! -d "$ZSH_CUSTOM/plugins/zsh-autosuggestions" ]; then
    echo "Installing zsh-autosuggestions..."
    git clone --depth=1 https://github.com/zsh-users/zsh-autosuggestions.git "$ZSH_CUSTOM/plugins/zsh-autosuggestions"
else
    echo "zsh-autosuggestions already installed, skipping."
fi

# ── 6. pyenv ──────────────────────────────────────────────────────────────────
if [ ! -d "$HOME/.pyenv" ]; then
    echo "Installing pyenv..."
    curl -fsSL https://pyenv.run | bash
else
    echo "pyenv already installed, skipping."
fi

# ── 7. Symlink dotfiles ───────────────────────────────────────────────────────
echo "Deploying wsl2-ubuntu dotfiles from $WSL_DIR"
symlink "$WSL_DIR/.zshrc"    "$HOME/.zshrc"
symlink "$WSL_DIR/.zprofile" "$HOME/.zprofile"
symlink "$WSL_DIR/.profile"  "$HOME/.profile"
if [ -f "$WSL_DIR/.p10k.zsh" ]; then
    symlink "$WSL_DIR/.p10k.zsh" "$HOME/.p10k.zsh"
fi

# ── 8. Set zsh as default shell ───────────────────────────────────────────────
ZSH_PATH="$(which zsh)"
if [ "$SHELL" != "$ZSH_PATH" ]; then
    echo "Setting zsh as default shell..."
    chsh -s "$ZSH_PATH"
    echo "  Done. Log out and back in for the change to take effect."
else
    echo "zsh is already the default shell."
fi

echo ""
echo "Done! Start a new zsh session or run: exec zsh"
echo "NOTE: Make sure your terminal uses a Nerd Font for the Powerlevel10k icons."
