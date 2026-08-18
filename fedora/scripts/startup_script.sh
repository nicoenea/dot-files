#!/usr/bin/env bash
set -e

# This script should be run as your regular user, NOT with sudo.
# It will prompt for sudo when needed.
if [ "$EUID" -eq 0 ]; then
    echo "ERROR: Do not run this script as root or with sudo."
    echo "Run it as your regular user: ./startup_script.sh"
    exit 1
fi

DOTFILES_DIR="$(cd "$(dirname "$0")/.." && pwd)"

# Visual Studio Code
sudo rpm --import https://packages.microsoft.com/keys/microsoft.asc
sudo sh -c 'echo -e "[code]\nname=Visual Studio Code\nbaseurl=https://packages.microsoft.com/yumrepos/vscode\nenabled=1\ngpgcheck=1\ngpgkey=https://packages.microsoft.com/keys/microsoft.asc" > /etc/yum.repos.d/vscode.repo'
sudo dnf check-update || true
sudo dnf install -y code

# pyenv dependencies + zsh plugin dependencies
sudo dnf install -y git gcc zlib-devel bzip2-devel readline-devel sqlite-devel openssl-devel python3-pygments

# pyenv (installed to user home)
if [ ! -d "$HOME/.pyenv" ]; then
    git clone https://github.com/pyenv/pyenv.git "$HOME/.pyenv"
else
    echo 'pyenv already installed, skipping.'
fi

# Zsh
echo 'Installing ZSH'
sudo dnf install -y zsh

if [ "$(getent passwd "$USER" | cut -d: -f7)" != "$(which zsh)" ]; then
    echo 'Making ZSH the default shell for '"$USER"'...'
    chsh -s "$(which zsh)"
else
    echo 'ZSH is already the default shell, skipping.'
fi

# Oh-My-ZSH (installed to user home)
if [ ! -d "$HOME/.oh-my-zsh" ]; then
    echo 'Installing Oh-My-Zsh'
    sh -c "$(curl -fsSL https://raw.github.com/ohmyzsh/ohmyzsh/master/tools/install.sh)" "" --unattended
else
    echo 'Oh-My-Zsh already installed, skipping.'
fi

# MesloLGS NF font (patched for Powerlevel10k)
FONT_DIR="$HOME/.local/share/fonts/MesloLGS"
if [ ! -d "$FONT_DIR" ]; then
    echo 'Installing MesloLGS Nerd Font...'
    mkdir -p "$FONT_DIR"
    BASE_URL="https://github.com/romkatv/powerlevel10k-media/raw/master"
    curl -fsSL "$BASE_URL/MesloLGS%20NF%20Regular.ttf"     -o "$FONT_DIR/MesloLGS NF Regular.ttf"
    curl -fsSL "$BASE_URL/MesloLGS%20NF%20Bold.ttf"        -o "$FONT_DIR/MesloLGS NF Bold.ttf"
    curl -fsSL "$BASE_URL/MesloLGS%20NF%20Italic.ttf"      -o "$FONT_DIR/MesloLGS NF Italic.ttf"
    curl -fsSL "$BASE_URL/MesloLGS%20NF%20Bold%20Italic.ttf" -o "$FONT_DIR/MesloLGS NF Bold Italic.ttf"
    fc-cache -fv "$FONT_DIR"
    echo 'MesloLGS NF installed. Set your terminal font to "MesloLGS NF" to enable all Powerlevel10k symbols.'
else
    echo 'MesloLGS NF already installed, skipping.'
fi

# Powerlevel10k
if [ ! -d "${ZSH_CUSTOM:-$HOME/.oh-my-zsh/custom}/themes/powerlevel10k" ]; then
    git clone --depth=1 https://github.com/romkatv/powerlevel10k.git "${ZSH_CUSTOM:-$HOME/.oh-my-zsh/custom}/themes/powerlevel10k"
else
    echo 'Powerlevel10k already installed, skipping.'
fi

# ohmyzsh plugins
if [ ! -d "${ZSH_CUSTOM:-$HOME/.oh-my-zsh/custom}/plugins/zsh-autosuggestions" ]; then
    git clone https://github.com/zsh-users/zsh-autosuggestions "${ZSH_CUSTOM:-$HOME/.oh-my-zsh/custom}/plugins/zsh-autosuggestions"
else
    echo 'zsh-autosuggestions already installed, skipping.'
fi

if [ ! -d "${ZSH_CUSTOM:-$HOME/.oh-my-zsh/custom}/plugins/zsh-syntax-highlighting" ]; then
    git clone https://github.com/zsh-users/zsh-syntax-highlighting.git "${ZSH_CUSTOM:-$HOME/.oh-my-zsh/custom}/plugins/zsh-syntax-highlighting"
else
    echo 'zsh-syntax-highlighting already installed, skipping.'
fi

# Deploy dotfiles
echo ''
echo 'Deploying dotfiles...'
"$DOTFILES_DIR/fedora/deploy.sh"

echo ''
echo 'Done! Restart your shell or run: exec zsh'
