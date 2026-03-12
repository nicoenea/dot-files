#!/usr/bin/env bash
set -e

# Deploy fedora dotfiles by symlinking them into $HOME.
# Run as your regular user: ./deploy.sh

if [ "$EUID" -eq 0 ]; then
    echo "ERROR: Do not run this script as root or with sudo."
    exit 1
fi

FEDORA_DIR="$(cd "$(dirname "$0")" && pwd)"
DOTFILES_DIR="$(cd "$FEDORA_DIR/.." && pwd)"

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

echo "Deploying fedora dotfiles from $FEDORA_DIR"
symlink "$FEDORA_DIR/.zshrc"                          "$HOME/.zshrc"
symlink "$FEDORA_DIR/.zprofile"                       "$HOME/.zprofile"
symlink "$FEDORA_DIR/.profile"                        "$HOME/.profile"
symlink "$FEDORA_DIR/.config/fontconfig/fonts.conf"   "$HOME/.config/fontconfig/fonts.conf"
if [ -f "$FEDORA_DIR/.p10k.zsh" ]; then
    symlink "$FEDORA_DIR/.p10k.zsh" "$HOME/.p10k.zsh"
fi

echo "Deploying shared dotfiles from $DOTFILES_DIR"
symlink "$DOTFILES_DIR/.config/i3/config"             "$HOME/.config/i3/config"
symlink "$DOTFILES_DIR/.config/i3/picom.conf"         "$HOME/.config/i3/picom.conf"
symlink "$DOTFILES_DIR/.config/i3/scripts/mouse.sh"   "$HOME/.config/i3/scripts/mouse.sh"

echo "Done."
