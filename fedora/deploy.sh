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
symlink "$FEDORA_DIR/.config/ghostty/config.ghostty"  "$HOME/.config/ghostty/config.ghostty"
symlink "$FEDORA_DIR/.config/ghostty/shaders/balatro.glsl" "$HOME/.config/ghostty/shaders/balatro.glsl"
symlink "$FEDORA_DIR/.config/ghostty/shaders/crt.glsl"     "$HOME/.config/ghostty/shaders/crt.glsl"
symlink "$FEDORA_DIR/.config/ghostty/shaders/test.glsl"    "$HOME/.config/ghostty/shaders/test.glsl"
if [ -f "$FEDORA_DIR/.p10k.zsh" ]; then
    symlink "$FEDORA_DIR/.p10k.zsh" "$HOME/.p10k.zsh"
fi

echo "Deploying Hyprland dotfiles"
symlink "$FEDORA_DIR/.config/hypr/hyprland.conf"      "$HOME/.config/hypr/hyprland.conf"

echo "Deploying caelestia shell config"
symlink "$FEDORA_DIR/.config/caelestia/shell.json"    "$HOME/.config/caelestia/shell.json"
symlink "$FEDORA_DIR/.local/bin/caelestia-gamemode"   "$HOME/.local/bin/caelestia-gamemode"
chmod +x "$FEDORA_DIR/.local/bin/caelestia-gamemode"

# Waybar and Mako configs are kept for reference but are no longer autostarted.
# caelestia-shell handles the bar, notifications, OSD, launcher, and lock screen.
# To restore the old setup: uncomment below and revert hyprland.conf autostart.
#symlink "$FEDORA_DIR/.config/waybar/config.jsonc"    "$HOME/.config/waybar/config.jsonc"
#symlink "$FEDORA_DIR/.config/waybar/style.css"       "$HOME/.config/waybar/style.css"
#symlink "$FEDORA_DIR/.config/mako/config"            "$HOME/.config/mako/config"

echo "Deploying shared dotfiles from $DOTFILES_DIR"
symlink "$DOTFILES_DIR/.config/i3/config"             "$HOME/.config/i3/config"
symlink "$DOTFILES_DIR/.config/i3/picom.conf"         "$HOME/.config/i3/picom.conf"
symlink "$DOTFILES_DIR/.config/i3/scripts/mouse.sh"   "$HOME/.config/i3/scripts/mouse.sh"

echo "Done."
