#!/usr/bin/env bash
set -e

# Deploy fedora dotfiles by symlinking them into $HOME.
# Run as your regular user: ./deploy.sh

if [ "$EUID" -eq 0 ]; then
    echo "ERROR: Do not run this script as root or with sudo."
    exit 1
fi

FEDORA_DIR="$(cd "$(dirname "$0")" && pwd)"

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

# Flatpak sandboxes can't follow a symlink into ~/repos, so these get copied.
copy() {
    install -Dm644 "$1" "$2"
    echo "  Copied $1 -> $2"
}

echo "Deploying fedora dotfiles from $FEDORA_DIR"
symlink "$FEDORA_DIR/.zshrc"                          "$HOME/.zshrc"
symlink "$FEDORA_DIR/.zprofile"                       "$HOME/.zprofile"
symlink "$FEDORA_DIR/.profile"                        "$HOME/.profile"
symlink "$FEDORA_DIR/.config/fontconfig/fonts.conf"   "$HOME/.config/fontconfig/fonts.conf"
symlink "$FEDORA_DIR/.config/ghostty/config.ghostty"  "$HOME/.config/ghostty/config.ghostty"
symlink "$FEDORA_DIR/.config/ghostty/shaders"          "$HOME/.config/ghostty/shaders"
if [ -f "$FEDORA_DIR/.p10k.zsh" ]; then
    symlink "$FEDORA_DIR/.p10k.zsh" "$HOME/.p10k.zsh"
fi

symlink "$FEDORA_DIR/.config/environment.d/gpu-wayland.conf" "$HOME/.config/environment.d/gpu-wayland.conf"
copy "$FEDORA_DIR/.var/app/com.brave.Browser/config/brave-flags.conf" "$HOME/.var/app/com.brave.Browser/config/brave-flags.conf"
symlink "$FEDORA_DIR/.config/fuzzel/fuzzel.ini"       "$HOME/.config/fuzzel/fuzzel.ini"

echo "Deploying GTK settings (icon theme, cursor, fonts)"
symlink "$FEDORA_DIR/.config/gtk-3.0/settings.ini"    "$HOME/.config/gtk-3.0/settings.ini"
symlink "$FEDORA_DIR/.config/gtk-4.0/settings.ini"    "$HOME/.config/gtk-4.0/settings.ini"

echo "Deploying Hyprland dotfiles"
symlink "$FEDORA_DIR/.config/hypr/hyprland.lua"       "$HOME/.config/hypr/hyprland.lua"
symlink "$FEDORA_DIR/.config/hypr/hyprland"           "$HOME/.config/hypr/hyprland"
# Hyprland >= 0.55 uses the Lua config; drop the old hyprlang symlink so the two never conflict.
if [ -L "$HOME/.config/hypr/hyprland.conf" ]; then
    rm "$HOME/.config/hypr/hyprland.conf"
    echo "  Removed stale $HOME/.config/hypr/hyprland.conf symlink"
fi
symlink "$FEDORA_DIR/.config/hypr/xdph.conf"          "$HOME/.config/hypr/xdph.conf"

echo "Deploying caelestia shell config"
symlink "$FEDORA_DIR/.config/caelestia/shell.json"    "$HOME/.config/caelestia/shell.json"
symlink "$FEDORA_DIR/.local/bin/caelestia-gamemode"   "$HOME/.local/bin/caelestia-gamemode"
chmod +x "$FEDORA_DIR/.local/bin/caelestia-gamemode"
symlink "$FEDORA_DIR/.local/bin/update-caelestia"     "$HOME/.local/bin/update-caelestia"
chmod +x "$FEDORA_DIR/.local/bin/update-caelestia"
symlink "$FEDORA_DIR/.local/bin/skydimo-init"          "$HOME/.local/bin/skydimo-init"
chmod +x "$FEDORA_DIR/.local/bin/skydimo-init"
symlink "$FEDORA_DIR/.local/bin/hyperhdr-picker"       "$HOME/.local/bin/hyperhdr-picker"
chmod +x "$FEDORA_DIR/.local/bin/hyperhdr-picker"

# Rapid-fire mouse scripts (bound to Super+R / Super+M in hypr/hyprland/keybinds.lua; needs python3-evdev)
symlink "$FEDORA_DIR/.local/bin/rapidfire-daemon"      "$HOME/.local/bin/rapidfire-daemon"
chmod +x "$FEDORA_DIR/.local/bin/rapidfire-daemon"
symlink "$FEDORA_DIR/.local/bin/rapidfire-toggle"      "$HOME/.local/bin/rapidfire-toggle"
chmod +x "$FEDORA_DIR/.local/bin/rapidfire-toggle"

# Login/daily update check (flatpak auto-update, dnf notify, dead-copr detector)
echo "Deploying update-check timer"
symlink "$FEDORA_DIR/.local/bin/update-check"          "$HOME/.local/bin/update-check"
chmod +x "$FEDORA_DIR/.local/bin/update-check"
symlink "$FEDORA_DIR/.config/systemd/user/update-check.service" "$HOME/.config/systemd/user/update-check.service"
symlink "$FEDORA_DIR/.config/systemd/user/update-check.timer"   "$HOME/.config/systemd/user/update-check.timer"
systemctl --user daemon-reload
systemctl --user enable --now update-check.timer

echo "Done."
echo ""
echo "NOTE: To revert GRUB to Fedora defaults (removes any custom theme), run:"
echo "  sudo ./scripts/grub_setup.sh"
