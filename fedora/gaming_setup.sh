#!/usr/bin/env bash
set -e

# Fedora gaming setup for NVIDIA + Wayland
# Run as your regular user: ./gaming_setup.sh
# Re-running is safe — all steps are idempotent and will upgrade GE-Proton if a newer version is available.

if [ "$EUID" -eq 0 ]; then
    echo "ERROR: Do not run as root or with sudo."
    echo "Run as your regular user: ./gaming_setup.sh"
    exit 1
fi

echo ""
echo "========================================"
echo "  Fedora Gaming Setup (NVIDIA/Wayland)"
echo "========================================"
echo ""

# ---------------------------------------------------------------------------
# RPM Fusion
# ---------------------------------------------------------------------------
echo "==> Enabling RPM Fusion repos..."
if ! dnf repolist 2>/dev/null | grep -q "rpmfusion-free"; then
    sudo dnf install -y \
        "https://mirrors.rpmfusion.org/free/fedora/rpmfusion-free-release-$(rpm -E %fedora).noarch.rpm" \
        "https://mirrors.rpmfusion.org/nonfree/fedora/rpmfusion-nonfree-release-$(rpm -E %fedora).noarch.rpm"
else
    echo "  RPM Fusion already enabled, skipping."
fi

# ---------------------------------------------------------------------------
# Packages
# ---------------------------------------------------------------------------
echo ""
echo "==> Swapping ffmpeg-free for full ffmpeg (patented codecs)..."
if rpm -q ffmpeg-free &>/dev/null; then
    sudo dnf swap -y ffmpeg-free ffmpeg --allowerasing
else
    echo "  ffmpeg-free not installed (already swapped or using full ffmpeg), skipping."
fi

echo ""
echo "==> Installing gaming packages..."
sudo dnf install -y \
    steam \
    gamemode \
    gamemode.i686 \
    gamescope \
    libva-nvidia-driver \
    libva-utils \
    mangohud \
    mangohud.i686 \
    vulkan-tools

# ---------------------------------------------------------------------------
# Flatpak + ProtonUp-Qt
# ---------------------------------------------------------------------------
echo ""
echo "==> Setting up Flatpak and Flathub..."
if ! flatpak remotes | grep -q flathub; then
    flatpak remote-add --if-not-exists flathub https://dl.flathub.org/repo/flathub.flatpakrepo
else
    echo "  Flathub already configured, skipping."
fi

echo ""
echo "==> Installing ProtonUp-Qt (GUI manager for GE-Proton updates)..."
if ! flatpak list | grep -q "net.davidotek.pupgui2"; then
    flatpak install -y flathub net.davidotek.pupgui2
else
    echo "  ProtonUp-Qt already installed, skipping."
fi

# ---------------------------------------------------------------------------
# GE-Proton (latest release from GitHub)
# ---------------------------------------------------------------------------
echo ""
echo "==> Installing latest GE-Proton..."
COMPAT_DIR="$HOME/.steam/steam/compatibilitytools.d"
mkdir -p "$COMPAT_DIR"

RELEASE=$(curl -s "https://api.github.com/repos/GloriousEggroll/proton-ge-custom/releases/latest")
TAG=$(echo "$RELEASE" | grep '"tag_name"' | cut -d'"' -f4)
TARBALL_URL=$(echo "$RELEASE" | grep "browser_download_url.*\.tar\.gz" | cut -d'"' -f4)

if [ -z "$TAG" ]; then
    echo "  WARNING: Could not fetch latest GE-Proton release from GitHub. Skipping."
elif [ -d "$COMPAT_DIR/$TAG" ]; then
    echo "  GE-Proton $TAG already installed, skipping."
else
    echo "  Downloading GE-Proton $TAG..."
    curl -L --progress-bar "$TARBALL_URL" -o "/tmp/${TAG}.tar.gz"
    echo "  Extracting to $COMPAT_DIR..."
    tar -xzf "/tmp/${TAG}.tar.gz" -C "$COMPAT_DIR"
    rm "/tmp/${TAG}.tar.gz"
    echo "  GE-Proton $TAG installed."
fi

# ---------------------------------------------------------------------------
# NVIDIA Wayland modprobe
# ---------------------------------------------------------------------------
echo ""
echo "==> Writing NVIDIA Wayland modprobe config..."
MODPROBE_FILE="/etc/modprobe.d/nvidia-wayland.conf"
MODPROBE_CONTENT="options nvidia-drm modeset=1 fbdev=1"
NEEDS_REBOOT=0
if [ "$(cat "$MODPROBE_FILE" 2>/dev/null)" != "$MODPROBE_CONTENT" ]; then
    echo "$MODPROBE_CONTENT" | sudo tee "$MODPROBE_FILE" > /dev/null
    echo "  Written. Rebuilding initramfs (this may take a moment)..."
    sudo dracut --force
    NEEDS_REBOOT=1
else
    echo "  Already configured, skipping."
fi

# ---------------------------------------------------------------------------
# Wayland environment variables
# ---------------------------------------------------------------------------
echo ""
echo "==> Writing Wayland environment variables..."
ENV_DIR="$HOME/.config/environment.d"
ENV_FILE="$ENV_DIR/nvidia-wayland.conf"
mkdir -p "$ENV_DIR"
cat > "$ENV_FILE" <<'EOF'
LIBVA_DRIVER_NAME=nvidia
ELECTRON_OZONE_PLATFORM_HINT=auto
EOF
echo "  Written to $ENV_FILE"

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------
echo ""
echo "========================================"
if [ "$NEEDS_REBOOT" -eq 1 ]; then
    echo "  Setup complete! Reboot required to"
    echo "  apply kernel/modprobe changes."
else
    echo "  Setup complete! No reboot needed."
fi
echo "========================================"
echo ""
echo "Manual steps:"
echo ""
echo "  1. In Steam, right-click a game > Properties > Compatibility"
echo "     > Force compatibility tool > select GE-Proton"
echo ""
echo "  2. Add to each game's Steam launch options:"
echo "       gamemoderun PROTON_ENABLE_NVAPI=1 __GL_SYNC_TO_VBLANK=0 %command%"
echo ""
echo "  3. (Optional) Enable MangoHud overlay for any game:"
echo "       MANGOHUD=1 gamemoderun PROTON_ENABLE_NVAPI=1 __GL_SYNC_TO_VBLANK=0 %command%"
echo ""
echo "  4. To update GE-Proton later: re-run this script or use ProtonUp-Qt:"
echo "       flatpak run net.davidotek.pupgui2"
echo ""
echo "  5. Verify VA-API hardware decode is working:"
echo "       vainfo"
echo ""
echo "  6. (AoE2 DE) gamescope fixes fullscreen on multi-monitor Hyprland setups."
echo "     Set this as the launch option for Age of Empires II DE (app 813780):"
echo "       gamescope -w 3840 -h 2160 -f -e -- %command%"
echo "     Renders at 4K (downscaled to 2K monitor) for maximum map zoom-out."
echo "     Inside the game, use Fullscreen Windowed — gamescope will confine the mouse."
echo ""
