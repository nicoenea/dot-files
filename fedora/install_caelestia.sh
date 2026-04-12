#!/usr/bin/env bash
set -e

# Install caelestia-dots on Fedora.
# Ports the Arch-based caelestia setup (https://github.com/caelestia-dots/caelestia)
# to Fedora by building the two AUR-only components from source.
#
# Run as your regular user: ./install_caelestia.sh
# Optional flags:
#   --vscode    also symlink VS Code / VSCodium caelestia settings
#   --spotify   also set up Spicetify (you must install Spotify separately)

if [ "$EUID" -eq 0 ]; then
    echo "ERROR: Do not run as root. Run as your regular user."
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="$HOME/.local/src"
CAELESTIA_DIR="$HOME/.local/share/caelestia"

OPT_VSCODE=false
OPT_SPOTIFY=false
for arg in "$@"; do
    case "$arg" in
        --vscode)  OPT_VSCODE=true ;;
        --spotify) OPT_SPOTIFY=true ;;
    esac
done

info()    { echo -e "\e[1;34m==>\e[0m $*"; }
success() { echo -e "\e[1;32m  ✓\e[0m $*"; }
skip()    { echo -e "\e[1;33m  -\e[0m $* (already done, skipping)"; }

mkdir -p "$BUILD_DIR"

# ────────────────────────────────────────────────────────────
# 1. DNF PACKAGES
# ────────────────────────────────────────────────────────────
info "Installing DNF packages..."

# Core packages — all confirmed available in Fedora 43
sudo dnf install -y --skip-unavailable \
    `# Build tools` \
    cmake ninja-build clang clang-devel gcc-c++ git \
    `# Qt6 — needed by quickshell + caelestia-shell` \
    `# quickshell build deps (all discovered from its CMakeLists)` \
    cli11-devel \
    libdrm-devel mesa-libgbm-devel libglvnd-devel \
    wayland-devel wayland-protocols-devel \
    libxcb-devel xcb-util-devel xcb-util-wm-devel \
    jemalloc-devel \
    pipewire-devel \
    polkit-devel glib2-devel pam-devel \
    vulkan-headers vulkan-loader-devel \
    `# cpptrace not in Fedora repos — crash handling disabled via cmake flag` \
    `# Qt6` \
    qt6-qtbase-devel qt6-qtbase-private-devel \
    qt6-qtdeclarative-devel qt6-qtdeclarative-private-devel \
    qt6-qtwayland-devel \
    qt6-qtsvg-devel qt6-qtmultimedia-devel qt6-qttranslations \
    qt6-qtshadertools-devel \
    `# Wayland / Hyprland ecosystem` \
    hyprland xdg-desktop-portal-hyprland xdg-desktop-portal-gtk \
    hyprpicker wireplumber \
    `# Hardware control (power-profiles-daemon conflicts with tuned-ppd on Fedora 43;` \
    `# tuned-ppd is already installed and serves the same role)` \
    ddcutil brightnessctl lm_sensors \
    `# Audio / visualisation` \
    aubio-devel pipewire-libs libqalculate-devel \
    `# libcava build deps (cava from Fedora repos has no pkg-config file)` \
    fftw-devel iniparser-devel alsa-lib-devel libtool autoconf automake \
    `# Fish shell (required by caelestia-cli scripts — your default shell stays zsh)` \
    fish \
    `# Utilities` \
    fastfetch btop jq inotify-tools trash-cli \
    fuzzel \
    wl-clipboard cliphist \
    grim slurp swappy \
    `# Python (caelestia-cli needs 3.13)` \
    python3 python3-pip python3-pillow python3-devel \
    `# Fonts infra` \
    fontconfig \
    `# GTK / icon themes` \
    papirus-icon-theme adw-gtk3-theme \
    `# Session / auth (polkit-gnome was dropped; the exec in hyprland.conf` \
    `# already falls back to /usr/libexec/polkit-gnome-authentication-agent-1)` \
    gnome-keyring uwsm \
    `# Misc` \
    playerctl dconf

success "DNF packages installed"

# ────────────────────────────────────────────────────────────
# 1b. PACKAGES NOT IN FEDORA REPOS
# ────────────────────────────────────────────────────────────

# starship — not in Fedora repos, install via official script
if ! command -v starship &>/dev/null; then
    info "Installing starship prompt..."
    curl -sS https://starship.rs/install.sh | sh -s -- --yes
    success "starship installed"
else
    skip "starship"
fi

# eza — not in Fedora 43 repos, install via cargo
if ! command -v eza &>/dev/null; then
    info "Installing eza (modern ls)..."
    if ! command -v cargo &>/dev/null; then
        sudo dnf install -y cargo
    fi
    cargo install eza
    success "eza installed (~/.cargo/bin/eza)"
else
    skip "eza"
fi

# ────────────────────────────────────────────────────────────
# 2. FONTS
# ────────────────────────────────────────────────────────────
info "Installing fonts..."

FONT_DIR="$HOME/.local/share/fonts"
mkdir -p "$FONT_DIR"

# JetBrains Mono Nerd Font
if [ ! -f "$FONT_DIR/JetBrainsMonoNerdFont-Regular.ttf" ]; then
    info "  Downloading JetBrains Mono Nerd Font..."
    JB_VER="$(curl -s https://api.github.com/repos/ryanoasis/nerd-fonts/releases/latest | jq -r '.tag_name')"
    curl -fsSL "https://github.com/ryanoasis/nerd-fonts/releases/download/${JB_VER}/JetBrainsMono.tar.xz" \
        | tar -xJ -C "$FONT_DIR" --wildcards '*.ttf'
    success "  JetBrains Mono Nerd Font installed"
else
    skip "  JetBrains Mono Nerd Font"
fi

# Material Symbols Variable font (required by caelestia-shell UI)
if [ ! -f "$FONT_DIR/MaterialSymbolsRounded[FILL,GRAD,opsz,wght].ttf" ]; then
    info "  Downloading Material Symbols..."
    MS_BASE="https://github.com/google/material-design-icons/raw/master/variablefont"
    curl -fsSL "${MS_BASE}/MaterialSymbolsRounded%5BFILL%2CGRAD%2Copsz%2Cwght%5D.ttf" \
        -o "$FONT_DIR/MaterialSymbolsRounded[FILL,GRAD,opsz,wght].ttf"
    success "  Material Symbols installed"
else
    skip "  Material Symbols"
fi

# Rubik variable font
if [ ! -f "$FONT_DIR/Rubik-VariableFont_wght.ttf" ]; then
    info "  Downloading Rubik..."
    curl -fsSL "https://github.com/googlefonts/rubik/raw/main/fonts/variable/Rubik%5Bwght%5D.ttf" \
        -o "$FONT_DIR/Rubik-VariableFont_wght.ttf"
    success "  Rubik installed"
else
    skip "  Rubik"
fi

fc-cache -fv "$FONT_DIR" > /dev/null
success "Font cache updated"

# ────────────────────────────────────────────────────────────
# 3. app2unit  (shell script wrapper, no AUR needed)
# ────────────────────────────────────────────────────────────
if ! command -v app2unit &>/dev/null; then
    info "Installing app2unit..."
    APP2UNIT_DIR="$BUILD_DIR/app2unit"
    git clone --depth=1 https://github.com/Vladimir-csp/app2unit.git "$APP2UNIT_DIR" 2>/dev/null || \
        git -C "$APP2UNIT_DIR" pull
    sudo install -Dm755 "$APP2UNIT_DIR/app2unit" /usr/local/bin/app2unit
    success "app2unit installed"
else
    skip "app2unit"
fi

# ────────────────────────────────────────────────────────────
# 4. dart-sass  (binary release — caelestia-cli uses it for theming)
# ────────────────────────────────────────────────────────────
if ! command -v sass &>/dev/null; then
    info "Installing dart-sass..."
    SASS_VER="$(curl -s https://api.github.com/repos/sass/dart-sass/releases/latest | jq -r '.tag_name')"
    curl -fsSL "https://github.com/sass/dart-sass/releases/download/${SASS_VER}/dart-sass-${SASS_VER}-linux-x64.tar.gz" \
        | tar -xz -C /tmp
    sudo cp /tmp/dart-sass/sass /usr/local/bin/sass
    sudo cp /tmp/dart-sass/src/dart /usr/local/bin/dart-sass-dart
    rm -rf /tmp/dart-sass
    success "dart-sass installed"
else
    skip "dart-sass"
fi

# ────────────────────────────────────────────────────────────
# 5. libcava  (Fedora's cava package has no pkg-config file; must build as lib)
# ────────────────────────────────────────────────────────────
if ! pkg-config --exists cava 2>/dev/null; then
    info "Building libcava from source..."
    CAVA_DIR="$BUILD_DIR/cava"
    git clone --depth=1 https://github.com/karlstav/cava.git "$CAVA_DIR" 2>/dev/null || \
        git -C "$CAVA_DIR" pull

    # Compile cavacore.c as a shared library (what Arch's libcava package does)
    gcc -shared -fPIC -o "$CAVA_DIR/libcava.so" \
        "$CAVA_DIR/cavacore.c" \
        $(pkg-config --cflags --libs fftw3) -lm

    sudo install -Dm755 "$CAVA_DIR/libcava.so"    /usr/local/lib/libcava.so
    sudo install -Dm644 "$CAVA_DIR/cavacore.h"    /usr/local/include/cava/cavacore.h

    # Create the pkg-config file caelestia-shell looks for
    sudo mkdir -p /usr/lib64/pkgconfig
    sudo tee /usr/lib64/pkgconfig/cava.pc > /dev/null <<EOF
prefix=/usr/local
exec_prefix=\${prefix}
libdir=\${exec_prefix}/lib
includedir=\${prefix}/include

Name: cava
Description: Console-based Audio Visualizer for Alsa (library)
Version: 0.10.0
Libs: -L\${libdir} -lcava
Cflags: -I\${includedir}
EOF
    # Register /usr/local/lib so the dynamic linker can find libcava.so
    echo '/usr/local/lib' | sudo tee /etc/ld.so.conf.d/usr-local.conf > /dev/null
    sudo ldconfig
    success "libcava installed"
else
    skip "libcava"
fi

# ────────────────────────────────────────────────────────────
# 6. quickshell  (the QML shell framework — AUR only, must build from source)
# ────────────────────────────────────────────────────────────
if ! command -v qs &>/dev/null; then
    info "Building quickshell from source (this will take a while)..."
    QS_DIR="$BUILD_DIR/quickshell"
    git clone --depth=1 --recurse-submodules \
        https://github.com/quickshell-mirror/quickshell.git "$QS_DIR" 2>/dev/null || \
        git -C "$QS_DIR" pull
    rm -rf "$QS_DIR/build"
    cmake -B "$QS_DIR/build" "$QS_DIR" \
        -GNinja \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DINSTALL_QMLDIR=/usr/local/lib64/qt6/qml \
        -DCRASH_HANDLER=OFF
    ninja -C "$QS_DIR/build"
    sudo ninja -C "$QS_DIR/build" install
    success "quickshell installed"
else
    skip "quickshell (qs already in PATH)"
fi

# ────────────────────────────────────────────────────────────
# 6. caelestia-cli  (Python package — pip install from source)
# ────────────────────────────────────────────────────────────
if ! command -v caelestia &>/dev/null; then
    info "Installing caelestia-cli..."
    CLI_DIR="$BUILD_DIR/caelestia-cli"
    git clone --depth=1 https://github.com/caelestia-dots/cli.git "$CLI_DIR" 2>/dev/null || \
        git -C "$CLI_DIR" pull
    pip install --user "$CLI_DIR"
    # Ensure ~/.local/bin is in PATH for this session
    export PATH="$HOME/.local/bin:$PATH"
    success "caelestia-cli installed"
else
    skip "caelestia-cli"
fi

# ────────────────────────────────────────────────────────────
# 7. caelestia-shell  (C++/QML — must build from source)
# ────────────────────────────────────────────────────────────
if ! command -v caelestia-shell &>/dev/null && ! [ -f /usr/local/lib/quickshell/caelestia/shell.qml ]; then
    info "Building caelestia-shell from source (this will take a while)..."
    SHELL_DIR="$BUILD_DIR/caelestia-shell"
    git clone --depth=1 https://github.com/caelestia-dots/shell.git "$SHELL_DIR" 2>/dev/null || \
        git -C "$SHELL_DIR" pull
    cmake -B "$SHELL_DIR/build" "$SHELL_DIR" \
        -GNinja \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DVERSION=0.1.0
    ninja -C "$SHELL_DIR/build"
    sudo ninja -C "$SHELL_DIR/build" install
    # quickshell looks in ~/.config/quickshell/ — symlink the installed config there
    mkdir -p "$HOME/.config/quickshell"
    ln -sfn /usr/local/etc/xdg/quickshell/caelestia "$HOME/.config/quickshell/caelestia"
    # Fix doubled install path: caelestia cmake stacks /usr/lib/qt6/qml onto our prefix
    sudo ln -sfn /usr/local/usr/lib/qt6/qml/Caelestia /usr/local/lib64/qt6/qml/Caelestia
    success "caelestia-shell installed"
else
    skip "caelestia-shell"
fi

# Ensure the quickshell symlink exists even if caelestia-shell was already built
mkdir -p "$HOME/.config/quickshell"
ln -sfn /usr/local/etc/xdg/quickshell/caelestia "$HOME/.config/quickshell/caelestia"

# ────────────────────────────────────────────────────────────
# 8. caelestia dotfiles  (configs for btop, fastfetch, fish, foot, etc.)
# ────────────────────────────────────────────────────────────
info "Cloning caelestia dotfiles..."
if [ ! -d "$CAELESTIA_DIR/.git" ]; then
    git clone --depth=1 https://github.com/caelestia-dots/caelestia.git "$CAELESTIA_DIR"
    success "Cloned caelestia dotfiles to $CAELESTIA_DIR"
else
    git -C "$CAELESTIA_DIR" pull
    success "Updated caelestia dotfiles"
fi

# Symlink caelestia configs — but NOT hyprland (we manage that separately)
info "Symlinking caelestia configs..."

symlink_caelestia() {
    local name="$1"
    local src="$CAELESTIA_DIR/$name"
    local dst="$HOME/.config/$name"
    [ -d "$src" ] || return 0
    if [ -e "$dst" ] && [ ! -L "$dst" ]; then
        mv "$dst" "${dst}.bak"
        echo "  Backed up existing $dst -> ${dst}.bak"
    fi
    mkdir -p "$(dirname "$dst")"
    ln -sfn "$src" "$dst"
    echo "  Linked ~/.config/$name -> $src"
}

symlink_caelestia btop
symlink_caelestia fastfetch
symlink_caelestia fish
# foot is caelestia's default terminal, but we use Ghostty — skip foot config.
symlink_caelestia starship.toml
symlink_caelestia micro
symlink_caelestia uwsm

# caelestia shell config dir (not the same as the shell repo)
mkdir -p "$HOME/.config/caelestia"
if [ ! -f "$HOME/.config/caelestia/shell.json" ]; then
    ln -sfn "$SCRIPT_DIR/.config/caelestia/shell.json" "$HOME/.config/caelestia/shell.json"
    echo "  Linked ~/.config/caelestia/shell.json"
fi

success "caelestia configs symlinked"

# ────────────────────────────────────────────────────────────
# 9. OPTIONAL: VS Code settings
# ────────────────────────────────────────────────────────────
if $OPT_VSCODE; then
    info "Symlinking VS Code caelestia settings..."
    VSCODE_CFG="$HOME/.config/Code/User"
    mkdir -p "$VSCODE_CFG"
    [ -f "$VSCODE_CFG/settings.json" ] && mv "$VSCODE_CFG/settings.json" "$VSCODE_CFG/settings.json.bak"
    ln -sfn "$CAELESTIA_DIR/vscode/settings.json" "$VSCODE_CFG/settings.json"
    ln -sfn "$CAELESTIA_DIR/vscode/keybindings.json" "$VSCODE_CFG/keybindings.json"
    success "VS Code settings linked"
fi

# ────────────────────────────────────────────────────────────
# 10. OPTIONAL: Spicetify
# ────────────────────────────────────────────────────────────
if $OPT_SPOTIFY; then
    if command -v spicetify &>/dev/null && command -v spotify &>/dev/null; then
        info "Setting up Spicetify caelestia theme..."
        mkdir -p "$HOME/.config/spicetify/Themes"
        ln -sfn "$CAELESTIA_DIR/spicetify/Themes/caelestia" \
            "$HOME/.config/spicetify/Themes/caelestia"
        spicetify config current_theme caelestia
        spicetify apply
        success "Spicetify theme applied"
    else
        echo "  Skipping Spicetify: spotify or spicetify not installed."
    fi
fi

# ────────────────────────────────────────────────────────────
echo ""
echo -e "\e[1;32m✓ caelestia installation complete!\e[0m"
echo ""
echo "Next steps:"
echo "  1. Add ~/.local/bin to your PATH (if not already):"
echo "       echo 'export PATH=\"\$HOME/.local/bin:\$PATH\"' >> ~/.zshrc"
echo "  2. Restart Hyprland — caelestia shell will autostart."
echo "  3. Customise ~/.config/caelestia/shell.json to your liking."
echo "  4. Set a wallpaper: place images in ~/Pictures/Wallpapers"
echo "  5. Set a profile picture: cp <photo> ~/.face"
echo ""
echo "  Tip: run 'caelestia scheme' after a wallpaper change to regenerate"
echo "       the Material You color theme across all apps."
