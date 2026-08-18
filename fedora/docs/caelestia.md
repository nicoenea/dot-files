# caelestia on Fedora 43

Complete reference for the caelestia-dots Hyprland rice port to Fedora 43.
Upstream: https://github.com/caelestia-dots/caelestia

---

## What is caelestia?

A Material Design 3 / Material You Hyprland desktop rice consisting of:

| Component | Purpose | Replaces |
|---|---|---|
| **caelestia-shell** | C++/QML desktop shell (bar, notifications, OSD, launcher, lock screen) | Waybar + Mako |
| **caelestia-cli** | Python control tool — theming, wallpaper, color scheme generation | – |
| **quickshell** | QML shell framework that caelestia-shell runs on top of | – |
| **caelestia dotfiles** | Configs for btop, fastfetch, fish, starship, micro, uwsm | – |

**Shell stays Zsh.** Fish is installed only because caelestia-cli's theming scripts require it.
**Terminal stays Ghostty** with your custom shader config. The foot terminal config (caelestia's default) is skipped.

---

## Files added / changed in this repo

| File | What changed |
|---|---|
| `install_caelestia.sh` | Full Fedora port install script — builds everything from source |
| `.config/caelestia/shell.json` | caelestia-shell config — bar excluded from DP-3 (right monitor) |
| `.config/hypr/hyprland.conf` | Replaced waybar/mako with caelestia-shell; added QML_IMPORT_PATH env |
| `deploy.sh` | Added caelestia shell.json symlink; commented out waybar/mako symlinks |

---

## Install

```bash
cd ~/repos/dot-files/fedora
./deploy.sh            # symlink configs first
./install_caelestia.sh # build and install everything
```

Optional flags:
- `--vscode` — symlink VS Code / VSCodium caelestia settings
- `--spotify` — set up Spicetify theme (Spotify must be installed separately)

The script is **idempotent** — safe to re-run; already-completed steps are skipped.

After install, **restart Hyprland**. caelestia-shell autostarts via hyprland.conf.

---

## Dual-monitor setup

| Monitor | Connector | Position | Workspaces |
|---|---|---|---|
| ASUS VG27AC1A | DP-2 | Left (main) | 1–5 (ALT+1-5) |
| Gigabyte M27Q | DP-3 | Right | 6–10 (ALT+Q-T) |

The bar runs only on DP-2 (main monitor). DP-3 is excluded via `shell.json`:

```json
{
  "bar": {
    "excludedScreens": ["DP-3"]
  }
}
```

---

## Customisation

### shell.json — `~/.config/caelestia/shell.json`

Symlinked from `dot-files/fedora/.config/caelestia/shell.json`.

Key options:
- `profile` — path to profile picture (shown in lock screen / session menu)
- `wallpaperDir` — directory caelestia picks wallpapers from
- `bar.excludedScreens` — list of connector names where the bar should NOT appear
- `bar.showOnHover` — hide bar until mouse moves to the edge
- `bar.persistent` — always show bar (overrides showOnHover)

Full bar config with all defaults is at `/usr/local/etc/xdg/quickshell/caelestia/config/BarConfig.qml`.

### Wallpaper & color scheme

1. Place images in `~/Pictures/Wallpapers`
2. Use the caelestia launcher or run `caelestia scheme` to regenerate the Material You theme from the current wallpaper
3. Run `caelestia wallpaper set <file>` to switch wallpaper

### Launcher / lock / session

- `Super+Space` — toggle app launcher (`caelestia shell drawers toggle launcher`)
- `Super+Escape` — lock screen
- Power menu is in the bar (bottom-right)

---

## Architecture: why things live where they do

### Quickshell install paths on Fedora

Built with `-DCMAKE_INSTALL_PREFIX=/usr/local -DINSTALL_QMLDIR=/usr/local/lib64/qt6/qml`.

- Binary: `/usr/local/bin/qs`
- QML modules: `/usr/local/lib64/qt6/qml/`

Fedora's Qt6 installation uses `/usr/lib64/qt6/` (not `/usr/lib/qt6/`), so `QML_IMPORT_PATH` must be set to tell Qt where to look. This is set in `hyprland.conf`:

```ini
env = QML_IMPORT_PATH,/usr/local/lib64/qt6/qml
```

### caelestia-shell install paths on Fedora

Built with `-DCMAKE_INSTALL_PREFIX=/usr/local`.

- QML config: `/usr/local/etc/xdg/quickshell/caelestia/`
- Caelestia QML plugin: `/usr/local/usr/lib/qt6/qml/Caelestia` ← doubled path (see bugs below)

Two symlinks are created during install:
1. `~/.config/quickshell/caelestia` → `/usr/local/etc/xdg/quickshell/caelestia` (quickshell only searches `~/.config/quickshell/`)
2. `/usr/local/lib64/qt6/qml/Caelestia` → `/usr/local/usr/lib/qt6/qml/Caelestia` (fixes doubled path)

### libcava on Fedora

Arch has a `libcava` AUR package that compiles cava's `cavacore.c` as a shared library. Fedora's `cava` package does not ship `libcava.so` or a pkg-config file, so we build it manually:

```
libcava.so   → /usr/local/lib/
cavacore.h   → /usr/local/include/cava/
cava.pc      → /usr/lib64/pkgconfig/
ldconfig     → /etc/ld.so.conf.d/usr-local.conf  (registers /usr/local/lib)
```

---

## All problems encountered and how they were fixed

### DNF package issues

| Problem | Fix |
|---|---|
| `power-profiles-daemon` conflicts with `tuned-ppd` | Removed — Fedora 43 ships `tuned-ppd` which serves the same purpose |
| `starship` not in Fedora repos | Installed via `curl -sS https://starship.rs/install.sh \| sh -s -- --yes` |
| `eza` not in Fedora 43 repos | Installed via `cargo install eza` |
| `fuzzel-fish-completion` doesn't exist | Removed from package list |
| `adw-gtk3` wrong package name | Renamed to `adw-gtk3-theme` |
| `polkit-gnome` not in Fedora 43 (dropped upstream) | Removed — `hyprland.conf` already has fallback path `/usr/libexec/polkit-gnome-authentication-agent-1` |

### quickshell build (CMake) failures

| Error | Fix |
|---|---|
| `Qt6ShaderTools not found` | Added `qt6-qtshadertools-devel` |
| `Qt6CorePrivate / Qt6QuickPrivate not found` | Added `qt6-qtbase-private-devel`, `qt6-qtdeclarative-private-devel` |
| `qt6-qtwayland-private-devel` doesn't exist | Removed — private headers not split out for wayland on Fedora |
| `CLI11 not found` | Added `cli11-devel` |
| `libdrm not found` | Added `libdrm-devel` |
| `cpptrace not in Fedora repos` | Added `-DCRASH_HANDLER=OFF` to cmake flags |
| `wayland-protocols < 1.41` | Added `wayland-protocols-devel` |
| `gbm (mesa) not found` | Added `mesa-libgbm-devel` |
| `polkit-agent-1 not found` | Added `polkit-devel glib2-devel` |
| `security/_pam_types.h: No such file` | Added `pam-devel` |
| `xcb-util not found` | Added `xcb-util-devel` |

### caelestia-shell build (CMake) failures

| Error | Fix |
|---|---|
| `VERSION not set` (shallow clone has no git tags) | Added `-DVERSION=0.1.0` to cmake flags |
| `aubio not found` | Added `aubio-devel` |
| `cava pkg-config not found` | Built `libcava.so` manually (see libcava section above) |
| `libtoolize not found` during libcava build | Added `libtool autoconf automake` to DNF |
| `tee: /usr/local/lib/pkgconfig/cava.pc: No such file` | Used `/usr/lib64/pkgconfig/` instead (this path exists on Fedora) |
| `pkg-config still can't find cava` after writing .pc file | `/usr/local/lib/pkgconfig` was in the .pc but `/usr/lib64/pkgconfig` is where it was written — fixed path consistency |

### Runtime failures (after install)

| Error | Fix |
|---|---|
| Black screen after Hyprland restart — `caelestia shell -d` says "Could not find caelestia config directory" | quickshell searches `~/.config/quickshell/`, but files are in `/usr/local/etc/xdg/quickshell/caelestia/`. Fix: `ln -sfn /usr/local/etc/xdg/quickshell/caelestia ~/.config/quickshell/caelestia` |
| `module "Caelestia" is not installed` | caelestia-shell cmake stacks its hardcoded `/usr/lib/qt6/qml` on top of the `/usr/local` prefix → plugin lands in `/usr/local/usr/lib/qt6/qml/Caelestia`. Fix: `sudo ln -sfn /usr/local/usr/lib/qt6/qml/Caelestia /usr/local/lib64/qt6/qml/Caelestia` |
| Caelestia QML module still not found even with symlink | Qt doesn't search `/usr/local/lib64/qt6/qml` by default on Fedora. Fix: set `QML_IMPORT_PATH=/usr/local/lib64/qt6/qml` at runtime and in `hyprland.conf` |
| `libcava.so: cannot open shared object file: No such file or directory` | `/usr/local/lib` not registered with the dynamic linker. Fix: `echo '/usr/local/lib' \| sudo tee /etc/ld.so.conf.d/usr-local.conf && sudo ldconfig` |

---

## Rebuilding after upstream updates

To rebuild any component, delete its build sentinel and re-run the script:

```bash
# Rebuild quickshell
rm -rf ~/.local/src/quickshell/build
sudo rm -f /usr/local/bin/qs
./install_caelestia.sh

# Rebuild caelestia-shell
rm -rf ~/.local/src/caelestia-shell/build
./install_caelestia.sh

# Rebuild libcava
sudo rm -f /usr/local/lib/libcava.so
sudo rm -f /usr/lib64/pkgconfig/cava.pc
./install_caelestia.sh
```

---

## Useful commands

```bash
# Restart caelestia-shell without restarting Hyprland
pkill -f "qs.*caelestia" && caelestia shell -d

# Check caelestia-shell logs
journalctl --user -u caelestia-shell -f
# or if running directly:
# check output from `caelestia shell -d` in your terminal

# Regenerate Material You color theme from current wallpaper
caelestia scheme

# Change wallpaper
caelestia wallpaper set ~/Pictures/Wallpapers/myimage.jpg

# Toggle launcher (NOTE: caelestia toggle launcher is WRONG — that toggles a Hyprland special workspace)
caelestia shell drawers toggle launcher

# Other drawers: bar, osd, session, dashboard, utilities, sidebar
caelestia shell drawers toggle session

# Lock screen
caelestia lock
```
