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
| **quickshell** | QML shell framework that caelestia-shell runs on top of (`quickshell-git` package, copr `errornointernet/quickshell`) | – |
| **m3shapes** | QML module (Material 3 shapes) the shell imports; commit pinned by the shell's `flake.nix` | – |
| **caelestia dotfiles** | Configs for btop, fastfetch, fish, starship, micro, uwsm | – |

**Shell stays Zsh.** Fish is installed only because caelestia-cli's theming scripts require it.
**Terminal stays Ghostty** with your custom shader config. The foot terminal config (caelestia's default) is skipped.

---

## Files added / changed in this repo

| File | What changed |
|---|---|
| `scripts/install_caelestia.sh` | Full Fedora port install script — packages, fonts, libcava, then `update-caelestia` |
| `.local/bin/update-caelestia` | Builds shell/cli/m3shapes at pinned tags into `/usr/local`; `--check` compares with the latest releases |
| `.config/caelestia/shell.json` | caelestia-shell config — `excludedScreens` empty, bar on both monitors |
| `.config/hypr/hyprland.lua` + `hyprland/*.lua` | Lua config (Hyprland ≥ 0.55): caelestia-shell autostart, QML_IMPORT_PATH env |
| `deploy.sh` | Added caelestia shell.json symlink; commented out waybar/mako symlinks |

---

## Install

```bash
cd ~/repos/dot-files/fedora
./deploy.sh            # symlink configs first
./scripts/install_caelestia.sh # build and install everything
```

Optional flags:
- `--vscode` — symlink VS Code / VSCodium caelestia settings
- `--spotify` — set up Spicetify theme (Spotify must be installed separately)

The script is **idempotent** — safe to re-run; already-completed steps are skipped.

After install, **restart Hyprland**. caelestia-shell autostarts via `hyprland/execs.lua`.

---

## Dual-monitor setup

| Monitor | Position | Workspaces |
|---|---|---|
| ASUS VG27AC1A | Left (main) — Skydimo strip is on this one | 1–5 (ALT+1-5) |
| Gigabyte M27Q | Right | 6–10 (ALT+Q-T) |

Connector names are deliberately **not** listed here, because they are not stable —
they shift when the GPU changes or the DisplayPort cables are reseated. The ASUS
has been `DP-2` and `DP-1`; the Gigabyte has been `DP-3`, `DP-1`, and `DP-2`.
`hyprland.lua` pins monitors by EDID `desc:` and `hyperhdr-picker` resolves by
EDID model at runtime, so neither needs updating when this happens.

`excludedScreens` is the one place that still takes a raw connector name — the bar
currently runs on **both** monitors because it is empty in `shell.json`. If you
restrict it, look the name up first and expect to redo it after any cable change:

```bash
hyprctl monitors -j | jq -r '.[] | "\(.name)\t\(.model)"'
```

```json
{
  "bar": {
    "excludedScreens": ["<current connector of the monitor to skip>"]
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

### Quickshell

`quickshell-git` from copr `errornointernet/quickshell` (the shell requires git quickshell,
not the tagged release): binary `/usr/bin/qs`, QML modules under `/usr/lib64/qt6/qml/`.
It updates with `sudo dnf upgrade`.

### caelestia-shell / m3shapes install paths on Fedora

Built by `update-caelestia` with `-DCMAKE_INSTALL_PREFIX=/usr/local -DINSTALL_QMLDIR=lib64/qt6/qml
-DINSTALL_LIBDIR=lib64/caelestia`. Both projects default to `usr/lib/...` relative to the prefix,
which used to produce a doubled `/usr/local/usr/lib/...` path (and a symlink hack to cover it).

- QML config: `/usr/local/etc/xdg/quickshell/caelestia/`
- QML plugins: `/usr/local/lib64/qt6/qml/Caelestia`, `/usr/local/lib64/qt6/qml/M3Shapes`
- Version helper: `/usr/local/lib64/caelestia/version`

Qt and the shell only look under `/usr` by default, so `hypr/hyprland/env.lua` sets:

```lua
hl.env("QML_IMPORT_PATH", "/usr/local/lib64/qt6/qml")
hl.env("CAELESTIA_LIB_DIR", "/usr/local/lib64/caelestia")
```

quickshell only searches `~/.config/quickshell/`, so `update-caelestia` links
`~/.config/quickshell/caelestia` → `/usr/local/etc/xdg/quickshell/caelestia`.

### Pinned versions (held for Qt 6.11)

`update-caelestia` pins **shell v2.0.2** with the cli commit its `flake.lock` pins (`d1c8c8f`).
Newer shells don't load on Fedora 43's Qt 6.10: v2.0.3+ has `id: char` in
`modules/lock/center/InputField.qml`, which Qt 6.10's QML parser rejects as a reserved word,
and v2.2.0+ uses `DoubleSpinBox` (new in Qt 6.11). Fedora 44 ships Qt 6.11. When bumping, take
the cli rev from the new shell tag's `flake.lock`.

### Hyprland Lua config caveats

Hyprland ≥ 0.55 with a Lua config rejects `hyprctl keyword` ("keyword can't work with non-legacy
parsers"). caelestia v2.0.2 uses it in two places: the Caps/Num Lock refresh binds (now declared
in `hyprland/keybinds.lua`) and game mode's Hyprland tweaks (animations/blur/gaps off), which
currently don't apply.

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

> History from the original March 2026 port. quickshell is now the `quickshell-git` package
> and the shell is built by `update-caelestia`, so the build-failure tables below no longer
> describe the current setup; they're kept for when something similar resurfaces.

### DNF package issues

| Problem | Fix |
|---|---|
| `power-profiles-daemon` conflicts with `tuned-ppd` | Removed — Fedora 43 ships `tuned-ppd` which serves the same purpose |
| `starship` not in Fedora repos | Installed via `curl -sS https://starship.rs/install.sh \| sh -s -- --yes` |
| `eza` not in Fedora 43 repos | Installed via `cargo install eza` |
| `fuzzel-fish-completion` doesn't exist | Removed from package list |
| `adw-gtk3` wrong package name | Renamed to `adw-gtk3-theme` |
| `polkit-gnome` not in Fedora 43 (dropped upstream) | Removed — the autostart in `hyprland/execs.lua` already has fallback path `/usr/libexec/polkit-gnome-authentication-agent-1` |

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
| `module "Caelestia" is not installed` | caelestia-shell cmake stacks its hardcoded `/usr/lib/qt6/qml` on top of the `/usr/local` prefix → plugin lands in `/usr/local/usr/lib/qt6/qml/Caelestia`. Fix: `update-caelestia` passes `-DINSTALL_QMLDIR=lib64/qt6/qml -DINSTALL_LIBDIR=lib64/caelestia` (the old symlink must be gone before installing, or the install writes through it) |
| Caelestia QML module still not found even with symlink | Qt doesn't search `/usr/local/lib64/qt6/qml` by default on Fedora. Fix: set `QML_IMPORT_PATH=/usr/local/lib64/qt6/qml` at runtime and in `hyprland/env.lua` |
| `libcava.so: cannot open shared object file: No such file or directory` | `/usr/local/lib` not registered with the dynamic linker. Fix: `echo '/usr/local/lib' \| sudo tee /etc/ld.so.conf.d/usr-local.conf && sudo ldconfig` |

---

## Rebuilding after upstream updates

```bash
sudo dnf upgrade            # Hyprland (sdegler copr) and quickshell-git (errornointernet copr)
update-caelestia --check    # installed vs latest caelestia-shell / caelestia-cli releases
update-caelestia            # rebuild + reinstall at the pinned tags (asks for sudo), restarts the shell
```

`update-check` (daily timer) flags dnf updates and new caelestia releases. To upgrade the shell,
bump `SHELL_TAG`/`CLI_TAG` in `.local/bin/update-caelestia`, commit, and run it.

To rebuild libcava:

```bash
sudo rm -f /usr/local/lib/libcava.so
sudo rm -f /usr/lib64/pkgconfig/cava.pc
./scripts/install_caelestia.sh
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
