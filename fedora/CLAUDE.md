# System Configuration — Nic's Fedora Workstation

## OS & Environment
- **Distro:** Fedora 43
- **Desktop:** Hyprland on **Wayland**
- **GPU:** AMD RX 9070 XT 16 GB (RDNA4/gfx1201) — in-kernel `amdgpu`, Mesa/RADV, no driver install needed
- **Session type:** Wayland — do not suggest X11-only solutions
- **Package manager:** `dnf`

## Hardware
- **Monitors:**
  - ASUS VG27AC1A, left/main, 2560x1440 @ 170Hz — workspaces 1–5
  - Gigabyte M27Q, right, 2560x1440 @ 170Hz — workspaces 6–10
  - **Never hardcode a `DP-N` connector name anywhere.** They are not stable: the
    Gigabyte went `DP-3` → `DP-1` during the Aug 2026 GPU swap (silently dropping to
    60 Hz because its `monitor=DP-3` line stopped matching), then the two swapped
    again when the cables were reseated. `hyprland.conf` matches by EDID `desc:` and
    `hyperhdr-picker` resolves by EDID model at runtime. Run `hyprctl monitors` if
    you ever need the current name.
- **Resolution:** 2K (2560x1440)

## Desktop Environment

### Hyprland config
**File:** `.config/hypr/hyprland.conf`
- Layout: dwindle
- Mod key: `Super`
- Terminal: `ghostty`
- Launcher: `caelestia shell drawers toggle launcher`

### caelestia-shell
Replaces both Waybar (bar) and Mako (notifications). Also provides: OSD overlays, app launcher, lock screen, session menu.
- Autostarted via `exec-once = caelestia shell -d` in hyprland.conf
- Config: `~/.config/caelestia/shell.json` (symlinked from dotfiles)
- `shell.json` `excludedScreens` is empty — the bar runs on both monitors
- See `docs/caelestia.md` for full details on the caelestia setup

### Wallpaper
`swww-daemon` runs at startup. Use `caelestia wallpaper set <file>` to change wallpaper.
Wallpapers live in `~/Pictures/Wallpapers`.

### Clipboard
`cliphist` manages clipboard history. Access via `Super+;`.

## Keybindings (Hyprland)

### Core
| Key | Action |
|-----|--------|
| `Super+Return` | Open terminal (ghostty) |
| `Super+Space` | Toggle app launcher |
| `Super+C` | Kill active window |
| `Super+F` | Fullscreen |
| `Super+Shift+F` | Fullscreen (maximize, keep bar) |
| `Super+V` | Toggle floating |
| `Super+G` | Game mode (suppresses notifications) |
| `Super+N` | Notification center (caelestia sidebar) |
| `Super+Escape` / `Super+L` | Lock screen |
| `Super+Shift+Ctrl+Q` | Exit Hyprland |

### Focus & Window Management
| Key | Action |
|-----|--------|
| `Alt+H/J/K/L` (or arrows) | Move focus left/down/up/right |
| `Alt+Shift+H/J/K/L` | Move window left/down/up/right |
| `Super+Ctrl+Arrows` | Resize active window (40px steps) |
| `Super+,` / `Super+.` | Focus left/right monitor |
| `Super+Shift+,` / `Super+Shift+.` | Move workspace to left/right monitor |

### Workspaces
| Key | Action |
|-----|--------|
| `Alt+1-5` | Switch to workspace 1–5 (left monitor, ASUS) |
| `Alt+Q/W/E/R/T` | Switch to workspace 6–10 (right monitor, Gigabyte) |
| `Shift+Alt+1-5` | Move window to workspace 1–5 |
| `Shift+Alt+Q/W/E/R/T` | Move window to workspace 6–10 |
| `Super+Scroll` | Cycle through workspaces |

### Mouse
| Action | Binding |
|--------|---------|
| Move window | `Super+LMB drag` |
| Resize window | `Super+RMB drag` |

### Utilities
| Key | Action |
|-----|--------|
| `Ctrl+Shift+3` | Screenshot full screen → clipboard |
| `Ctrl+Shift+4` | Screenshot region → clipboard |
| `Super+;` | Clipboard history (cliphist + wofi) |
| Media keys | Volume (wpctl), playback (playerctl), brightness (brightnessctl) |

## Shell & Terminal
- **Shell:** zsh with oh-my-zsh + Powerlevel10k theme
- **Terminal:** Ghostty (custom shaders: balatro.glsl, crt.glsl)
- **Prompt:** Powerlevel10k (`~/.p10k.zsh`)
- **Notable tools:** `eza` (modern ls, via cargo), `btop`, `fastfetch`

## Dotfile Management
**Script:** `deploy.sh` — symlinks all configs into `$HOME`. Run as regular user.
**Install caelestia:** `scripts/install_caelestia.sh` — builds quickshell + caelestia-shell from source on Fedora.
**Other scripts:** in `scripts/` — `startup_script.sh` (fresh-install bootstrap), `gaming_setup.sh` (kept separate from main setup), `grub_setup.sh` (revert GRUB to defaults).
**Update check:** `.local/bin/update-check`, run by the `update-check.timer` systemd user unit (3 min after login + daily). It auto-updates flatpaks and sends one notification when dnf has updates (flagging Hyprland/Mesa/kernel/caelestia ones) or an enabled copr is dead/stale for this release. It **never installs dnf packages**: `sudo dnf upgrade` stays manual.
Log: `~/.local/state/update-check.log`. Run by hand with `update-check --dry-run` (prints instead of notifying, doesn't apply flatpak updates).
**Docs:** in `docs/` — `caelestia.md`, `TODO.md`, `lights_setup.md` (local only, gitignored).

## Notes for Claude Code
- This is a **Wayland** system — always prefer Wayland-compatible solutions
- **AMD + Wayland**: only `LIBVA_DRIVER_NAME=radeonsi` is set (in `hyprland.conf` and `.config/environment.d/gpu-wayland.conf`). Never reintroduce `__GLX_VENDOR_LIBRARY_NAME` or `NVD_BACKEND` — they are NVIDIA-only and break GLX/VA-API on this box
- Package manager is `dnf` — confirm package names exist in Fedora repos before recommending
- caelestia-shell is the bar/notification daemon — don't suggest waybar or mako as replacements
- `docs/caelestia.md` has detailed notes on the caelestia install, architecture, and known Fedora-specific issues
